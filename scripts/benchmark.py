#!/usr/bin/env python3
"""設定駆動でL別の時間・RSS・出力容量を測定する（Linux/WSL）。"""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
import uuid

from benchmark_resources import (MIB, BudgetStop, atomic_json, available_memory,
                                 digest, directory_bytes, disk_guard, environment, measure)
from run_demo import git_output, snapshot_sources

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'configs/benchmarks/pc-pilot.json'
SETTINGS = {'schema_version', 'sizes', 'target_phis', 'seeds', 'repeats', 'arm_length',
            'a_probability', 'attempts_per_molecule', 'min_attempts', 'postprocess',
            'postprocess_max_molecules', 'plot', 'plot_max_molecules', 'plot_dpi',
            'stage_timeout_seconds', 'max_rss_mib', 'max_case_output_mib',
            'max_suite_output_mib', 'min_free_disk_mib', 'retention', 'compression_level'}
PHASES = ('c', 'postprocess', 'plot')
BASE_FIELDS = ['case_id', 'L', 'target_phi', 'seed', 'repeat', 'target_molecules',
               'placed_molecules', 'actual_phi', 'attempts', 'center_occupied', 'trapped']
STAGE_FIELDS = BASE_FIELDS + ['stage', 'status', 'stop_reason', 'exit_code', 'wall_seconds',
                             'cpu_user_seconds', 'cpu_system_seconds', 'peak_rss_bytes',
                             'sampled_peak_rss_bytes', 'read_validate_seconds', 'render_save_seconds']
CASE_FIELDS = BASE_FIELDS + ['status', 'generation_status', 'raw_csv_bytes', 'derived_csv_bytes',
                            'gzip_csv_bytes', 'gzip_derived_csv_bytes', 'compression_wall_seconds',
                            'saved_bytes', 'largest', 'largest_fraction', 'directory']


def no_duplicate_keys(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f'JSONキーの重複: {key}')
        result[key] = value
    return result


def load_config(path):
    config = json.loads(Path(path).read_text(), object_pairs_hook=no_duplicate_keys)
    if set(config) != SETTINGS:
        raise ValueError(f'設定キーの過不足: missing={SETTINGS-set(config)}, unknown={set(config)-SETTINGS}')
    if type(config['schema_version']) is not int or config['schema_version'] != 1:
        raise ValueError('schema_versionは1のみ対応')
    for key in ('sizes', 'seeds', 'target_phis'):
        if not isinstance(config[key], list) or not config[key]:
            raise ValueError(f'{key}は空でない配列にしてください')
        if len(set(map(str, config[key]))) != len(config[key]):
            raise ValueError(f'{key}に重複があります')
    if any(type(v) is not int or not 2 <= v <= 1290 for v in config['sizes']):
        raise ValueError('sizesは2〜1290の整数')
    if config['sizes'] != sorted(config['sizes']):
        raise ValueError('sizesは小さい順にしてください')
    if any(type(v) is not int or not 0 <= v < 2**64 for v in config['seeds']):
        raise ValueError('seedは符号なし64ビット整数')
    for value in config['target_phis'] + [config['a_probability']]:
        number = Decimal(str(value))
        if not number.is_finite() or not 0 <= number <= 1:
            raise ValueError('占有率とA型確率は有限な0〜1の数値')
    if len({Decimal(str(v)) for v in config['target_phis']}) != len(config['target_phis']):
        raise ValueError('同じ値の目標占有率が重複しています')
    if any(Decimal(str(v)) == 0 for v in config['target_phis']):
        raise ValueError('目標占有率は正にしてください')
    for key in SETTINGS - {'schema_version', 'sizes', 'seeds', 'target_phis', 'a_probability',
                           'postprocess', 'plot', 'retention', 'compression_level'}:
        if type(config[key]) not in (int, float) or config[key] <= 0 or not float(config[key]) < float('inf'):
            raise ValueError(f'{key}は正の有限数')
    for key in ('repeats', 'arm_length', 'attempts_per_molecule', 'min_attempts',
                'postprocess_max_molecules', 'plot_max_molecules', 'plot_dpi'):
        if type(config[key]) is not int:
            raise ValueError(f'{key}は整数')
    if any(type(config[k]) is not bool for k in ('postprocess', 'plot')):
        raise ValueError('postprocess/plotは真偽値')
    if config['plot'] and not config['postprocess']:
        raise ValueError('plotにはpostprocessが必要です')
    if config['retention'] not in ('gzip', 'raw'):
        raise ValueError('retentionはgzipまたはraw')
    if type(config['compression_level']) is not int or not 0 <= config['compression_level'] <= 9:
        raise ValueError('compression_levelは0〜9')
    return config


def append_csv(path, fields, row):
    exists = path.exists()
    with path.open('a', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def first_row(path):
    with path.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f'1行の集計を想定: {path}')
    return rows[0]


def gzip_file(path, level, retention, guard):
    """生成したgzipを復号して検証してから、明示方針に従い片方を保持。"""
    compressed = path.with_name(path.name + '.gz')
    temporary = path.with_name(path.name + '.gz.tmp')
    if compressed.exists() or temporary.exists():
        raise ValueError(f'圧縮先が既にあります: {path}')
    original_size, original_hash = path.stat().st_size, digest(path)
    try:
        # 最悪ほぼ元サイズの一時領域が必要。圧縮率を過信して開始しない。
        guard(original_size + 1024)
        with path.open('rb') as source, temporary.open('xb') as target:
            with gzip.GzipFile(filename='', mode='wb', fileobj=target, compresslevel=level, mtime=0) as archive:
                for block in iter(lambda: source.read(MIB), b''):
                    archive.write(block)
                    guard(0)
        actual = hashlib.sha256()
        with gzip.open(temporary, 'rb') as stream:
            for block in iter(lambda: stream.read(MIB), b''):
                actual.update(block)
        if actual.hexdigest() != original_hash:
            raise ValueError(f'gzip復号ハッシュ不一致: {path}')
        size = temporary.stat().st_size
        if retention == 'gzip':
            temporary.replace(compressed)
            path.unlink()  # このベンチマークが作った当該CSVだけ。照合済みで復元可能。
        else:
            temporary.unlink()
        return {'raw_bytes': original_size, 'gzip_bytes': size, 'sha256': original_hash,
                'retained': retention, 'compressed_file': compressed.name}
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def restore(run):
    result = json.loads((run / 'performance.json').read_text())
    files = result.get('csv_files', [])
    required = sum(f['raw_bytes'] for f in files if not (run / f['path']).exists())
    if shutil.disk_usage(run).free < required + 256*MIB:
        raise ValueError('復元に必要な空き容量がありません（256 MiBの余裕も必要）')
    for entry in files:
        path = run / entry['path']
        if not path.resolve().is_relative_to(run.resolve()) or path.is_symlink():
            raise ValueError('復元先が実行フォルダ外です')
        if path.exists():
            if digest(path) != entry['sha256']:
                raise ValueError(f'既存ファイルが保存時のハッシュと違います: {path}')
            continue
        archive = path.with_name(path.name + '.gz')
        if archive.is_symlink():
            raise ValueError('圧縮ファイルのシンボリックリンクは使用できません')
        with archive.open('rb') as source, path.open('xb') as target:
            try:
                with gzip.GzipFile(fileobj=source) as stream:
                    written = 0
                    for block in iter(lambda: stream.read(MIB), b''):
                        written += len(block)
                        if written > entry['raw_bytes']:
                            raise ValueError('復元サイズが記録された容量を超えています')
                        target.write(block)
            except BaseException:
                path.unlink(missing_ok=True)
                raise
        if path.stat().st_size != entry['raw_bytes'] or digest(path) != entry['sha256']:
            path.unlink()
            raise ValueError(f'復元ハッシュ不一致: {archive}')
    print(f'復元・照合済み: {run}')


def skipped(reason):
    return {'status': 'skipped', 'stop_reason': reason}


def _execute_case(config, suite, case_id, L, phi, seed, repeat, storage_policy="full"):
    run = suite / 'runs' / case_id
    run.mkdir(parents=True, exist_ok=False)
    cfg = run / 'config.cfg'
    cfg.write_text(f'dim=3\nL={L}\nseed={seed}\ntarget_phi={phi}\n'
                   f'arm_length={config["arm_length"]}\na_probability={config["a_probability"]}\n')
    result = {'case_id': case_id, 'L': L, 'target_phi': str(phi), 'seed': seed, 'repeat': repeat,
              'directory': str(run.relative_to(suite)), 'status': 'running', 'stages': {},
              'generation_status': 'not_started', 'csv_files': [], 'storage_policy': storage_policy}
    atomic_json(run / 'performance.json', result)
    checked = subprocess.run([str(ROOT / 'build/check_network_config'), str(cfg)],
                             capture_output=True, text=True, check=False)
    (run / 'config-check.log').write_text(checked.stdout + checked.stderr)
    if checked.returncode:
        raise ValueError('C設定検査に失敗: ' + checked.stderr.strip())
    conditions = dict(token.split('=', 1) for token in checked.stdout.split() if '=' in token)
    n = int(conditions['target_molecules'])
    result['target_molecules'] = n
    result['base_memory_estimate_bytes'] = 12*L**3 + n*(16 + 36*config['arm_length'])
    free = available_memory()
    if (result['base_memory_estimate_bytes'] > config['max_rss_mib']*MIB or
            (free is not None and result['base_memory_estimate_bytes'] + 256*MIB > free)):
        result['stages']['c'] = {'status': 'blocked', 'stop_reason': 'base_memory_preflight'}
        result['status'] = 'limited'
        result['stop_reason'] = 'base_memory_preflight'
    else:
        attempts = max(config['min_attempts'], config['attempts_per_molecule'] * n)
        result['max_placement_attempts'] = attempts
        command = [str(ROOT / 'build/demo_network'), str(cfg)] + [str(Path(os.devnull) if storage_policy == 'compact' and name == 'candidates.csv' else run / name) for name in
                   ('candidates.csv', 'bonds.csv', 'cluster-members.csv', 'cluster-summary.csv',
                    'network-summary.csv', 'placement-summary.csv')] + [str(attempts)]
        trajectory = Path(os.devnull) if storage_policy == 'compact' else run / 'trajectories.csv'
        stage = measure(command, trajectory, run / 'metadata.log', config, suite, run)
        result['stages']['c'] = stage
        if stage.get('stop_reason'):
            result['stop_reason'] = stage['stop_reason']
        if stage.get('exit_code') in (0, 3) and stage['status'] not in ('limited', 'interrupted', 'blocked'):
            placement = first_row(run / 'placement-summary.csv')
            result['placement'] = placement
            for key in ('placed_molecules', 'actual_phi', 'attempts', 'center_occupied', 'trapped'):
                result[key] = placement[key]
            completed = stage['exit_code'] == 0 and placement['completed'] == '1'
            result['generation_status'] = 'completed' if completed else 'incomplete'
            result['status'] = result['generation_status']
            stage['status'] = result['generation_status']
            if completed:
                net = first_row(run / 'network-summary.csv')
                result.update({key: net[key] for key in ('largest', 'largest_fraction')})
        else:
            result['generation_status'] = 'limited' if stage.get('stop_reason') else 'failed'
            result['status'] = 'interrupted' if stage['status'] == 'interrupted' else result['generation_status']
    completed = result['generation_status'] == 'completed'
    # analyze_graph用の最小互換情報。性能データの正本はperformance.json。
    atomic_json(run / 'run.json', {'schema_version': 2, 'mode': 'network',
                'status': 'completed' if completed else result['status'],
                'generation_completed': completed, 'analysis_completed': completed,
                'reported_conditions': {'L': str(L)}, 'entrypoint': 'benchmark.py'})
    for name in ('postprocess', 'plot'):
        enabled = config[name]
        limit = config[name + '_max_molecules']
        if not completed:
            stage = skipped('generation_not_completed')
        elif not enabled:
            stage = skipped('disabled')
        elif n > limit:
            stage = skipped('molecule_limit')
        elif name == 'plot' and result['stages']['postprocess']['status'] != 'completed':
            stage = skipped('postprocess_not_completed')
        elif result['status'] in ('limited', 'interrupted', 'failed'):
            stage = skipped('previous_stage_failed')
        else:
            command = [sys.executable, str(ROOT / 'scripts/benchmark_stage.py'), name, str(run),
                       '--dpi', str(config['plot_dpi'])]
            stage = measure(command, run / f'{name}.stdout.log', run / f'{name}.stderr.log',
                            config, suite, run)
            if stage['status'] != 'completed':
                if stage.get('stop_reason'):
                    result['stop_reason'] = stage['stop_reason']
                result['status'] = 'limited' if stage.get('stop_reason') else 'failed'
                if stage['status'] == 'interrupted':
                    result['status'] = 'interrupted'
            elif name == 'plot':
                stage.update(json.loads((run / 'plot-timing.json').read_text()))
        result['stages'][name] = stage
        atomic_json(run / 'performance.json', result)
    csv_paths = sorted(run.rglob('*.csv'))
    result['raw_csv_bytes'] = sum(p.stat().st_size for p in csv_paths if p.parent == run)
    result['derived_csv_bytes'] = sum(p.stat().st_size for p in csv_paths if p.parent != run)
    result['gzip_csv_bytes'] = result['gzip_derived_csv_bytes'] = 0
    start = time.perf_counter()
    if result['status'] not in ('limited', 'interrupted', 'failed'):
        def guard(extra):
            reason = disk_guard(config, suite, run, extra)
            if reason:
                raise BudgetStop(reason)
        try:
            for path in csv_paths:
                entry = gzip_file(path, config['compression_level'], config['retention'], guard)
                entry['path'] = path.relative_to(run).as_posix()
                result['csv_files'].append(entry)
                key = 'gzip_csv_bytes' if path.parent == run else 'gzip_derived_csv_bytes'
                result[key] += entry['gzip_bytes']
                atomic_json(run / 'performance.json', result)
        except BudgetStop as error:
            result['status'] = 'limited'
            result['stop_reason'] = str(error)
    result['compression_wall_seconds'] = time.perf_counter() - start
    finalize_size(run / 'performance.json', result, run)
    return result



def finalize_size(path, result, directory):
    # 保存量にこのJSON自身も含め、桁数変化による自己参照を収束させる。
    for _ in range(8):
        atomic_json(path, result)
        size = directory_bytes(directory)
        if result.get('saved_bytes') == size:
            return
        result['saved_bytes'] = size
    atomic_json(path, result)


def run_case(config, suite, case_id, L, phi, seed, repeat, storage_policy="full"):
    run = suite / 'runs' / case_id
    try:
        if storage_policy not in ("full", "compact"):
            raise ValueError("Unknown storage_policy")
        return _execute_case(config, suite, case_id, L, phi, seed, repeat, storage_policy)
    except (Exception, KeyboardInterrupt) as error:
        path = run / 'performance.json'
        result = json.loads(path.read_text()) if path.exists() else {
            'case_id': case_id, 'L': L, 'target_phi': str(phi), 'seed': seed,
            'repeat': repeat, 'directory': str(run.relative_to(suite)), 'stages': {}}
        result['status'] = ('interrupted' if isinstance(error, KeyboardInterrupt) else
                            'limited' if isinstance(error, BudgetStop) else 'failed')
        result['stop_reason'] = str(error) or type(error).__name__
        result['stages'].setdefault('c', {'status': 'failed', 'stop_reason': result['stop_reason']})
        run.mkdir(parents=True, exist_ok=True)
        finalize_size(path, result, run)
        return result


def summarize(suite, results):
    groups = {}
    for case in results:
        for name, stage in case['stages'].items():
            key = (case['L'], case['target_phi'], name)
            groups.setdefault(key, []).append(stage)
    fields = ['L', 'target_phi', 'stage', 'attempted', 'completed', 'incomplete', 'skipped', 'other']
    metrics = ['wall_seconds', 'peak_rss_bytes', 'render_save_seconds']
    fields += [f'{metric}_{stat}' for metric in metrics for stat in ('median', 'min', 'max')]
    with (suite / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for (L, phi, name), stages in sorted(groups.items()):
            good = [stage for stage in stages if stage['status'] == 'completed']
            row = {'L': L, 'target_phi': phi, 'stage': name, 'attempted': len(stages),
                   'completed': len(good), 'incomplete': sum(s['status'] == 'incomplete' for s in stages),
                   'skipped': sum(s['status'] == 'skipped' for s in stages)}
            row['other'] = len(stages) - row['completed'] - row['incomplete'] - row['skipped']
            for metric in metrics:
                values = [s[metric] for s in good if metric in s]
                if values:
                    row.update({f'{metric}_median': statistics.median(values),
                                f'{metric}_min': min(values), f'{metric}_max': max(values)})
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path, help='新規フォルダのみ。省略時はdata/benchmarks/<実行ID>')
    parser.add_argument('--restore', type=Path, help='指定試料のgzipをハッシュ照合しCSVへ復元')
    args = parser.parse_args()
    if args.restore:
        try:
            restore(args.restore.resolve())
            return 0
        except (OSError, ValueError, KeyError) as error:
            parser.exit(1, f'{error}\n')
    if sys.platform != 'linux' or not hasattr(os, 'wait4'):
        parser.error('ベンチマークはLinux/WSLで実行してください')
    try:
        config = load_config(args.config)
    except (OSError, ValueError, TypeError, InvalidOperation) as error:
        parser.error(str(error))
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8]
    suite = (args.output or ROOT / 'data/benchmarks' / run_id).resolve()
    suite.mkdir(parents=True, exist_ok=False)
    atomic_json(suite / 'config.json', config)
    atomic_json(suite / 'environment.json', environment())
    snapshot_sources(suite / 'source.zip')
    manifest = {'schema_version': 1, 'status': 'running', 'run_id': run_id,
                'git_commit': git_output('rev-parse', 'HEAD'), 'git_status': git_output('status', '--porcelain'),
                'started_at_utc': datetime.now(timezone.utc).isoformat(), 'cases': []}
    atomic_json(suite / 'benchmark.json', manifest)
    results, exit_code = [], 0
    try:
        reason = disk_guard(config, suite)
        if reason:
            raise BudgetStop(reason)
        with (suite / 'build.log').open('w') as log:
            subprocess.run(['make', '-B', 'build/demo_network', 'build/check_network_config'],
                           cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
        manifest['binary_sha256'] = digest(ROOT / 'build/demo_network')
        manifest['source_zip_sha256'] = digest(suite / 'source.zip')
        compiler = os.environ.get('CC', 'cc')
        try:
            import shlex
            manifest['compiler_version'] = subprocess.run(shlex.split(compiler) + ['--version'],
                capture_output=True, text=True, check=True).stdout.splitlines()[0]
        except (OSError, subprocess.SubprocessError, IndexError):
            manifest['compiler_version'] = None
        total = len(config['sizes'])*len(config['target_phis'])*len(config['seeds'])*config['repeats']
        print(f'保存先: {suite}\n予定: {total}試料（逐次実行）', flush=True)
        for L in config['sizes']:
            for phi in config['target_phis']:
                for seed in config['seeds']:
                    for repeat in range(1, config['repeats'] + 1):
                        case_id = f'case{len(results)+1:05d}_L{L}_seed{seed}_r{repeat}'
                        result = run_case(config, suite, case_id, L, phi, seed, repeat)
                        results.append(result)
                        append_csv(suite / 'cases.csv', CASE_FIELDS, result)
                        for name, stage in result['stages'].items():
                            append_csv(suite / 'stages.csv', STAGE_FIELDS, {**result, **stage, 'stage': name})
                        for entry in result['csv_files']:
                            append_csv(suite / 'files.csv', ['case_id', 'path', 'raw_bytes', 'gzip_bytes',
                                      'sha256', 'retained'], {'case_id': case_id, **entry})
                        manifest['cases'].append({'case_id': case_id, 'status': result['status']})
                        atomic_json(suite / 'benchmark.json', manifest)
                        ctime = result['stages']['c'].get('wall_seconds')
                        print(f'[{len(results)}/{total}] L={L} phi={phi} seed={seed} r={repeat} '
                              f'status={result["status"]} C={ctime}', flush=True)
                        if result['status'] in ('limited', 'interrupted', 'failed'):
                            raise BudgetStop(result.get('stop_reason') or result['status'])
        manifest['status'] = 'completed'
    except (BudgetStop, OSError, ValueError, subprocess.SubprocessError, KeyboardInterrupt) as error:
        manifest['status'] = 'stopped'
        manifest['stop_reason'] = str(error) or 'interrupted'
        print(f'停止: {manifest["stop_reason"]}', file=sys.stderr)
        exit_code = 1
    finally:
        summarize(suite, results)
        manifest['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        finalize_size(suite / 'benchmark.json', manifest, suite)
        print(f'状態: {manifest["status"]}\n集計: {suite / "summary.csv"}', flush=True)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())

