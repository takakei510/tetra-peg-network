"""Linux/WSLで独立プロセスを測る。外部Pythonパッケージは不要。"""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import threading
import time

MIB = 1024 ** 2


class BudgetStop(RuntimeError):
    """データを壊さず、以後の試料を開始しないための停止。"""


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(MIB), b''):
            result.update(block)
    return result.hexdigest()


def directory_bytes(path):
    size = 0
    for file in Path(path).rglob('*'):
        try:
            if file.is_file() and not file.is_symlink():
                size += file.stat().st_size
        except FileNotFoundError:
            pass  # atomic_jsonの一時ファイルがrenameされた場合。
    return size


def proc_rss(pid):
    try:
        for line in Path(f'/proc/{pid}/status').read_text().splitlines():
            if line.startswith('VmRSS:'):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return 0


def available_memory():
    """Linuxの空きメモリと、検出できたcgroup v2残量の小さい方。"""
    values = []
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'):
                values.append(int(line.split()[1]) * 1024)
    except (OSError, ValueError):
        pass
    try:
        maximum = Path('/sys/fs/cgroup/memory.max').read_text().strip()
        if maximum != 'max':
            current = int(Path('/sys/fs/cgroup/memory.current').read_text())
            values.append(max(0, int(maximum) - current))
    except (OSError, ValueError):
        pass
    return min(values) if values else None


def environment():
    result = {'platform': platform.platform(), 'machine': platform.machine(),
              'python': platform.python_version(), 'cpu_count': os.cpu_count(),
              'available_memory_bytes': available_memory(),
              'compiler_environment': {k: os.environ.get(k) for k in ('CC', 'CFLAGS', 'CPPFLAGS')},
              'measurement': 'Linux wait4 ru_maxrss, wall perf_counter, sequential processes',
              'monitor_interval_seconds': 0.02, 'warmup_runs': 0,
              'page_cache': 'not cleared', 'process_startup_included': True}
    for key, path in [('cpuinfo', '/proc/cpuinfo'), ('meminfo', '/proc/meminfo'),
                      ('cgroup_memory_max', '/sys/fs/cgroup/memory.max'),
                      ('cgroup_memory_current', '/sys/fs/cgroup/memory.current')]:
        try:
            value = Path(path).read_text()
            result[key] = next((line.split(':', 1)[1].strip() for line in value.splitlines()
                                if line.startswith('model name')), 'unknown') if key == 'cpuinfo' else value
        except OSError:
            result[key] = None
    try:
        result['matplotlib_version'] = importlib.metadata.version('matplotlib')
    except importlib.metadata.PackageNotFoundError:
        result['matplotlib_version'] = None
    return result


def disk_guard(config, suite, case=None, extra=0):
    if shutil.disk_usage(suite).free - extra < config['min_free_disk_mib'] * MIB:
        return 'free_disk_limit'
    if directory_bytes(suite) + extra > config['max_suite_output_mib'] * MIB:
        return 'suite_output_limit'
    if case and directory_bytes(case) + extra > config['max_case_output_mib'] * MIB:
        return 'case_output_limit'
    return None


def measure(command, stdout, stderr, config, suite, case):
    """wait4は主スレッドで直ちに回収し、短い実行を監視周期で量子化しない。"""
    reason = disk_guard(config, suite, case)
    if reason:
        return {'status': 'blocked', 'stop_reason': reason, 'command': command}
    stop = threading.Event()
    state = {'reason': None, 'sampled_peak': 0}
    with Path(stdout).open('wb') as output, Path(stderr).open('wb') as errors:
        started = time.perf_counter()
        process = subprocess.Popen(command, stdout=output, stderr=errors,
                                   stdin=subprocess.DEVNULL, start_new_session=True)

        def monitor():
            while not stop.is_set():
                rss = proc_rss(process.pid)
                state['sampled_peak'] = max(state['sampled_peak'], rss)
                reason = None
                if time.perf_counter() - started > config['stage_timeout_seconds']:
                    reason = 'time_limit'
                elif rss > config['max_rss_mib'] * MIB:
                    reason = 'rss_limit'
                else:
                    reason = disk_guard(config, suite, case)
                if reason:
                    state['reason'] = reason
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    return
                stop.wait(0.02)

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()
        interrupted = False
        try:
            _, status, usage = os.wait4(process.pid, 0)
        except KeyboardInterrupt:
            interrupted = True
            state['reason'] = 'interrupted'
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, status, usage = os.wait4(process.pid, 0)
        finally:
            elapsed = time.perf_counter() - started
            stop.set()
            thread.join()
        process.returncode = os.waitstatus_to_exitcode(status)
    # 終了直前の容量増加・短時間のRSS増加も判定する。
    peak = int(usage.ru_maxrss) * 1024
    reason = state['reason'] or disk_guard(config, suite, case)
    if peak > config['max_rss_mib'] * MIB:
        reason = reason or 'rss_limit_after_exit'
    return {'status': 'interrupted' if interrupted else 'limited' if reason else
                      'completed' if process.returncode == 0 else 'nonzero',
            'stop_reason': reason, 'exit_code': process.returncode,
            'wall_seconds': elapsed, 'cpu_user_seconds': usage.ru_utime,
            'cpu_system_seconds': usage.ru_stime, 'peak_rss_bytes': peak,
            'sampled_peak_rss_bytes': state['sampled_peak'], 'command': command}
