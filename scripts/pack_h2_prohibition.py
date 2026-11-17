#!/usr/bin/env python3
"""Package a completed no-multiple-bond study for transfer to Google Drive."""
import argparse
import csv
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
import shutil
from zipfile import ZipFile, ZIP_STORED

from benchmark_resources import atomic_json, digest
from run_h2_validation import read

RAW_FILES = (
    'bonds.csv.gz', 'cluster-members.csv.gz', 'cluster-summary.csv.gz',
    'network-summary.csv.gz', 'placement-summary.csv.gz', 'config.cfg',
    'metadata.log', 'run.json',
)


def pack(study, export, destination):
    manifest = read(export / 'manifest.json')
    if (manifest['status'] != 'completed' or manifest['completed_samples'] != 18300 or
            manifest['requested_samples'] != 18300 or
            manifest['study']['allow_multiple_bonds'] is not False or
            manifest['study'] != read(study / 'study.json')):
        raise ValueError('Expected a complete 18,300-sample export from this prohibition study')

    with (export / 'samples.csv').open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    cases = {row['case_id'] for row in rows}
    keys = Counter((int(row['L']), Decimal(row['target_phi'])) for row in rows)
    expected = {(L, Decimal(i) / 100): 100 for L in (32, 64, 128) for i in range(10, 71)}
    if (len(rows) != len(cases) or keys != expected or
            any(int(row['allow_multiple_bonds']) or any(int(row[k]) for k in ('H2', 'H3', 'H4'))
                for row in rows)):
        raise ValueError('Unexpected samples, counts, or multiple bonds in export')

    groups = defaultdict(list)
    for row in rows:
        group = (int(Decimal(row['target_phi']) * 100) - 10) // 10
        groups[int(row['L']), group].append(row['case_id'])

    destination.mkdir(parents=True, exist_ok=False)
    archives = []
    for (L, group), case_ids in sorted(groups.items()):
        name = f'L{L}_phi{10+10*group:02d}-{min(19+10*group,70):02d}.zip'
        path = destination / name
        temporary = path.with_suffix('.partial')
        with ZipFile(temporary, 'w', compression=ZIP_STORED) as archive:
            for case_id in sorted(case_ids):
                case = study / 'cases' / case_id
                run = case / 'runs' / case_id
                record = read(case / 'result.json')
                if record['identity'] != manifest['identity'] or record['case']['case_id'] != case_id:
                    raise ValueError(f'Run identity mismatch for {case_id}')
                inputs = [run / name for name in RAW_FILES]
                inputs += [case / 'clusters.csv.gz', case / 'result.json']
                for source in inputs:
                    relative = source.relative_to(case)
                    expected_hash = record['hashes'].get(str(relative))
                    if expected_hash is not None and digest(source) != expected_hash:
                        raise ValueError(f'Hash mismatch for {source}')
                    archive.write(source, f'cases/{case_id}/{relative.as_posix()}')
        with ZipFile(temporary) as archive:
            bad = archive.testzip()
            if bad is not None:
                raise ValueError(f'CRC mismatch for {bad}')
        temporary.replace(path)
        archives.append(dict(file=name, samples=len(case_ids), sha256=digest(path),
                             bytes=path.stat().st_size))
        print(f'{name}: {len(case_ids)} samples', flush=True)

    tables = ('samples.csv', 'summary.csv', 'manifest.json', 'convergence.csv',
              'size_relation.csv', 'storage_report.json')
    for name in tables:
        shutil.copy2(export / name, destination / name)
    for name in ('study.json', 'environment.json', 'source.zip'):
        shutil.copy2(study / name, destination / name)
    atomic_json(destination / 'drive_manifest.json',
                dict(samples=len(rows), archives=archives, study_identity=manifest['identity'],
                     tables={name: digest(destination / name) for name in tables}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('export', type=Path)
    parser.add_argument('destination', type=Path, help='New local folder')
    args = parser.parse_args()
    pack(args.study, args.export, args.destination)


if __name__ == '__main__':
    main()
