#!/usr/bin/env python3
"""保存済みC出力から表示専用JSONを作成する。元CSVは変更しない。"""
import argparse
import csv
import json
from pathlib import Path


def convert(trajectory, bonds=None):
    molecules = {}
    count = 0
    with Path(trajectory).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not {'molecule', 'type', 'arm', 'step', 'x', 'y', 'z'}.issubset(reader.fieldnames or []):
            raise ValueError('軌跡CSVの必要な列がありません')
        for row in reader:
            count += 1
            if count > 500000:
                raise ValueError('ビューアの軌跡50万行上限を超えています')
            identity, arm, step = (int(row[k]) for k in ('molecule', 'arm', 'step'))
            point = [int(row[k]) for k in ('x', 'y', 'z')]
            if identity < 0 or min(point) < 0 or row['type'] not in ('A', 'B'):
                raise ValueError('分子ID・座標・型が不正です')
            molecule = molecules.setdefault(identity, [identity, row['type'], None, None, None, [[], [], [], []]])
            if molecule[1] != row['type']:
                raise ValueError('同じ分子の型が一致しません')
            if arm == -1 and step == 0 and molecule[2] is None:
                molecule[2:5] = point
            elif 0 <= arm < 4 and step == len(molecule[5][arm]) // 3 + 1:
                molecule[5][arm].extend(point)
            else:
                raise ValueError('中心・腕・stepが不正です。C出力順のCSVを使用してください')
    if not molecules:
        raise ValueError('空試料は表示できません')
    length = None
    for molecule in molecules.values():
        if molecule[2] is None:
            raise ValueError('中心がありません')
        for arm in molecule[5]:
            if length is None:
                length = len(arm)
            if not length or len(arm) != length:
                raise ValueError('腕長が一致しません')
    edges = []
    used = set()
    if bonds:
        with Path(bonds).open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle)
            if not {'molecule_a', 'molecule_b', 'arm_a', 'arm_b'}.issubset(reader.fieldnames or []):
                raise ValueError('結合CSVの必要な列がありません')
            for row in reader:
                a, b, aa, ab = (int(row[k]) for k in ('molecule_a', 'molecule_b', 'arm_a', 'arm_b'))
                if a not in molecules or b not in molecules or a == b or not (0 <= aa < 4 and 0 <= ab < 4):
                    raise ValueError('結合先が不正です')
                if molecules[a][1] == molecules[b][1] or (a, aa) in used or (b, ab) in used:
                    raise ValueError('結合型・末端重複が不正です')
                # 元CSVの末端ID・座標も照合し、違う試料を混ぜない。
                for side, identity, arm in [('a', a, aa), ('b', b, ab)]:
                    if 'endpoint_' + side in row and int(row['endpoint_' + side]) != 4 * identity + arm:
                        raise ValueError('末端IDが一致しません')
                    point = molecules[identity][5][arm][-3:]
                    for k, axis in enumerate('xyz'):
                        if axis + '_' + side in row and int(row[axis + '_' + side]) != point[k]:
                            raise ValueError('軌跡と結合の座標が一致しません')
                used.update([(a, aa), (b, ab)])
                edges.append([a, b, aa, ab])
                if len(edges) > 200000:
                    raise ValueError('ビューアの結合20万行上限を超えています')
    return {'format': 'ab-saw-viewer-1', 'hasBonds': bonds is not None,
            'molecules': [molecules[k] for k in sorted(molecules)], 'bonds': edges}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trajectories', required=True, type=Path)
    parser.add_argument('--bonds', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    try:
        data = convert(args.trajectories, args.bonds)
        # xモードで既存成果物の上書きを防ぐ。出力途中の失敗は削除する。
        with args.output.open('x', encoding='utf-8') as handle:
            try:
                json.dump(data, handle, separators=(',', ':'), ensure_ascii=False)
            except BaseException:
                args.output.unlink(missing_ok=True)
                raise
        print(f'{args.output}: {args.output.stat().st_size} bytes')
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
