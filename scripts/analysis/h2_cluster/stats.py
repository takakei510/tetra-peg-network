"""CSV export and per-sample statistical summaries."""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

def write_csv(path, rows):
    if not rows:
        raise ValueError("No rows to export")
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def ratio(a, b):
    return a / b if b else None


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row['L'], row['target_phi'])].append(row)
    skip = {'L', 'target_phi', 'seed', 'case_id', 'cluster_source', 'bond_source', 'selected_largest_id'}
    metrics = [k for k in rows[0] if k not in skip]
    result = []
    for (size, phi), group in sorted(groups.items()):
        r = dict(L=size, target_phi=phi, n_samples=len(group))
        for k in metrics:
            values = [x[k] for x in group if isinstance(x.get(k), (int, float))]
            if not values:
                continue
            r.update({k+'_n': len(values), k+'_mean': statistics.mean(values),
                      k+'_sd': statistics.stdev(values) if len(values)>1 else None,
                      k+'_min': min(values), k+'_max': max(values)})
        result.append(r)
    return result


