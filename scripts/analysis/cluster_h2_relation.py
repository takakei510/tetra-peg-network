"""Independent graph reconstruction and within-sample cluster/H2 measurements."""
from collections import Counter, defaultdict, deque
from h2_cluster.raw import iter_csv, locate, count_bonds
from h2_cluster.stats import ratio


def graph_rows(n, bonds):
    """BFS is independent of the simulation's Union-Find; IDs follow min molecule."""
    pairs = Counter(tuple(sorted((int(r['molecule_a']), int(r['molecule_b'])))) for r in bonds)
    adjacent = [set() for _ in range(n)]
    for a, b in pairs:
        adjacent[a].add(b); adjacent[b].add(a)
    owner = [-1]*n
    rows = []
    for start in range(n):
        if owner[start] >= 0:
            continue
        cid = len(rows); owner[start] = cid
        queue = deque([start]); members = []
        while queue:
            a = queue.popleft(); members.append(a)
            for b in adjacent[a]:
                if owner[b] < 0:
                    owner[b] = cid; queue.append(b)
        rows.append(dict(cluster_id=cid, size=len(members), bonds=0, h2=0))
    for (a, b), multiplicity in pairs.items():
        r = rows[owner[a]]
        r['bonds'] += multiplicity
        r['h2'] += multiplicity == 2  # exactly two, not >= two
    for r in rows:
        r['endpoint_fraction_all'] = r['h2']/r['size']
        r['endpoint_fraction_bonded'] = ratio(2*r['h2'], r['bonds'])
        r['p'] = r['bonds']/(2*r['size'])
    return rows, owner


def inspect(folder, sample):
    bonds = list(iter_csv(locate(folder, 'bonds.csv')))
    count_bonds(bonds, sample['N'], sample['bonds'])
    rows, owner = graph_rows(sample['N'], bonds)
    expected = {int(r['cluster_id']):(int(r['molecule_count']), int(r['bond_count']))
                for r in iter_csv(locate(folder, 'cluster-summary.csv'))}
    observed = {r['cluster_id']:(r['size'], r['bonds']) for r in rows}
    if observed != expected:
        raise ValueError('BFS disagrees with original cluster-summary.csv')
    seen = set(); degrees = Counter()
    for b in bonds:
        degrees[int(b['molecule_a'])] += 1; degrees[int(b['molecule_b'])] += 1
    for r in iter_csv(locate(folder, 'cluster-members.csv')):
        m = int(r['molecule'])
        if m in seen or not 0 <= m < len(owner):
            raise ValueError('Invalid/duplicate membership')
        seen.add(m)
        if (int(r['cluster_id']), int(r['cluster_size']), int(r['bond_degree'])) != (
                owner[m], rows[owner[m]]['size'], degrees[m]):
            raise ValueError('BFS disagrees with original membership/degree')
    if len(seen) != sample['N'] or sum(r['h2'] for r in rows) != sample['H2']:
        raise ValueError('Membership/H2 total mismatch')
    return rows


def size_bin(s):
    # 1, 2, 3..4, 5..8, ...; fixed bins, not chosen after viewing results.
    high = 1 << (s-1).bit_length()
    return (1 if high == 1 else high//2+1), high


def bin_rows(rows):
    groups = defaultdict(list)
    largest = min(rows, key=lambda r:(-r['size'], r['cluster_id']))['cluster_id']
    for scope in ('all', 'exclude_one_largest'):
        for r in rows:
            if scope == 'exclude_one_largest' and r['cluster_id'] == largest:
                continue
            groups[(scope, *size_bin(r['size']))].append(r)
    result = []
    for (scope, low, high), group in sorted(groups.items()):
        result.append(dict(scope=scope, size_low=low, size_high=high,
            cluster_count=len(group), molecule_count=sum(r['size'] for r in group),
            h2=sum(r['h2'] for r in group),
            fraction_mean=sum(r['endpoint_fraction_all'] for r in group)/len(group),
            fraction_pooled=sum(r['h2'] for r in group)/sum(r['size'] for r in group)))
    return result
