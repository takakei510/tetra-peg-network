"""Read and validate stored network CSVs without loading trajectories."""
import csv
import gzip
import hashlib
from collections import Counter
from .stats import ratio

def locate(folder, name):
    options=[folder/name,folder/(name+'.gz')]
    found=[p for p in options if p.is_file()]
    if not found:raise FileNotFoundError(f'{folder}: {name}(.gz) がありません')
    if len(found)>1:raise ValueError(f'{folder}: 圧縮前後の両方があります。使用するコピーを一方にしてください')
    return found[0]


def iter_csv(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as f:
        yield from csv.DictReader(f)


def cluster_metrics(rows, n, bonds, cluster_count, largest):
    data=[];ids=set()
    for r in rows:
        cid,size,b=int(r['cluster_id']),int(r['molecule_count']),int(r['bond_count'])
        if cid<0 or cid in ids or size<1 or not (size-1<=b<=2*size):
            raise ValueError('Invalid or duplicate connected cluster row')
        if size==1 and b!=0:raise ValueError('Self-bonds are excluded in AB model')
        ids.add(cid);data.append((cid,size,b))
    if (len(data),sum(r[1] for r in data),sum(r[2] for r in data),max((r[1] for r in data),default=0)) != (cluster_count,n,bonds,largest):
        raise ValueError('Cluster counts disagree with saved sample N/B/K/N1')
    tied=sorted(cid for cid,s,b in data if s==largest)
    selected=tied[0]
    s2=sum(s*s for cid,s,b in data)
    count=len(data)
    return dict(cluster_number_mean=n/count, cluster_weighted_all=s2/n,
                cluster_weighted_excl_one=ratio(s2-largest*largest,n-largest),
                cluster_weighted_excl_all_tied=ratio(s2-len(tied)*largest*largest,n-len(tied)*largest),
                largest_tie_count=len(tied),selected_largest_id=selected,
                cluster_s1_excl_one=n-largest,cluster_s2_excl_one=s2-largest*largest),Counter(s for cid,s,b in data)


def count_bonds(rows, n, expected_bonds):
    pairs=Counter();endpoints=set();ids=set()
    for r in rows:
        a,b=int(r['molecule_a']),int(r['molecule_b'])
        ea,eb=int(r['endpoint_a']),int(r['endpoint_b'])
        bid=int(r['bond_id'])
        if not (0<=a<n and 0<=b<n and a!=b and 0<=ea<4*n and 0<=eb<4*n):
            raise ValueError('Invalid molecule or endpoint ID')
        if ea//4!=a or eb//4!=b or ea==eb or ea in endpoints or eb in endpoints or bid<0 or bid in ids:
            raise ValueError('Duplicate endpoint/bond or endpoint owner mismatch')
        if 'arm_a' in r and (int(r['arm_a'])!=ea%4 or int(r['arm_b'])!=eb%4):
            raise ValueError('Arm ID mismatch')
        if 'type_a' in r and {r['type_a'],r['type_b']}!={'A','B'}:
            raise ValueError('Not an A-B bond')
        if 'manhattan_distance' in r and int(r['manhattan_distance']) not in (1,2):
            raise ValueError('Invalid Manhattan distance')
        if all(k in r for k in ['x_a','y_a','z_a','x_b','y_b','z_b']):
            distance=sum(abs(int(r[c+'_a'])-int(r[c+'_b'])) for c in 'xyz')
            if distance not in (1,2) or distance!=int(r['manhattan_distance']):
                raise ValueError('Coordinate/distance mismatch')
        endpoints.update((ea,eb));ids.add(bid)
        pairs[tuple(sorted((a,b)))]+=1
    if len(ids)!=expected_bonds or any(k>4 for k in pairs.values()):
        raise ValueError('Bond count mismatch or more than four bonds per pair')
    hist=Counter(pairs.values())
    return [hist[k] for k in range(1,5)]


def filehash(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


