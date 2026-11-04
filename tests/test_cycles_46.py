import sys
from pathlib import Path
import unittest
import itertools
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
from analyze_cycles_46 import count_cycles, aggregate, summarize_case


def bonds(edges):
    used={};out=[]
    for i,(a,b) in enumerate(edges):
        ea=4*a+used.get(a,0);eb=4*b+used.get(b,0)
        used[a]=used.get(a,0)+1;used[b]=used.get(b,0)+1
        out.append(dict(bond_id=str(i),molecule_a=str(a),molecule_b=str(b),endpoint_a=str(ea),endpoint_b=str(eb)))
    return out


def brute(n,edges,k):
    edge_set={tuple(sorted(e)) for e in edges};cycles=set()
    for vertices in itertools.combinations(range(n),k):
        start=min(vertices)
        for tail in itertools.permutations([v for v in vertices if v!=start]):
            cycle=(start,)+tail
            if all(tuple(sorted((cycle[i],cycle[(i+1)%k]))) in edge_set for i in range(k)):
                cycles.add(min(cycle,(start,)+tuple(reversed(tail))))
    return len(cycles)


class Tests(unittest.TestCase):
    def test_tree_has_no_cycles(self):
        rows,_,_=count_cycles(5,bonds([(0,1),(1,2),(2,3),(3,4)]))
        self.assertEqual((rows[0]['C4_mol'],rows[0]['C6_mol']),(0,0))

    def test_square_parallel_multiplicity_and_endpoint_union(self):
        rows,_,_=count_cycles(4,bonds([(0,1),(1,2),(2,3),(3,0),(0,1)]))
        r=rows[0]
        self.assertEqual((r['C4_mol'],r['C4_edge'],r['B4'],r['H2']),(1,2,5,1))
        self.assertEqual((r['q4'],r['qH2']),(.625,.25))

    def test_shared_edge_not_double_counted_and_chorded_hexagon(self):
        e=[(0,1),(1,2),(2,3),(3,0),(1,4),(4,5),(5,0)]
        r=count_cycles(6,bonds(e))[0][0]
        self.assertEqual((r['C4_mol'],r['C6_mol'],r['B4'],r['B6']),(2,1,7,6))
        self.assertAlmostEqual(r['q4'],7/12)

    def test_k33_known_counts_and_brute_oracle(self):
        e=[(a,b) for a in range(3) for b in range(3,6)]
        r=count_cycles(6,bonds(e))[0][0]
        self.assertEqual((r['C4_mol'],r['C6_mol']),(9,6))
        for k in (4,6):self.assertEqual(r[f'C{k}_mol'],brute(6,e,k))

    def test_k44_counts(self):
        e=[(a,b) for a in range(4) for b in range(4,8)]
        r=count_cycles(8,bonds(e))[0][0]
        self.assertEqual((r['C4_mol'],r['C6_mol']),(36,96))
        self.assertEqual((r['q4'],r['q6']),(1,1))
        self.assertGreater(r['C6_per_molecule'],1)

    def test_vertex_relabeling_and_components(self):
        e=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4)]
        mapping={i:7-i for i in range(8)}
        rows,_,_=count_cycles(9,bonds([(mapping[a],mapping[b]) for a,b in reversed(e)]))
        self.assertEqual(sum(r['C4_mol'] for r in rows),2)
        self.assertEqual([r['cluster_id'] for r in rows if r['is_largest']],[0])

    def test_invalid_reuse_and_odd_cycle(self):
        bb=bonds([(0,1),(1,2)]);bb[1]['endpoint_a']=bb[0]['endpoint_b']
        with self.assertRaises(ValueError):count_cycles(3,bb)
        with self.assertRaises(ValueError):count_cycles(3,bonds([(0,1),(1,2),(2,0)]))

    def test_missing_bins_excluded_but_existing_zero_included(self):
        case=dict(case_id='a',L=128,target_phi=.43)
        r=count_cycles(5,bonds([(0,1),(1,2),(2,3),(3,0)]))[0]
        bins=summarize_case(case,r)
        # In all scope, a present isolated cluster has zero q4; it contributes.
        s=aggregate([{'bins':bins}],2)
        isolated=next(x for x in s if x['scope']=='all' and x['size_low']==1)
        self.assertEqual((isolated['q4_mean'],isolated['n_contributing'],isolated['cluster_count_mean']),(0,1,.5))

    def test_many_small_graphs_against_independent_brute_force(self):
        import random
        rng=random.Random(27)
        for _ in range(25):
            e=[(a,b) for a in range(3) for b in range(3,6) if rng.random()<.6]
            rows=count_cycles(6,bonds(e))[0]
            for k in (4,6):self.assertEqual(sum(r[f'C{k}_mol'] for r in rows),brute(6,e,k))

if __name__=='__main__':unittest.main()
