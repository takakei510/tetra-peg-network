"""多重枝・孤立・同率最大・サイクル・壊れた所属を検証する。描画依存なし。"""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/visualization"))
from analyze_graph import validate_graph


class GraphTests(unittest.TestCase):
    def test_triple_bond_and_isolate(self):
        bonds = [{"molecule_a": 0, "molecule_b": 1} for _ in range(3)]
        members = [dict(molecule=0, cluster_id=0, cluster_size=2, bond_degree=3),
                   dict(molecule=1, cluster_id=0, cluster_size=2, bond_degree=3),
                   dict(molecule=2, cluster_id=1, cluster_size=1, bond_degree=0)]
        _, groups, tied = validate_graph(range(3), bonds, members)
        self.assertEqual(groups[0], {0, 1})
        self.assertEqual(tied, [0])
        members[0]["bond_degree"] = 1
        with self.assertRaisesRegex(ValueError, "bond_degree"):
            validate_graph(range(3), bonds, members)

    def test_ties_and_isolates(self):
        members = [dict(molecule=m, cluster_id=m, cluster_size=1, bond_degree=0)
                   for m in range(3)]
        self.assertEqual(validate_graph(range(3), [], members)[2], [0, 1, 2])

    def test_cycle(self):
        bonds = [{"molecule_a": u, "molecule_b": v}
                 for u, v in [(0, 1), (1, 2), (2, 3), (3, 0)]]
        members = [dict(molecule=m, cluster_id=0, cluster_size=4, bond_degree=2)
                   for m in range(4)]
        self.assertEqual(validate_graph(range(4), bonds, members)[1][0], set(range(4)))

    def test_disconnected_same_label_rejected(self):
        members = [dict(molecule=m, cluster_id=0, cluster_size=2, bond_degree=0)
                   for m in range(2)]
        with self.assertRaisesRegex(ValueError, "BFS"):
            validate_graph(range(2), [], members)

    def test_cross_cluster_bond_rejected(self):
        members = [dict(molecule=m, cluster_id=m, cluster_size=1, bond_degree=1)
                   for m in range(2)]
        with self.assertRaisesRegex(ValueError, "異なるクラスタ"):
            validate_graph(range(2), [{"molecule_a": 0, "molecule_b": 1}], members)


if __name__ == "__main__":
    unittest.main()
