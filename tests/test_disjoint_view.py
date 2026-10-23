from copy import deepcopy
from genome_kit._disjoint_view import DisjointView
import numpy as np
import unittest


class TestDisjointView(unittest.TestCase):

    def setUp(self):
        self.parent = np.arange(40).reshape(10, 4)
        self.disjoint_view = DisjointView(self.parent, [self.parent[0:2], self.parent[5:8]], axis=0)
        self.expected = np.concatenate([self.parent[0:2], self.parent[5:8]])
        self.untouched = np.ones(10, dtype=bool)
        self.untouched[0:2] = False
        self.untouched[5:8] = False

    def test_read(self):
        np.testing.assert_array_equal(self.expected, self.disjoint_view[...])
        np.testing.assert_array_equal(self.expected[3], self.disjoint_view[3])
        np.testing.assert_array_equal(self.expected, np.asarray(self.disjoint_view))
        self.assertEqual(self.expected.shape, self.disjoint_view.shape)
        self.assertEqual(self.expected.dtype, self.disjoint_view.dtype)
        self.assertEqual(len(self.expected), len(self.disjoint_view))
        self.assertEqual(self.expected.sum(), self.disjoint_view.sum())
        np.testing.assert_array_equal(self.expected, list(self.disjoint_view))

    def test_operators(self):
        np.testing.assert_array_equal(self.expected + 1, self.disjoint_view + 1)
        np.testing.assert_array_equal(self.expected == 4, self.disjoint_view == 4)
        np.testing.assert_array_equal(np.add.reduce(self.expected), np.add.reduce(self.disjoint_view))

    def test_set_all(self):
        before = self.parent.copy()
        self.disjoint_view[...] = -1
        np.testing.assert_array_equal(-1, self.parent[0:2])
        np.testing.assert_array_equal(-1, self.parent[5:8])
        np.testing.assert_array_equal(before[self.untouched], self.parent[self.untouched])

    def test_set_index(self):
        before = self.parent.copy()
        self.disjoint_view[2] = -1
        np.testing.assert_array_equal(-1, self.parent[5])
        before[5] = -1
        np.testing.assert_array_equal(before, self.parent)

    def test_set_slice_across_regions(self):
        before = self.parent.copy()
        self.disjoint_view[1:3] = np.full((2, 4), -1)
        before[1] = -1
        before[5] = -1
        np.testing.assert_array_equal(before, self.parent)

    def test_in_place_operator(self):
        before = self.parent.copy()
        view = self.disjoint_view
        view += 1
        self.assertIs(self.disjoint_view, view)
        before[0:2] += 1
        before[5:8] += 1
        np.testing.assert_array_equal(before, self.parent)

    def test_ufunc_out(self):
        before = self.parent.copy()
        np.negative(self.disjoint_view, out=self.disjoint_view)
        before[0:2] *= -1
        before[5:8] *= -1
        np.testing.assert_array_equal(before, self.parent)

    def test_reflects_parent_changes(self):
        self.parent[6] = -1
        np.testing.assert_array_equal(-1, self.disjoint_view[3])

    def test_reads_are_read_only(self):
        with self.assertRaises(ValueError):
            self.disjoint_view.concatenated_data[...] = 1
        with self.assertRaises(ValueError):
            np.asarray(self.disjoint_view)[...] = 1
        with self.assertRaises(ValueError):
            self.disjoint_view.view()[...] = 1

    def test_slice_returns_disjoint_view(self):
        keys = [
            np.s_[...], np.s_[1:4], np.s_[::-1], np.s_[4:0:-2], np.s_[::2], np.s_[3:3],
            np.s_[1:4, 2], np.s_[:, 1:3], np.s_[None, 1:4], np.s_[..., 1:], np.s_[1:-1, None, ::-1],
        ]
        for key in keys:
            with self.subTest(key=key):
                sliced = self.disjoint_view[key]
                self.assertIsInstance(sliced, DisjointView)
                np.testing.assert_array_equal(self.expected[key], sliced)
                self.assertEqual(self.expected[key].shape, sliced.shape)

    def test_index_on_axis_returns_view(self):
        for index in [0, 2, 4, -1, -5]:
            with self.subTest(index=index):
                row = self.disjoint_view[index]
                self.assertIsInstance(row, np.ndarray)
                np.testing.assert_array_equal(self.expected[index], row)
        self.disjoint_view[2][1] = -1
        self.assertEqual(-1, self.parent[5, 1])
        with self.assertRaises(IndexError):
            self.disjoint_view[5]

    def test_index_off_axis(self):
        parent = np.arange(40).reshape(4, 10)
        view = DisjointView(parent, [parent[:, 0:2], parent[:, 5:8]], axis=1)
        expected = np.concatenate([parent[:, 0:2], parent[:, 5:8]], axis=1)
        sliced = view[1]
        self.assertIsInstance(sliced, DisjointView)
        self.assertEqual(0, sliced.axis)
        np.testing.assert_array_equal(expected[1], sliced)
        np.testing.assert_array_equal(expected[1, 3], view[1, 3])

    def test_chained_set(self):
        before = self.parent.copy()
        self.disjoint_view[1:4][::-1][0:2] = -1
        before[5:7] = -1
        np.testing.assert_array_equal(before, self.parent)

        self.disjoint_view[1:][1][2] = -2
        before[5, 2] = -2
        np.testing.assert_array_equal(before, self.parent)

    def test_chained_in_place_operator(self):
        before = self.parent.copy()
        sliced = self.disjoint_view[1:3][:, 1:]
        sliced += 10
        before[1, 1:] += 10
        before[5, 1:] += 10
        np.testing.assert_array_equal(before, self.parent)

    def test_invalid_index(self):
        with self.assertRaisesRegex(IndexError, "too many indices"):
            self.disjoint_view[0, 0, 0]
        with self.assertRaisesRegex(IndexError, "single ellipsis"):
            self.disjoint_view[..., ...]

    def test_ufunc_at(self):
        before = self.parent.copy()
        np.add.at(self.disjoint_view, [0, 3, 3], 100)
        before[0] += 100
        before[6] += 200
        np.testing.assert_array_equal(before, self.parent)

    def test_bool(self):
        parent = np.arange(10)
        self.assertFalse(DisjointView(parent, [parent[0:1]], axis=0))
        self.assertTrue(DisjointView(parent, [parent[5:6]], axis=0))
        with self.assertRaises(ValueError):  # numpy behaviour
            bool(self.disjoint_view)

    def test_advanced_index_unsupported(self):
        before = self.parent.copy()
        for key in [[0, 3], np.array([0, 3]), self.expected > 10, True, (slice(None), [0, 1])]:
            with self.subTest(key=key):
                with self.assertRaises(NotImplementedError):
                    self.disjoint_view[key]
                with self.assertRaises(NotImplementedError):
                    self.disjoint_view[key] = -1
        np.testing.assert_array_equal(before, self.parent)

    def test_copy(self):
        copied = np.array(self.disjoint_view, copy=True)
        copied[...] = -1
        np.testing.assert_array_equal(self.expected, self.disjoint_view[...])

    def test_deepcopy(self):
        copied = deepcopy(self.disjoint_view)
        np.testing.assert_array_equal(self.expected, copied[...])

    def test_reversed_regions(self):
        parent = np.arange(10)
        view = DisjointView(parent, [parent[8:6:-1], parent[2::-1]], axis=0)
        np.testing.assert_array_equal([8, 7, 2, 1, 0], view[...])
        view[...] = [-8, -7, -2, -1, -0]
        np.testing.assert_array_equal([0, -1, -2, 3, 4, 5, 6, -7, -8, 9], parent)

    def test_axis(self):
        parent = np.arange(40).reshape(4, 10)
        view = DisjointView(parent, [parent[:, 0:2], parent[:, 5:8]], axis=1)
        np.testing.assert_array_equal(
            np.concatenate([parent[:, 0:2], parent[:, 5:8]], axis=1), view[...]
        )
        view[:, 2] = -1
        np.testing.assert_array_equal(-1, parent[:, 5])


if __name__ == '__main__':
    unittest.main()
