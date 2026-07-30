from copy import deepcopy
from genome_kit import Genome, Interval, IntervalData
from genome_kit.diseq import DisjointIntervalSequence
from tests import MiniGenome
import numpy as np
import unittest

nucleotide_complement = str.maketrans('ACGT', 'TGCA')
REFG = "hg19.mini"


def _make_intervals(specs, refg=REFG):
    """Build Intervals from a list of ``(chrom, strand, start, end)`` specs."""
    return [
        Interval(chrom, strand, start, end, refg) for chrom, strand, start, end in specs
    ]


def _dis_cases(length):
    """Return ``(label, dis)`` pairs realizing a ``length``-base segment in
    several :py:class:`~genome_kit.DisjointIntervalSequence` configurations.
    """
    cases = []
    for strand in ("+", "-"):
        # Contiguous segment interior to a single coord interval, with room to
        # expand/contract on both sides.
        single = _make_intervals([("chr1", strand, 100, 200)])
        cases.append(
            (
                f"single_{strand}_coord_iv",
                DisjointIntervalSequence(single, start=50, end=50 + length),
            )
        )
        cases.append(
            (
                f"single_{strand}_coord_iv_dis_off_coord_strand",
                DisjointIntervalSequence(single, start=50, end=50 + length, on_coordinate_strand=False),
            )
        )
    if length >= 2:
        for strand in ("+", "-"):
            # Two coord intervals with a gap. The first is 10 bases, so a segment
            # starting at index 9 straddles the internal boundary at index 10.
            double_coord_iv = _make_intervals(
                [("chr1", strand, 100, 110), ("chr1", strand, 200, 210)]
            )
            cases.append(
                (
                    f"spans_2_{strand}_coord_ivs",
                    DisjointIntervalSequence(double_coord_iv, start=9, end=9 + length),
                )
            )
    return cases


def _slice_index_cases(rank):
    """Slice-index cases exercised by the ``[start:stop:step]`` tests."""
    return [('empty', dict(start=0, stop=0)),
            ('none', dict()),
            ('full', dict(start=0, stop=rank, step=1)),
            ('stop_full', dict(stop=rank)),
            ('stop_middle', dict(stop=1)),
            ('stop_after_range', dict(stop=rank * 2)),
            ('stop_negative', dict(stop=-1)),
            ('stop_negative_before_range', dict(stop=-(rank * 2))),
            ('start_full', dict(start=0)),
            ('start_middle', dict(start=1)),
            ('start_after_range', dict(start=rank * 2)),
            ('start_negative', dict(start=-1)),
            ('start_negative_before_range', dict(start=-(rank * 2))),
            ('step_forward', dict(step=1)),
            ('step_reverse', dict(step=-1)),
            ('step_discontinuous', dict(step=2, throws=True)),
            ('reverse_segment', dict(start=1, step=-1)),
            ('reverse_negative_segment', dict(start=-1, stop=1, step=-1)),
            ('reverse_outside_range', dict(start=rank * 2, stop=-(rank * 2), step=-1))]


class TestIntervalData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = MiniGenome(REFG)

    def test_same_dimension(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(3, 10)
        interval_data = IntervalData(interval, data)
        self.assertEqual(len(interval), len(interval_data))
        self.assertEqual(interval, interval_data._interval)
        np.testing.assert_equal(data, [x for x in interval_data])

    def test_same_dimension_dis(self):
        data = np.arange(0, 30).reshape(3, 10)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                self.assertEqual(len(dis), len(interval_data))
                self.assertEqual(dis, interval_data._interval)
                np.testing.assert_equal(data, [x for x in interval_data])

    def test_empty(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval_data = IntervalData(position, [])
        self.assertEqual(0, len(interval_data))
        self.assertEqual(0, len([x for x in interval_data]))

    def test_empty_dis(self):
        for name, dis in _dis_cases(0):
            with self.subTest(name):
                interval_data = IntervalData(dis, [])
                self.assertEqual(0, len(interval_data))
                self.assertEqual(0, len([x for x in interval_data]))

    def test_invalid_dimension(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        with self.assertRaises(ValueError):
            IntervalData(position.expand(1), [1])

    def test_invalid_dimension_dis(self):
        for name, dis in _dis_cases(2):
            with self.subTest(name):
                with self.assertRaises(ValueError):
                    IntervalData(dis, [1])

    def test_outside_range(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(3, 10)
        interval_data = IntervalData(interval, data)
        with self.assertRaises(IndexError):
            interval_data[len(data)]

    def test_outside_range_dis(self):
        data = np.arange(0, 30).reshape(3, 10)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                with self.assertRaises(IndexError):
                    interval_data[len(data)]

    def test_set(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(3, 10)
        interval_data = IntervalData(interval, deepcopy(data))
        interval_data[0] = interval_data[0].transpose()
        np.testing.assert_equal(data[0].transpose(), interval_data[0])

    def test_set_dis(self):
        data = np.arange(0, 30).reshape(3, 10)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, deepcopy(data))
                interval_data[0] = interval_data[0].transpose()
                np.testing.assert_equal(data[0].transpose(), interval_data[0])


class TestIntervalDataSliceIndex(unittest.TestCase):
    rank = 3

    @classmethod
    def setUpClass(cls):
        cls.genome = MiniGenome(REFG)

    def test_slice_index(self):
        rank = self.rank
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, rank)
        dna = self.genome.dna(interval)
        data = np.arange(0, 10 * rank).reshape(rank, 10)
        interval_data = IntervalData(interval, data)

        for name, kwargs in _slice_index_cases(rank):
            with self.subTest(name):
                start = kwargs.get('start')
                stop = kwargs.get('stop')
                step = kwargs.get('step')
                if kwargs.get('throws'):
                    with self.assertRaises(KeyError):
                        interval_data[start:stop:step]
                else:
                    sliced = interval_data[start:stop:step]
                    expected_dna = dna[start:stop:step]
                    if step is not None and step < 0:
                        expected_dna = expected_dna.translate(nucleotide_complement)
                    self.assertEqual(expected_dna, self.genome.dna(sliced.interval))
                    np.testing.assert_equal(
                        [x for x in data[start:stop:step]], [x for x in sliced]
                    )

    def test_slice_index_dis(self):
        rank = self.rank
        data = np.arange(0, 10 * rank).reshape(rank, 10)

        for dis_name, dis in _dis_cases(rank):
            interval_data = IntervalData(dis, data)
            dna = dis.dna()
            for name, kwargs in _slice_index_cases(rank):
                with self.subTest(f"{dis_name}:{name}"):
                    start = kwargs.get('start')
                    stop = kwargs.get('stop')
                    step = kwargs.get('step')
                    if kwargs.get('throws'):
                        with self.assertRaises(KeyError):
                            interval_data[start:stop:step]
                    else:
                        sliced = interval_data[start:stop:step]
                        expected_dna = dna[start:stop:step]
                        if step is not None and step < 0:
                            expected_dna = expected_dna.translate(nucleotide_complement)
                        self.assertEqual(expected_dna, sliced.interval.dna())
                        np.testing.assert_equal(
                            [x for x in data[start:stop:step]], [x for x in sliced]
                        )

    def test_negative_strand_negative_step(self):
        rank = self.rank
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, rank).as_opposite_strand()
        data = np.arange(0, 10 * rank).reshape(rank, 10)
        start = None
        stop = 1
        step = -1  # Will flip the strand and reverse the slice
        interval_data = IntervalData(interval, data)[start:stop:step]
        dna = self.genome.dna(interval)[start:stop:step].translate(nucleotide_complement)
        expected_data = data[start:stop:step]
        self.assertEqual(dna, self.genome.dna(interval_data.interval))
        np.testing.assert_equal([x for x in expected_data], [x for x in interval_data])

    def test_negative_strand_negative_step_dis(self):
        rank = self.rank
        data = np.arange(0, 10 * rank).reshape(rank, 10)
        start = None
        stop = 1
        step = -1  # Will flip the strand and reverse the slice
        for name, dis in _dis_cases(rank):
            with self.subTest(name):
                backing = dis.as_opposite_strand()
                interval_data = IntervalData(backing, data)[start:stop:step]
                dna = backing.dna()[start:stop:step].translate(nucleotide_complement)
                expected_data = data[start:stop:step]
                self.assertEqual(dna, interval_data.interval.dna())
                np.testing.assert_equal(
                    [x for x in expected_data], [x for x in interval_data]
                )

    def test_negative_strand(self):
        rank = self.rank
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, rank).as_opposite_strand()
        data = np.arange(0, 10 * rank).reshape(rank, 10)
        start = 2
        stop = None
        step = 1
        interval_data = IntervalData(interval, data)[start:stop:step]
        dna = self.genome.dna(interval)[start:stop:step]
        expected_data = data[start:stop:step]
        self.assertEqual(dna, self.genome.dna(interval_data.interval))
        np.testing.assert_equal([x for x in expected_data], [x for x in interval_data])

    def test_negative_strand_dis(self):
        rank = self.rank
        data = np.arange(0, 10 * rank).reshape(rank, 10)
        start = 2
        stop = None
        step = 1
        for name, dis in _dis_cases(rank):
            with self.subTest(name):
                backing = dis.as_opposite_strand()
                interval_data = IntervalData(backing, data)[start:stop:step]
                dna = backing.dna()[start:stop:step]
                expected_data = data[start:stop:step]
                self.assertEqual(dna, interval_data.interval.dna())
                np.testing.assert_equal(
                    [x for x in expected_data], [x for x in interval_data]
                )


class TestIntervalDataSliceInterval(unittest.TestCase):
    rank = 3

    @classmethod
    def setUpClass(cls):
        cls.genome = MiniGenome(REFG)

    def test_empty(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(interval, data)
        sliced = interval_data[position]
        self.assertEqual(position, sliced.interval)
        self.assertEqual(0, len([x for x in sliced]))

    def test_empty_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                sliced = interval_data[dis.end5]
                self.assertEqual(dis.end5, sliced.interval)
                self.assertEqual(0, len([x for x in sliced]))

    def test_full(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(interval, data)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data, [x for x in sliced])

    def test_full_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                sliced = interval_data[dis]
                self.assertEqual(dis, sliced.interval)
                np.testing.assert_equal(data, [x for x in sliced])

    def test_outside_range(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(interval, data)
        with self.assertRaises(IndexError):
            interval_data[interval.shift(1)]

    def test_outside_range_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                with self.assertRaises(IndexError):
                    interval_data[dis.shift(1)]

    def test_start_middle(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(base, data)
        interval = base.expand(-1, 0)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[1:], [x for x in sliced])

    def test_start_middle_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                key = dis.expand(-1, 0)
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[1:], [x for x in sliced])

    def test_end_middle(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(base, data)
        interval = base.expand(0, -1)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[:-1], [x for x in sliced])

    def test_end_middle_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                key = dis.expand(0, -1)
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[:-1], [x for x in sliced])

    def test_reverse_full(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(base, data)
        interval = base.as_opposite_strand()
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[::-1], [x for x in sliced])

    def test_reverse_full_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                key = dis.as_opposite_strand()
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[::-1], [x for x in sliced])

    def test_reverse_start_middle(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(base, data)
        interval = base.as_opposite_strand().expand(-1, 0)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[-2::-1], [x for x in sliced])

    def test_reverse_start_middle_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                key = dis.as_opposite_strand().expand(-1, 0)
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[-2::-1], [x for x in sliced])

    def test_reverse_end_middle(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(base, data)
        interval = base.as_opposite_strand().expand(0, -1)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[:0:-1], [x for x in sliced])

    def test_reverse_end_middle_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, data)
                key = dis.as_opposite_strand().expand(0, -1)
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[:0:-1], [x for x in sliced])

    def test_negative_strand(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, self.rank).as_opposite_strand()
        data = range(len(interval))
        interval_data = IntervalData(interval, data)
        np.testing.assert_equal(interval_data[interval.end5.expand(0, 1)].data, data[:1])
        np.testing.assert_equal(interval_data[interval].data, data)

    def test_negative_strand_dis(self):
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                backing = dis.as_opposite_strand()
                data = range(len(backing))
                interval_data = IntervalData(backing, data)
                np.testing.assert_equal(
                    interval_data[backing.end5.expand(0, 1)].data, data[:1]
                )
                np.testing.assert_equal(interval_data[backing].data, data)

    def test_set_data(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, self.rank)
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        interval_data = IntervalData(interval, deepcopy(data))
        interval_data[interval.expand(-1, 0)] = 0
        np.testing.assert_array_equal(
            interval_data._data[1:, :], np.zeros_like(interval_data._data[1:, :])
        )
        np.testing.assert_array_equal(interval_data._data[:1, :], data[:1, :])

    def test_set_data_dis(self):
        data = np.arange(0, 10 * self.rank).reshape(self.rank, 10)
        for name, dis in _dis_cases(self.rank):
            with self.subTest(name):
                interval_data = IntervalData(dis, deepcopy(data))
                interval_data[dis.expand(-1, 0)] = 0
                np.testing.assert_array_equal(
                    interval_data._data[1:, :],
                    np.zeros_like(interval_data._data[1:, :]),
                )
                np.testing.assert_array_equal(interval_data._data[:1, :], data[:1, :])


class TestIntervalDataAxisAlign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = MiniGenome(REFG)

    def test_same_dimension(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, data, axis=1)
        self.assertEqual(len(interval), len(interval_data))
        self.assertEqual(interval, interval_data._interval)
        np.testing.assert_equal(data, [x for x in interval_data])

    def test_same_dimension_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                self.assertEqual(len(dis), len(interval_data))
                self.assertEqual(dis, interval_data._interval)
                np.testing.assert_equal(data, [x for x in interval_data])

    def test_negative_axis(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        with self.assertRaises(IndexError):
            IntervalData(interval, data, -1)

    def test_negative_axis_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                with self.assertRaises(IndexError):
                    IntervalData(dis, data, -1)

    def test_invalid_dimensions(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        with self.assertRaises(ValueError):
            IntervalData(interval, data)
        with self.assertRaises(ValueError):
            IntervalData(interval, data, 2)

    def test_invalid_dimensions_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                with self.assertRaises(ValueError):
                    IntervalData(dis, data)
                with self.assertRaises(ValueError):
                    IntervalData(dis, data, 2)

    def test_outside_range(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, data, axis=1)
        with self.assertRaises(IndexError):
            interval_data[0, data.shape[1], 0]

    def test_outside_range_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                with self.assertRaises(IndexError):
                    interval_data[0, data.shape[1], 0]

    def test_invalid_slice(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, data, axis=1)
        with self.assertRaises(IndexError):
            interval_data[1:]

    def test_invalid_slice_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                with self.assertRaises(IndexError):
                    interval_data[1:]

    def test_invalid_slice_set(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, data, axis=1)
        with self.assertRaises(IndexError):
            interval_data[1:] = 0

    def test_invalid_slice_set_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                with self.assertRaises(IndexError):
                    interval_data[1:] = 0

    def test_slice_index(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, data, axis=1)
        start = 1
        sliced = interval_data[:, start:, ...]
        self.assertEqual(interval.expand(-1, 0), sliced.interval)
        np.testing.assert_equal(data[:, start:, ...], [x for x in sliced])

    def test_slice_index_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        start = 1
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                sliced = interval_data[:, start:, ...]
                self.assertEqual(dis.expand(-1, 0), sliced.interval)
                np.testing.assert_equal(data[:, start:, ...], [x for x in sliced])

    def test_slice_interval(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        base = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(base, data, axis=1)
        interval = base.expand(-1, 0)
        sliced = interval_data[interval]
        self.assertEqual(interval, sliced.interval)
        np.testing.assert_equal(data[:, 1:, ...], [x for x in sliced])

    def test_slice_interval_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, data, axis=1)
                key = dis.expand(-1, 0)
                sliced = interval_data[key]
                self.assertEqual(key, sliced.interval)
                np.testing.assert_equal(data[:, 1:, ...], [x for x in sliced])

    def test_set_data(self):
        position = Interval('chr1', '+', 100, 100, self.genome)
        interval = position.expand(0, 3)
        data = np.arange(0, 30).reshape(5, 3, 2)
        interval_data = IntervalData(interval, deepcopy(data), axis=1)
        interval_data[interval.expand(-1, 0)] = 0
        np.testing.assert_array_equal(
            interval_data._data[:, 1:, :], np.zeros_like(interval_data._data[:, 1:, :])
        )
        np.testing.assert_array_equal(interval_data._data[:, :1, :], data[:, :1, :])

    def test_set_data_dis(self):
        data = np.arange(0, 30).reshape(5, 3, 2)
        for name, dis in _dis_cases(3):
            with self.subTest(name):
                interval_data = IntervalData(dis, deepcopy(data), axis=1)
                interval_data[dis.expand(-1, 0)] = 0
                np.testing.assert_array_equal(
                    interval_data._data[:, 1:, :],
                    np.zeros_like(interval_data._data[:, 1:, :]),
                )
                np.testing.assert_array_equal(interval_data._data[:, :1, :], data[:, :1, :])


if __name__ == '__main__':
    unittest.main()
