#!/usr/bin/env python3

"""
Unit test Overlap.
"""

__author__ = "Scott Teresi, Michael Teresi"

import pytest

import numpy as np

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import Overlap

# NOTE these should be similar to the density tests, you are refactoring them sto that
# the scope is smaller so that we can do the pseudo split merge pattern more easily


class MockData(object):

    def __init__(self, g_start, g_stop, t_start, t_stop, window, expected_overlap_left,
                 expected_overlap_right, expected_overlap_intra, description='Test'):

        self.Gene = GeneData.mock(np.array([[g_start, g_stop]]))
        self.Transposon = TransposonData.mock(np.array([[t_start, t_stop]]))
        self.gene_name = 'gene_0'  # TODO scott, add str

        self.window = window
        self.expected_overlap_left = expected_overlap_left
        self.expected_overlap_intra = expected_overlap_intra
        self.expected_overlap_right = expected_overlap_right
        self.description = description

    @property
    def gene_start(self):
        return self.gene.start(self.gene_name)

    @property
    def gene_stop(self):
        return self.gene.stop(self.gene_name)

    @property
    def transposon_start(self):
        return self.Transposon.starts[0]

    @property
    def transposon_stop(self):
        return self.Transposon.stops[0]

Left_Win_Only_FLOAT = MockData(
    g_start=float(1000),
    g_stop=float(2000),
    t_start=float(700),
    t_stop=float(800),
    window=500,
    expected_overlap_left=float(101),
    expected_overlap_intra=float(0),
    expected_overlap_right=float(0),
    description='In left window only')

Left_Win_Only_INT = MockData(
    g_start=int(1000),
    g_stop=int(2000),
    t_start=int(700),
    t_stop=int(800),
    window=500,
    expected_overlap_left=int(101),
    expected_overlap_intra=int(0),
    expected_overlap_right=int(0),
    description='In left window only')

Left_Win_Only_heterogenous_INT_FLOAT = MockData(
    g_start=int(1000),
    g_stop=int(2000),
    t_start=int(700),
    t_stop=int(800),
    window=500,
    expected_overlap_left=float(101),
    expected_overlap_intra=float(0),
    expected_overlap_right=float(0),
    description='In left window only')

Left_Win_Only_heterogenous_FLOAT_INT = MockData(
    g_start=float(1000),
    g_stop=float(2000),
    t_start=float(700),
    t_stop=float(800),
    window=500,
    expected_overlap_left=int(101),
    expected_overlap_intra=int(0),
    expected_overlap_right=int(0),
    description='In left window only')




Left_Win_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=700,
    t_stop=800,
    window=500,
    expected_overlap_left=101,
    expected_overlap_intra=0,
    expected_overlap_right=0,
    description='In left window only')

Right_Win_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2700,
    t_stop=2800,
    window=1000,
    expected_overlap_left=0,
    expected_overlap_intra=0,
    expected_overlap_right=101,
    description='In right window only')

Left_Outside_Win_End_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=400,
    t_stop=800,
    window=500,
    expected_overlap_left=302,
    expected_overlap_intra=0,
    expected_overlap_right=0,
    description='''Left: Start outside window & end
        inside window''')

Right_Inside_Win_End_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2400,
    t_stop=2700,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=0,
    expected_overlap_right=102,
    description='''Right: Begin inside window & end
        outside window''')
# 4
Left_Inside_Win_End_Inside_Gene = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=800,
    t_stop=1500,
    window=500,
    expected_overlap_left=200,
    expected_overlap_intra=501,
    expected_overlap_right=0,
    description='Left: Start in window, end in gene')

Right_Inside_Gene_End_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1500,
    t_stop=2200,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=501,
    expected_overlap_right=200,
    description='Right: Start in gene, end in window')

Inside_Gene_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1500,
    t_stop=1700,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=201,
    expected_overlap_right=0,
    description='In gene only')

# 7
Overlap_Gene_Full = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1000,
    t_stop=2000,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=1001,
    expected_overlap_right=0,
    description='In gene full overlap')

# Exceptional Cases
Left_Outside_Win_End_Inside_Gene = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=1200,
    window=500,
    expected_overlap_left=501,
    expected_overlap_intra=201,
    expected_overlap_right=0,
    description='''Left: Start outside window, end
        inside gene''')
# 9
Right_Inside_Gene_End_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1250,
    t_stop=2600,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=751,
    expected_overlap_right=501,
    description='''Right: Start inside gene, end
        outside window''')

Left_Win_To_Right_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=800,
    t_stop=2200,
    window=500,
    expected_overlap_left=200,
    expected_overlap_intra=1001,
    expected_overlap_right=200,
    description='Cover gene, start and end in windows')

Left_Outside_Win_Right_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=2800,
    window=500,
    expected_overlap_left=501,
    expected_overlap_intra=1001,
    expected_overlap_right=501,
    description='Cover gene, start and end outside both windows')

# 12
Left_Inside_Win_Right_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=700,
    t_stop=2372,
    window=500,
    expected_overlap_left=300,
    expected_overlap_intra=1001,
    expected_overlap_right=372,
    description='Start inside left window, end inside right window')

Left_Inside_Win_Right_Inside_Win_V2 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=701,
    t_stop=2372,
    window=500,
    expected_overlap_left=299,
    expected_overlap_intra=1001,
    expected_overlap_right=372,
    description='Start inside left window, end inside right window')

Left_Outside_Win_Right_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=2248,
    window=500,
    expected_overlap_left=501,
    expected_overlap_intra=1001,
    expected_overlap_right=248,
    description='Start outside left window, end inside right window')

Left_Inside_Win_Right_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=879,
    t_stop=2800,
    window=500,
    expected_overlap_left=121,
    expected_overlap_intra=1001,
    expected_overlap_right=501,
    description='Start inside left window, end outside right window')

# Edge Cases
# 16
Left_On_WinStart_End_On_GeneLeft = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=500,
    t_stop=1000,
    window=500,
    expected_overlap_left=500,
    expected_overlap_intra=1,
    expected_overlap_right=0,
    description='Start inside left window, end outside right window')

Intra_Edge_V1 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=600,
    t_stop=1000,
    window=500,
    expected_overlap_left=400,
    expected_overlap_intra=1,
    expected_overlap_right=0,
    description='Test')

Intra_Edge_V2 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2000,
    t_stop=2200,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=1,
    expected_overlap_right=200,
    description='Test')

Intra_Edge_V3 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1001,
    t_stop=1999,
    window=500,
    expected_overlap_left=0,
    expected_overlap_intra=999,
    expected_overlap_right=0,
    description='Test')


@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only,
                            Right_Win_Only,
                            Left_Outside_Win_End_Inside_Win,
                            Right_Inside_Win_End_Outside_Win,
                            Left_Inside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Inside_Win,
                            Inside_Gene_Only,
                            Overlap_Gene_Full,
                            Left_Outside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Outside_Win,
                            Left_Win_To_Right_Win,
                            Left_Outside_Win_Right_Outside_Win,
                            Left_Inside_Win_Right_Inside_Win,
                            Left_Inside_Win_Right_Inside_Win_V2,
                            Left_Outside_Win_Right_Inside_Win,
                            Left_Inside_Win_Right_Outside_Win,
                            Left_On_WinStart_End_On_GeneLeft,
                            Intra_Edge_V1,
                            Intra_Edge_V2,
                            Intra_Edge_V3
                         ]
                         )
def test_overlap_left_window(MockData_Obj):
    """
    Test the left window
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.left(gene_datum,
                            MockData_Obj.Transposon,
                            MockData_Obj.window)
    expected_overlap_lefts = np.array([MockData_Obj.expected_overlap_left])
    assert np.all(overlaps == expected_overlap_lefts)


@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only,
                            Right_Win_Only,
                            Left_Outside_Win_End_Inside_Win,
                            Right_Inside_Win_End_Outside_Win,
                            Left_Inside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Inside_Win,
                            Inside_Gene_Only,
                            Overlap_Gene_Full,
                            Left_Outside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Outside_Win,
                            Left_Win_To_Right_Win,
                            Left_Outside_Win_Right_Outside_Win,
                            Left_Inside_Win_Right_Inside_Win,
                            Left_Inside_Win_Right_Inside_Win_V2,
                            Left_Outside_Win_Right_Inside_Win,
                            Left_Inside_Win_Right_Outside_Win,
                            Left_On_WinStart_End_On_GeneLeft,
                            Intra_Edge_V1,
                            Intra_Edge_V2,
                            Intra_Edge_V3
                         ]
                         )
def test_overlap_right_window(MockData_Obj):
    """
    Test the right window
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.right(gene_datum,
                             MockData_Obj.Transposon,
                             MockData_Obj.window)
    expected_overlap_rights = np.array([MockData_Obj.expected_overlap_right])
    assert np.all(overlaps == expected_overlap_rights)


@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only,
                            Right_Win_Only,
                            Left_Outside_Win_End_Inside_Win,
                            Right_Inside_Win_End_Outside_Win,
                            Left_Inside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Inside_Win,
                            Inside_Gene_Only,
                            Overlap_Gene_Full,
                            Left_Outside_Win_End_Inside_Gene,
                            Right_Inside_Gene_End_Outside_Win,
                            Left_Win_To_Right_Win,
                            Left_Outside_Win_Right_Outside_Win,
                            Left_Inside_Win_Right_Inside_Win,
                            Left_Outside_Win_Right_Inside_Win,
                            Left_Inside_Win_Right_Inside_Win_V2,
                            Left_Inside_Win_Right_Outside_Win,
                            Left_On_WinStart_End_On_GeneLeft,
                            Intra_Edge_V1,
                            Intra_Edge_V2,
                            Intra_Edge_V3
                         ]
                         )
def test_overlap_intra_window(MockData_Obj):
    """
    Test the right window
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.intra(gene_datum,
                             MockData_Obj.Transposon)
    expected_overlap_intra = np.array([MockData_Obj.expected_overlap_intra])
    assert np.all(overlaps == expected_overlap_intra)



@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only_FLOAT,
                            Left_Win_Only_INT,
                            Left_Win_Only_heterogenous_INT_FLOAT,
                            Left_Win_Only_heterogenous_FLOAT_INT
                         ]
                         )
def test_overlap_left_window_data_types(MockData_Obj):
    """
    LEFT
    Test the data types
    Test to see if uint32 and float32 produce the same number for overlap
    I have a hunch that numbers that are supposed to be 0 get returned
    incorrectly.
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.left(gene_datum,
                            MockData_Obj.Transposon,
                            MockData_Obj.window)
    expected_overlap_lefts = np.array([MockData_Obj.expected_overlap_left])
    assert np.all(overlaps == expected_overlap_lefts)


@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only_FLOAT,
                            Left_Win_Only_INT,
                            Left_Win_Only_heterogenous_INT_FLOAT,
                            Left_Win_Only_heterogenous_FLOAT_INT
                         ]
                         )
def test_overlap_intra_window_data_types(MockData_Obj):
    """
    INTRA
    Test the data types
    Test to see if uint32 and float32 produce the same number for overlap
    I have a hunch that numbers that are supposed to be 0 get returned
    incorrectly.
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.intra(gene_datum,
                            MockData_Obj.Transposon)
    expected_overlap_intra = np.array([MockData_Obj.expected_overlap_intra])
    assert np.all(overlaps == expected_overlap_intra)


@pytest.mark.parametrize("MockData_Obj",
                         [
                            Left_Win_Only_FLOAT,
                            Left_Win_Only_INT,
                            Left_Win_Only_heterogenous_INT_FLOAT,
                            Left_Win_Only_heterogenous_FLOAT_INT
                         ]
                         )
def test_overlap_right_window_data_types(MockData_Obj):
    """
    RIGHT
    Test the data types
    Test to see if uint32 and float32 produce the same number for overlap
    I have a hunch that numbers that are supposed to be 0 get returned
    incorrectly.
    """
    gene_datum = MockData_Obj.Gene.get_gene('gene_0')
    overlaps = Overlap.right(gene_datum,
                            MockData_Obj.Transposon,
                            MockData_Obj.window)
    expected_overlap_rights = np.array([MockData_Obj.expected_overlap_right])
    assert np.all(overlaps == expected_overlap_rights)

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
