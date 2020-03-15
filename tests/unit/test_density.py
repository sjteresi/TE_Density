#!/usr/bin/env python3

"""
Unit test density calculations.

"""

__author__ = "Michael Teresi, Scott Teresi"

import pytest
import numpy as np
import pandas as pd

from transposon.density import rho_intra
from transposon.density import rho_left_window
from transposon.density import rho_right_window
# from transposon.density import validate_window
from transposon.data import GeneData, TransposonData


def mock_te_data(start_stop):
    """Transposon data for given the start/stop indices
        Creates one

        Args:
            start_stop (np.array): N gene x (start_idx, stop_idx). 2D array
    """

    n_genes = start_stop.shape[0]
    data = []
    family = "Family_0"  # FUTURE may want to parametrize family name later
    # NB overall order is not important but the names are
    columns = ['Start', 'Stop', 'Length', 'Order', 'SuperFamily', 'Chromosome']
    for gi in range(n_genes):
        g0 = start_stop[gi, 0]
        g1 = start_stop[gi, 1]
        gL = g1 - g0 + 1
        # FUTURE may want to parametrize sub
        # family name later
        subfam_suffix = "A" if gi % 2 else "B"
        subfamily = "SubFamily_{}".format(subfam_suffix)
        chromosome = 'Chr_Test'
        datum = [g0, g1, gL, family, subfamily, chromosome]
        data.append(datum)
    frame = pd.DataFrame(data, columns=columns)
    return TransposonData(frame)


class MockData(object):

    def __init__(self, g_start, g_stop, t_start, t_stop, window, expected_rho_left,
                 expected_rho_right, expected_rho_intra, description='Test'):

        self.Gene = GeneData.mock(np.array([[g_start, g_stop]]))
        self.Transposon = mock_te_data(np.array([[t_start, t_stop]]))
        self.gene_name = 'gene_0'  # TODO scott, add str

        self.window = window
        self.expected_rho_left = expected_rho_left
        self.expected_rho_intra = expected_rho_intra
        self.expected_rho_right = expected_rho_right
        self.description = description

        # NOTE add expected overlap / density funcs?

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


Left_Win_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=700,
    t_stop=800,
    window=500,
    expected_rho_left=101/501,
    expected_rho_intra=0,
    expected_rho_right=0,
    description='In left window only')

Right_Win_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2700,
    t_stop=2800,
    window=1000,
    expected_rho_left=0,
    expected_rho_intra=0,
    expected_rho_right=101/1001,
    description='In right window only')

Left_Outside_Win_End_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=400,
    t_stop=800,
    window=500,
    expected_rho_left=302/501,
    expected_rho_intra=0,
    expected_rho_right=0,
    description='''Left: Start outside window & end
        inside window''')

Right_Inside_Win_End_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2400,
    t_stop=2700,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=0,
    expected_rho_right=np.divide(102, 501),
    description='''Right: Begin inside window & end
        outside window''')
# 4
Left_Inside_Win_End_Inside_Gene = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=800,
    t_stop=1500,
    window=500,
    expected_rho_left=200/501,
    expected_rho_intra=501/1001,
    expected_rho_right=0,
    description='Left: Start in window, end in gene')

Right_Inside_Gene_End_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1500,
    t_stop=2200,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=501/1001,
    expected_rho_right=200/501,
    description='Right: Start in gene, end in window')

Inside_Gene_Only = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1500,
    t_stop=1700,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=201/1001,
    expected_rho_right=0,
    description='In gene only')

# 7
Overlap_Gene_Full = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1000,
    t_stop=2000,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=1,
    expected_rho_right=0,
    description='In gene full overlap')

# Exceptional Cases
Left_Outside_Win_End_Inside_Gene = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=1200,
    window=500,
    expected_rho_left=1,
    expected_rho_intra=201/1001,
    expected_rho_right=0,
    description='''Left: Start outside window, end
        inside gene''')
# 9
Right_Inside_Gene_End_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1250,
    t_stop=2600,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=np.divide(751, 1001),
    expected_rho_right=1,
    description='''Right: Start inside gene, end
        outside window''')

Left_Win_To_Right_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=800,
    t_stop=2200,
    window=500,
    expected_rho_left=200/501,
    expected_rho_intra=1,
    expected_rho_right=200/501,
    description='Cover gene, start and end in windows')

Left_Outside_Win_Right_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=2800,
    window=500,
    expected_rho_left=1,
    expected_rho_intra=1,
    expected_rho_right=1,
    description='Cover gene, start and end outside both windows')
# 12
Left_Inside_Win_Right_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=700,
    t_stop=2372,
    window=500,
    expected_rho_left=300/501,
    expected_rho_intra=1,
    expected_rho_right=372/501,
    description='Start inside left window, end inside right window')

Left_Inside_Win_Right_Inside_Win_V2 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=701,
    t_stop=2372,
    window=500,
    expected_rho_left=299/501,
    expected_rho_intra=1,
    expected_rho_right=372/501,
    description='Start inside left window, end inside right window')

Left_Outside_Win_Right_Inside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=300,
    t_stop=2248,
    window=500,
    expected_rho_left=1,
    expected_rho_intra=1,
    expected_rho_right=np.divide(248, 501),
    description='Start outside left window, end inside right window')

Left_Inside_Win_Right_Outside_Win = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=879,
    t_stop=2800,
    window=500,
    expected_rho_left=121/501,
    expected_rho_intra=1,
    expected_rho_right=1,
    description='Start inside left window, end outside right window')

# Edge Cases
# 16
Left_On_WinStart_End_On_GeneLeft = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=500,
    t_stop=1000,
    window=500,
    expected_rho_left=500/501,
    expected_rho_intra=np.divide(1, 1001),
    expected_rho_right=0,
    description='Start inside left window, end outside right window')

Intra_Edge_V1 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=600,
    t_stop=1000,
    window=500,
    expected_rho_left=400/501,
    expected_rho_intra=np.divide(1, 1001),
    expected_rho_right=0,
    description='Test')

Intra_Edge_V2 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=2000,
    t_stop=2200,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=np.divide(1, 1001),
    expected_rho_right=200/501,
    description='Test')

Intra_Edge_V3 = MockData(
    g_start=1000,
    g_stop=2000,
    t_start=1001,
    t_stop=1999,
    window=500,
    expected_rho_left=0,
    expected_rho_intra=np.divide(999, 1001),
    expected_rho_right=0,
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
                            Intra_Edge_V2
                         ]
                         )
def test_rho_left_window(MockData_Obj):
    """
    Test the left window
    """
    rhos = rho_left_window(MockData_Obj.Gene,
                           'gene_0',  # fake name needed because of GeneData
                           MockData_Obj.Transposon,
                           MockData_Obj.window)
    expected_rho_lefts = np.array([MockData_Obj.expected_rho_left])
    assert np.all(rhos == expected_rho_lefts)


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
                            Left_Inside_Win_Right_Outside_Win,
                            Left_On_WinStart_End_On_GeneLeft,
                            Intra_Edge_V1,
                            Intra_Edge_V2,
                            Intra_Edge_V3
                         ]
                         )
def test_rho_right_window(MockData_Obj):
    """
    Test the right window
    """
    rhos = rho_right_window(MockData_Obj.Gene,
                            'gene_0',  # fake name needed because of GeneData
                            MockData_Obj.Transposon,
                            MockData_Obj.window)
    expected_rho_rights = np.array([MockData_Obj.expected_rho_right])
    assert np.all(rhos == expected_rho_rights)


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
                            Left_Inside_Win_Right_Outside_Win,
                            Left_On_WinStart_End_On_GeneLeft,
                            Intra_Edge_V1,
                            Intra_Edge_V2,
                            Intra_Edge_V3
                         ]
                         )
def test_rho_intra(MockData_Obj):
    """
    Test intra density
    """
    rhos = rho_intra(MockData_Obj.Gene,
                     'gene_0',  # fake name needed because of GeneData
                     MockData_Obj.Transposon)
    expected_rho_intra = np.array([MockData_Obj.expected_rho_intra])
    assert np.all(rhos == expected_rho_intra)


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
