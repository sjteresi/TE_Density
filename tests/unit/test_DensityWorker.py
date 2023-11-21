
import pytest

from transponson.density_worker2 import DensityWorker
from transponson.density_worker2 import DensityWorkerConfig


from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData





def test_start_stop(self):
    """Can we initialize, start, then stop?"""

    worker = DensityWorker(
