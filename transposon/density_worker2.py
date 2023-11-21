

@dataclass
class DensityJob:

    te_names: list(str)  # a subset of transposons to process
    group_type: str  # superfam or order, may not need this one?
    direction: Overlap.Direction


class DensityResult:

    job: DensityJob
    te_name: str
    density: numpy.ndarray  # n_win x n_genes
    direction: Overlap.Direction


@dataclass
class DensityWorkerConfig:

    self.gene_data_path
    self.te_data_path
    self.overlap_data_path
    self.density_data_path  # open this read only!

    self.direction: OverlapDirection =
    self.group_type =  # superfam or order
    self.windows = 

class DensityWorker(Process):

    def __init__(
        self,
        config: DensityWorkerConfig,
        job_queue: multiprocessing.Queue,
        result_queue: multiprocessing.Queue,
    ):
        # TODO check if density data exists, it should exist before this
        # b/c we want to read it from multiple workers, not create one

        # TODO OverlapData will need to be part of the ExitStack
        self.overlap_data = OverlapData.from_file(cfg.overlap_data_path)

        # TODO DensityData will need to be part of the ExitStack
        te_data = TransposonData.read(cfg.te_data_path)
        gene_data = GeneData.read(cfg.gene_data_path)
        self.density_data = DensityData(gene_data, te_data, self.density_data_path)
        self.density_subset = None

        self.job_queue = 
        self.result_queue = 
        self.stop_event = multiprocessing.Event()


    def __enter__(self):
        # self.start()
        # return self

    def __exit__(self,):
        # self.stop()

    def _get_job(self):
        # self.job_queue.get()

    def process_job(self, job: DensityJob):

        overlap_data = self.overlap_data
        overlap_gene_to_idx = self.overlap_data...  # new func
        overlap_window_to_idx = self.overlap_data... # new func
        overlap_sub = self.overlap_data.get_left_intra_right(..)  # new func, l|i|r

        gene_names = self.overlap_data.gene_names
        windows = self.overlap_data.windows

        gene_data = self.density_data.gene_data
        te_data = self.density_data.te_data
        density_data = self.density_data
        density_sub = self.density_data.get_sub(..)  # new func, super||order

        windows = density_data.windows
        if direction == Overlap.Direction.INTRA:
            windows = [0]

        te_array = te_data.get_group(..)  # new fun, super||order
        for te_name in job.te_names:
            group_match = np.equal(te_array, te_name)
            for g_idx, gene_name in enumerate(density_data.gene_names):
                gene_datum = gene_data.get_gene(gene_name)

                # find divisors here using gene datum and direction

                o_gene_idx = overlap_gene_to_idx(gene_name)
                overlap_sub_win_tes = overlap_sub[o_gene_idx, :, :]

                sum_out = numpy.zeros(len(gene_names), len(windows))
                for w_idx, window in enumerate(windows):
                    o_win_idx = overlap_window_to_idx(window, job.direction)
                    overlap = overlap_sub_win_tes[o_win_idx, :]
                    o_sum = np.sum(
                        overlap,
                        where=group_match,
                    )
                    sum_out[g_idx, w_idx] = o_sum


    def overlap_subset(self, overlap_data, window, gene_name):

        # need to get the overlap subarray
        # first, need to find the window_index, gene_index that
        # overlap_data uses



