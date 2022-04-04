#!/usr/bin/env python

"""
Sums and contains OverlapData with respect to the genes / transposons.
"""

__author__ = "Michael Teresi, Scott Teresi"

from collections import namedtuple
import logging
import os
from functools import partial
import random

import h5py
import numpy as np

import transposon  # TODO move all the functions in __init__ to utils and have an empty __init__ (and change all the imports)
from transposon.gene_datum import GeneDatum


_MergeConfigSink = namedtuple(
    "_MergeConfigSink",
    ["transposons", "gene_data", "gene_names", "windows", "filepath", "ram_bytes"],
)
_MergeConfigSource = namedtuple("_MergeConfigSource", ["filepath"])
_Density = namedtuple("_Density", ["left", "intra", "right"])
_SummationArgs = namedtuple(
    "_SummationArgs",
    [
        "input",
        "output",
        "windows",
        "te_idx_name",
        "slice_in",
        "slice_out",
        "where",
        "divisor_func",
    ],
)


# class MergeCommand:  # TODO refactor the list density args into a command
#
#    DIV_LEFT = GeneDatum.divisor_left
#    DIV_INTRA = GeneDatum.divisor_intra
#    DIV_RIGHT = GeneDatum.divisor_right
#
#    def __init__(self, overlap, transposons, gene_idx, window_idx):
#        """"""
#
#    def exec_superfam(self, superfam_idx):
#        """"""
#
#    def exec_order(self, order_idx):
#        """"""


class MergeData:
    """Contains density of superfamilies and orders.

    Merges from one or more OverlapData into two sets: superfamilies, orders.
    The result has the following dimensions (not necessarily in this order):
        superfamilies x windows x genes x 3
        orders        x windows x genes x 3
    Where the three is for the left, intra, and right densities.

    Usage:
        Use the factory methods to create an instance (`from_param` | `from_file`).
        Use the context manager to activate the buffer.
        Use the slice methods to provide slices for the array public instance variables.
        Use the sum
    """

    DTYPE = np.float64
    _GENE_NAMES = "GENE_NAMES"
    _WINDOWS = "WINDOWS"
    _CHROME_ID = "CHROMOSOME_ID"
    _SUPERFAMILY_NAMES = "SUPERFAMILY_NAMES"
    _ORDER_NAMES = "ORDER_NAMES"
    _RHO_SUPERFAMILY_DATASET_KEY = "RHO_SUPERFAMILIES"
    _S_LEFT = _RHO_SUPERFAMILY_DATASET_KEY + "_LEFT"
    _S_RIGHT = _RHO_SUPERFAMILY_DATASET_KEY + "_RIGHT"
    _S_INTRA = _RHO_SUPERFAMILY_DATASET_KEY + "_INTRA"
    _RHO_ORDER_DATASET_KEY = "RHO_ORDERS"
    _O_LEFT = _RHO_ORDER_DATASET_KEY + "_LEFT"
    _O_RIGHT = _RHO_ORDER_DATASET_KEY + "_RIGHT"
    _O_INTRA = _RHO_ORDER_DATASET_KEY + "_INTRA"

    def __init__(self, configuration, compression="lzf", logger=None):
        """Initializer.

        Args:
            configuration (tuple): _MergeConfigSink or _MergeConfigSource
            compression (str): hdf5 compression type
            logger (logging.Logger): instance for log messages
        """

        self._logger = logger or logging.getLogger(__name__)
        self._config = configuration
        self._compression = self._check_compression(compression, self._logger)
        self._h5_file = None

        self.chromosome_id = None
        self.windows = None
        self.gene_names = None
        self.superfamily_names = None
        self.order_names = None
        self._window_2_idx = None
        self._gene_2_idx = None
        self._superfam_2_idx = None
        self._order_2_idx = None

        self.superfamily = None
        self.order = None

    @staticmethod
    def _check_compression(filter, logger):
        """Raise if compression type is invalid."""

        LOSSLESS_FILTERS = ["gzip", "lzf", "szip"]  # SEE h5py datasets
        if filter not in LOSSLESS_FILTERS:
            msg = "'%s' compression invalid, options are: %s" % (
                filter,
                LOSSLESS_FILTERS,
            )
            logger.critical(msg)
            raise ValueError(msg)
        else:
            return filter

    def __add__(self, other):
        """Combine density files."""
        raise NotImplementedError()

    @property
    def filepath(self):
        """Return the active filepath, None if no open file."""

        return self._h5_file.filename if self._h5_file is not None else None

    @classmethod
    def from_param(
        cls, transposon_data, gene_data, windows, output_dir, ram=1, logger=None
    ):
        """Writable sink for a new file.

        Args:
            transposon_data (TransposonData): transposon container.
            gene_data (GeneData): gene data container
            gene_names (iterable(str)): gene names to process.
            windows (iterable(int)): window inputs to process (not the same as windows).
            output_dir (str): directory to output merge data files.
            ram (int): gigabytes RAM to cache during processing.
            logger (logging.Logger): logging instance for messages.
        """

        logger = logger or logging.getLogger(__name__)
        bytes2gigabytes = 1024.0 ** 3
        ram_bytes = int(ram * bytes2gigabytes)
        transposon.check_ram(ram_bytes, logger)
        chrome = str(transposon_data.chromosome_unique_id)
        genome = str(transposon_data.genome_id)
        # MAGIC, filename creation
        filename = genome + "_" + chrome + ".h5"
        filepath = os.path.join(output_dir, filename)
        config = _MergeConfigSink(
            transposons=transposon_data,
            gene_data=gene_data,
            gene_names=gene_data.names,
            windows=windows,
            filepath=filepath,
            ram_bytes=ram_bytes,
        )
        return cls(config, logger=logger)

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Read only source for an existing file."""

        raise NotImplementedError()

    @classmethod
    def left_right_slice(cls, group_idx=None, gene_idx=None, window_idx=None):
        """Slice to a left || right density to a superfamily|order, a window, gene(s).

        One can use a slice like so: np.array[my_slice].
        This is useful to obtain views into the internal data,
        without hard-coding the order of the arrays.

        Args:
            group_idx (int): superfamily | order index (required)
            gene_idx(int|Slice): gene name index, all genes if None
            window_idx(int|Slice): window index (required), all windows if None
        Returns:
            tuple(int|Slice): integers and/or Slice instances for numpy basic indexing
        """
        # NOTE this is confusing because input arguments are group_idx,
        # gene_idx, window_idx, and it returns (group_idx, window_idx,
        # gene_idx). Potential refactor?

        # SEE numpy basic indexing
        # TODO SCOTT test this pls
        # NB return at least one slice so that one can obtain a view into the array,
        # rather than just the value
        if isinstance(group_idx, int):
            group_idx = slice(group_idx, group_idx + 1, 1)
        if isinstance(window_idx, int):
            window_idx = slice(window_idx, window_idx + 1, 1)
        if isinstance(gene_idx, int):
            gene_idx = slice(gene_idx, gene_idx + 1, 1)
        return (group_idx, window_idx, gene_idx)

    @classmethod
    def intra_slice(cls, group_idx=None, gene_idx=None, window_idx=None):
        """Slice for an intra density to a superfamily|order / window, gene."""

        if window_idx is not None:
            # there is only 1 window for intra operations, and it is None
            raise ValueError("intra window must be None but is %s" % window_idx)

        # TODO SCOTT test this pls
        w_slice = slice(
            0, 1, 1
        )  # use a slice (not 0) so you can get a view of the array
        return cls.left_right_slice(
            group_idx=group_idx, gene_idx=gene_idx, window_idx=w_slice
        )

    def __enter__(self):
        """Context manager begin."""

        self._open_dispatcher()
        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        """Context manager end.

        Writes all the data to the disk.
        """

        self._h5_file.flush()
        self._h5_file.close()
        self._h5_file = None

        self.chromosome_id = None
        self.windows = None
        self.gene_names = None
        self.superfamily_names = None
        self.order_names = None
        self._window_2_idx = None
        self._gene_2_idx = None
        self._superfam_2_idx = None
        self._order_2_idx = None

        self.superfamily = None
        self.order = None

    def _open_dispatcher(self):
        """Open the file.

        NOTE when upgrading to python 3.8 refactor to use functools.singledispatchmethod.
        """

        if isinstance(self._config, _MergeConfigSink):
            self._open_new_file(self._config)
        elif isinstance(self._config, _MergeConfigSource):
            self._open_existing_file(self._config)
        else:
            raise TypeError(
                "expecting {} or {} but got {}".format(
                    type(_MergeConfigSink), type(_MergeConfigSource), type(self._config)
                )
            )

    def _open_new_file(self, cfg):
        """Opens a new H5 file for writing and initializes the data sets.

        Args:
            cfg(_MergeConfigSink): configuration parameters
        """

        self.windows = list(cfg.windows)
        self.gene_names = list(cfg.gene_names)
        self.chromosome_id = cfg.transposons.chromosome_unique_id
        self.superfamily_names = sorted(cfg.transposons.superfamily_name_set)
        self._superfam_2_idx = {s: i for i, s in enumerate(self.superfamily_names)}
        self.order_names = sorted(cfg.transposons.order_name_set)
        self._order_2_idx = {o: i for i, o in enumerate(self.order_names)}
        self._h5_file = h5py.File(cfg.filepath, "w", rdcc_nbytes=cfg.ram_bytes)
        self._create_sets(self._h5_file, cfg)
        transposon.write_vlen_str_h5py(self._h5_file, self.windows, self._WINDOWS)
        transposon.write_vlen_str_h5py(self._h5_file, self.gene_names, self._GENE_NAMES)
        transposon.write_vlen_str_h5py(
            self._h5_file, self.chromosome_id, self._CHROME_ID
        )
        transposon.write_vlen_str_h5py(
            self._h5_file, self.order_names, self._ORDER_NAMES
        )
        transposon.write_vlen_str_h5py(
            self._h5_file, self.superfamily_names, self._SUPERFAMILY_NAMES
        )
        self._window_2_idx = {w: i for i, w in enumerate(self.windows)}
        self._gene_2_idx = {g: i for i, g in enumerate(self.gene_names)}

    # TODO get the index of gene function in here too
    # Clean up/
    def _post_process_swap(self):
        """
        Swap the TE density values for antisense genes ONLY.

    #def _swap_strand_vals(self, gene_names):
        Switch density values for the genes in which it is antisense due to
        the fact that antisense genes point in the opposite direction to sense
        genes

        Args:
            gene_names(list of str):
        """
        for name in gene_names:
            index_to_switch = self._index_of_gene(name)

            # SWAP left and right superfamily values for antisense genes
            (
                self.data_frame["RHO_SUPERFAMILIES_LEFT"][:, :, index_to_switch],
                self.data_frame["RHO_SUPERFAMILIES_RIGHT"][:, :, index_to_switch],
            ) = (
                self.data_frame["RHO_SUPERFAMILIES_RIGHT"][:, :, index_to_switch],
                self.data_frame["RHO_SUPERFAMILIES_LEFT"][:, :, index_to_switch],
            )

            # SWAP left and right order values for antisense genes
            (
                self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch],
                self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch],
            ) = (
                self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch],
                self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch],
            )


    # TODO possibly refactor with DensityData
    def _open_existing_file(self, cfg):

        self._h5_file = h5py.File(cfg.filepath, "r")
        transposon.read_vlen_str_h5py(
            self._h5_file,
        )
        # TODO
        # read windows
        # read chromosome id
        # read gene names
        # read superfams
        # read orders

    def _create_sets(self, h5_file, cfg):
        """Add empty data sets to the file.

        Args:
            h5_file(h5py.File): open h5 file
            cfg(_MergeConfigSink): tuple of configuration parameters
        """

        # FUTURE resizable / unlimited datasets are possible, consider for streaming
        # OverlapData into the merge? (since it can be split across many OverlapData)
        # then one does not need to know at this point how many genes / windows there are
        # that could be more complicated but more robust...
        create_set = partial(
            h5_file.create_dataset, dtype=self.DTYPE, compression=self._compression
        )

        n_win = len(self.windows)
        n_genes = sum(1 for n in self.gene_names)
        lr_shape = (n_win, n_genes)
        i_shape = (1, n_genes)

        n_sfams = len(self.superfamily_names)
        s_left = create_set(self._S_LEFT, (n_sfams, *lr_shape))
        s_intra = create_set(self._S_INTRA, (n_sfams, *i_shape))
        s_right = create_set(self._S_RIGHT, (n_sfams, *lr_shape))
        self.superfamily = _Density(left=s_left, intra=s_intra, right=s_right)

        n_orders = len(self.order_names)
        o_left = create_set(self._O_LEFT, (n_orders, *lr_shape))
        o_intra = create_set(self._O_INTRA, (n_orders, *i_shape))
        o_right = create_set(self._O_RIGHT, (n_orders, *lr_shape))
        self.order = _Density(left=o_left, intra=o_intra, right=o_right)

    def n_updates(self, overlap):
        """Number of progress bar updates for processing density."""

        return len(self._list_density_args(overlap))

    def sum(self, overlap, gene_data, progress_bar=None):
        """Sum across the superfamily / order dimension.

        Args:
            overlap(OverlapData): container for the intermediate overlap values
            gene_data (GeneData): container for the gene data for one
            pseudomolecule
            progress_bar(?):
        """

        self._validate_chromosome(overlap)
        self._validate_windows(overlap)
        self._validate_gene_names(overlap)

        sums_ = self._list_density_args(overlap)
        n_genes = len(overlap.gene_names)
        iter_max = sum(
            len(arg.windows) * n_genes * len(arg.te_idx_name) for arg in sums_
        )
        random.shuffle(sums_)  # make progress bar update more evenly
        for args in sums_:
            self._process_sum(overlap, gene_data, args)
            if progress_bar is not None:
                progress_bar()

    def _process_sum(self, overlap, gene_data, sum_args):
        """Calculate the sum for one (left | intra | right) & (superfamily | order).

        Args:
            overlap (OverlapData): input data
            gene_data (GeneData): chromosome conatiner
            sum_args (_SummationArgs): parameters for calculations
            progress (callable): command to update progress
        """

        # N.B. loop order affects speed wrt order in data set (SEE self._create_sets)
        # N.B. also affected by memory layout, default row major (SEE numpy.array)
        # (process outer most variable in the inner-most loop)

        # FUTURE could use refactoring or a redesign, unfortunately it was rushed
        # the point of the _SummationArgs tuple was to have the same syntax to deal with
        # a) left / intra / right, and, b) superfamily / order
        # although there are so few intra calcs it might be easier to do that separately

        overlaps = sum_args.input[()]
        w_indices = [self._window_2_idx.get(w, None) for w in sum_args.windows]
        for gene_name in overlap.gene_names:
            gene_datum = gene_data.get_gene(gene_name)
            g_idx = self._gene_2_idx[gene_name]
            divisors = [sum_args.divisor_func(gene_datum, w) for w in sum_args.windows]

            te_indices, te_names = sum_args.te_idx_name
            for _i_te, _te_idx_name in enumerate(zip(te_indices, te_names)):
                te_idx, te_name = _te_idx_name
                superfam_or_order_match = sum_args.where(te_name)

                for w_idx, divisor in zip(w_indices, divisors):
                    slice_in = sum_args.slice_in(w_idx, g_idx)
                    # for one gene, one window, and all the TEs, we have overlap values
                    # select only the overlaps for the TEs that match the TE type
                    g_slice_in, w_slice_in, te_slice_in = slice_in
                    filtered_slice_in = (
                        g_slice_in,
                        w_slice_in,
                        superfam_or_order_match,
                    )
                    # sum all the entries for the gene/window at that TE type
                    overlap_sum = np.sum(
                        overlaps[g_slice_in, w_slice_in, slice(None)],
                        where=superfam_or_order_match,
                    )
                    slice_out = sum_args.slice_out(
                        window_idx=w_idx, gene_idx=g_idx, group_idx=te_idx
                    )
                    # NOTE np.divide w/ 'out' flag didn't work?
                    # NOTE must assign to the slice, rather than storing a reference to array
                    # and then assigning to it
                    sum_args.output[slice_out] = np.divide(overlap_sum, divisor)

    def _list_density_args(self, overlap):
        """List all arguments for calculating the densities."""

        superfam = self._list_sum_input_outputs(
            overlap,
            self.superfamily,
            self._config.transposons.superfamilies,
            self._config.transposons.superfamily_name_set,
            self._superfam_2_idx,
            overlap.windows,
        )
        order = self._list_sum_input_outputs(
            overlap,
            self.order,
            self._config.transposons.orders,
            self._config.transposons.order_name_set,
            self._order_2_idx,
            overlap.windows,
        )
        return superfam + order

    @classmethod
    def _list_sum_input_outputs(
        cls, overlap, density, te_group, te_set, te_idx_map, windows
    ):
        """

        Args:
            overlap(OverlapData): input overlap container
            density(_Density): density output container
            te_group(numpy.array): superfamily or order column of transposon container
            te_set(list(str)): set of superfamily or order identifiers
            windows(list(int)): window sizes

        Returns:
            iterable(_SummationArgs): arguments for calculating the sums.
        """

        # NOTE this could be (very) improved, but for now let's get a first iteration...
        arr_in = [overlap.left, overlap.intra, overlap.right]
        arr_out = [density.left, density.intra, density.right]
        win_idx_list = [
            overlap.windows,
            [None],  # intra has no window, this special case is handled in caller
            overlap.windows,
        ]
        # we need to know the index of the subset (order|family) and the name
        te_set_names = list(te_set)
        te_set_idx = [te_idx_map[t] for t in te_set_names]
        te_set_idx_name = zip([te_set_idx] * 3, [te_set_names] * 3)
        # use this lambda to make the interface for all directions the same
        slice_in = [
            overlap.left_right_slice,
            lambda w, g: overlap.intra_slice(g),
            overlap.left_right_slice,
        ]
        slice_out = [
            cls.left_right_slice,
            cls.intra_slice,
            cls.left_right_slice,
        ]
        where_out = [
            *(partial(np.equal, te_group),) * 3,
        ]  # compare ID to name

        divisor_func = [
            GeneDatum.divisor_left,
            GeneDatum.divisor_intra,
            GeneDatum.divisor_right,
        ]

        summation_args = []
        for i, o, w, te, si, so, wo, df in zip(
            arr_in,
            arr_out,
            win_idx_list,
            te_set_idx_name,
            slice_in,
            slice_out,
            where_out,
            divisor_func,
        ):
            s = _SummationArgs(
                input=i,
                output=o,
                windows=w,
                te_idx_name=te,
                slice_in=si,
                slice_out=so,
                where=wo,
                divisor_func=df,
            )
            summation_args.append(s)
        return summation_args

    def _validate_chromosome(self, overlap):
        """ValueError if the chromosome ID of the overlap data does not match this.

        The OverlapData cannot be marged into a container created for another chromosome.
        This is because the genes may be different, and more importantly the density is
            not defined between chromosomes as they are different locations.

        Args:
            overlap(OverlapData): container for the intermediate overlap values
        """

        if self.chromosome_id is None:
            raise ValueError(f"MergeData.chromosome_id (self) is None")
        elif overlap.chromosome_id is None:
            raise ValueError(f"overlap.chromosome_id is None")
        elif self.chromosome_id != overlap.chromosome_id:
            raise ValueError(
                f"""The chromosome identifier for MergeData,
                             {self.chromosome_id} does not match with the
                             chromosome identifier of overlap:
                             {overlap.chromosome_id}"""
            )

    def _validate_windows(self, overlap):
        """ValueError if the overlap windows are not in the destination.

        Args:
            overlap(OverlapData): container for the intermediate overlap values
        """

        if self.windows is None:
            raise ValueError(f"MergeData.windows (self) is None")
        elif overlap.windows is None:
            raise ValueError(f"overlap.windows is None")
        elif self.windows != overlap.windows:
            raise ValueError(
                f"""The windows for MergeData, {self.windows},
                             do not match with the windows of overlap:
                             {overlap.windows}.
                             Please delete the old overlap files, they
                             may have been calculated using a different
                             set of windows than the ones
                             provided for the current execution
                             of TE Density."""
            )

    def _validate_gene_names(self, overlap):
        """ValueError if the overlap genes are not in the destination.

        Args:
            overlap(OverlapData): container for the intermediate overlap values
        """

        if self.gene_names is None:
            raise ValueError(f"MergeData.gene_names (self) is None")
        elif overlap.gene_names is None:
            raise ValueError(f"overlap.gene_names is None")
        elif self.gene_names != overlap.gene_names:
            raise ValueError(
                f"""The gene_names of MergeData,
                             {self.gene_names}, do not match with the
                             gene_names of overlap: {overlap.gene_names}"""
            )
