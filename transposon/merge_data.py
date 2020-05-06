#!/usr/bin/env python

"""
Sums and contains OverlapData with respect to the genes / transposons.


"""

__author__ = "Michael Teresi"

from collections import namedtuple
import logging
import os
from functools import partial

import h5py
import numpy as np

import transposon


_MergeConfigSink = namedtuple(
    '_MergeConfigSink',
    ['transposons', 'gene_names', 'windows', 'filepath', 'ram_bytes']
)
_MergeConfigSource = namedtuple('_MergeConfigSource', ['filepath'])
_Density = namedtuple('_Density', ['left', 'intra', 'right'])
_SummationArgs = namedtuple(
    '_SummationArgs',
    ['input', 'output', 'windows', 'te_idx_name', 'slice_in', 'slice_out', 'where'])


class MergeData():
    """Contains the sum and density of the overlap wrt superfamilies and orders.

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

    DTYPE = np.float32  # MAGIC NUMBER experimental, depends on your data
    _GENE_NAMES = 'GENE_NAMES'
    _WINDOWS = 'WINDOWS'
    _CHROME_ID = 'CHROMOSOME_ID'
    _SUPERFAMILY_NAMES = 'SUPERFAMILY_NAMES'
    _ORDER_NAMES = 'ORDER_NAMES'
    _RHO_SUPERFAMILY_DATASET_KEY = 'RHO_SUPERFAMILIES'
    _S_LEFT = _RHO_SUPERFAMILY_DATASET_KEY + '_LEFT'
    _S_RIGHT = _RHO_SUPERFAMILY_DATASET_KEY + '_RIGHT'
    _S_INTRA = _RHO_SUPERFAMILY_DATASET_KEY + '_INTRA'
    _RHO_ORDER_DATASET_KEY = 'RHO_ORDERS'
    _O_LEFT = _RHO_ORDER_DATASET_KEY + '_LEFT'
    _O_RIGHT = _RHO_ORDER_DATASET_KEY + '_RIGHT'
    _O_INTRA = _RHO_ORDER_DATASET_KEY + '_INTRA'

    def __init__(self, configuration, compression='lzf', logger=None):
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

        LOSSLESS_FILTERS = ['gzip', 'lzf', 'szip']  # SEE h5py datasets
        if filter not in LOSSLESS_FILTERS:
            msg = "'%s' compression invalid, options are: %s" % (filter, LOSSLESS_FILTERS)
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
    def from_param(cls, transposon_data, gene_names, windows,
                   output_dir, ram=2, logger=None):
        """Writable sink for a new file.

        Args:
            transposon_data (TransposonData): transposon container.
            gene_names (iterable(str)): gene names to process.
            windows (iterable(int)): window inputs to process (not the same as windows).
            output_dir (str): directory to output merge data files.
            rame (int): gigabytes RAM to cache during processing.
            logger (logging.Logger): logging instance for messages.
        """

        logger = logger or logging.getLogger(__name__)
        bytes2gigabytes = 1024. ** 3
        ram_bytes = int(ram * bytes2gigabytes)
        transposon.check_ram(ram_bytes, logger)
        chrome = str(transposon_data.chromosome_unique_id)
        genome = str(transposon_data.genome_id)
        filename = genome + '_' + chrome + '.h5'
        filepath = os.path.join(output_dir, filename)
        config = _MergeConfigSink(
            transposons=transposon_data,
            gene_names=gene_names,
            windows=windows,
            filepath=filepath,
            ram_bytes=ram_bytes
        )
        return cls(config, logger=logger)

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Read only source for an existing file."""

        raise NotImplementedError()

    @classmethod
    def left_right_slice(cls, group_idx=None, gene_idx=None, window_idx=None):
        """Slice to a left || right density to a superfamily|order / window, gene."""

        if group_idx is None or gene_idx is None or window_idx is None:
            kwargs = {'group_idx': group_idx,
                      'gene_idx': gene_idx,
                      'window_idx': window_idx}
            raise ValueError("cannot slice with input of None: %s" % kwargs)

        # SEE numpy basic indexing
        return (group_idx, window_idx, slice(gene_idx, gene_idx+1, 1))

    @classmethod
    def intra_slice(cls, group_idx=None, gene_idx=None, window_idx=None):
        """Slice for an intra densitiy to a superfamily|order / window, gene."""

        if window_idx is not None:
            raise ValueError("intra window must be None but is %s" % window_idx)
        return cls.left_right_slice(group_idx=group_idx, gene_idx=gene_idx, window_idx=0)

    def __enter__(self):
        """Context manager begin."""

        self._open_dispatcher()
        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        """Context manager end."""

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
            raise TypeError("expecting {} or {} but got {}".format(
                type(_MergeConfigSink), type(_MergeConfigSource), type(self._config)))

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
        self._h5_file = h5py.File(cfg.filepath, 'w', rdcc_nbytes=cfg.ram_bytes)
        self._create_sets(self._h5_file, cfg)
        transposon.write_vlen_str_h5py(self._h5_file, self.windows, self._WINDOWS)
        transposon.write_vlen_str_h5py(self._h5_file, self.gene_names, self._GENE_NAMES)
        transposon.write_vlen_str_h5py(self._h5_file, self.chromosome_id, self._CHROME_ID)
        self._window_2_idx = {w: i for i, w in enumerate(self.windows)}
        self._gene_2_idx = {g: i for i, g in enumerate(self.gene_names)}

    def _open_existing_file(self, cfg):

        self._h5_file = h5py.File(cfg.filepath, 'r')
        transposon.read_vlen_str_h5py(self._h5_file, )
        # read windows
        # read chromosome id
        # read gene names
        # read superfams
        # read orders

    def _create_sets(self, h5_file, cfg):
        """Add data sets to the file.

        Args:
            h5_file(h5py.File): open h5 file
            cfg(_MergeConfigSink): tuple of configuration parameters
        """

        # FUTURE resizable / unlimited datasets are possible, consider for streaming
        # OverlapData into the merge? (since it can be split across many OverlapData)
        # then one does not need to know at this point how many genes / windows there are
        # that could be more complicated but more robust...
        create_set = partial(h5_file.create_dataset,
                             dtype=self.DTYPE,
                             compression=self._compression)

        n_win = len(self.windows)
        n_genes = sum(1 for n in self.gene_names)
        lr_shape = (n_win, n_genes)
        i_shape = (1, n_genes)

        n_sfams = len(self.superfamily_names)
        s_left = create_set(self._S_LEFT, (n_sfams, *lr_shape))
        s_intra = create_set(self._S_INTRA, (n_sfams, *i_shape))
        s_right = create_set(self._S_RIGHT, (n_sfams, *lr_shape))
        self.superfamily = _Density(left=s_left[:], intra=s_intra[:], right=s_right[:])

        n_orders = len(self.order_names)
        o_left = create_set(self._O_LEFT, (n_orders, *lr_shape))
        o_intra = create_set(self._O_INTRA, (n_orders, *i_shape))
        o_right = create_set(self._O_RIGHT, (n_orders, *lr_shape))
        self.order = _Density(left=o_left[:], intra=o_intra[:], right=o_right[:])

    def sum(self, overlap, progress_bar):
        """Sum across the superfamily / order dimension."""

        self._validate_chromosome(overlap)
        self._validate_windows(overlap)
        self._validate_gene_names(overlap)

        sums_ = self._list_sum_args(overlap)
        n_genes = len(overlap.gene_names)
        iter_max = sum(len(arg.windows) * n_genes * len(arg.te_idx_name) for arg in sums_)
        # TODO set progress bar
        for args in sums_:
            self._process_sum(overlap, args, progress_bar)

    def _process_sum(self, overlap, sum_args, progress):
        """Calculate the sum for one (left | intra | right) & (superfamily | order).

        Args:
            overlap (OverlapData): input data
            sum_args (_SummationArgs): container for arguments wrt sum
            progress (tqdm.std.tqdm): progress bar
        """

        # N.B. loop order affects speed wrt order in data set (SEE self._create_sets)
        # N.B. also affected by memory layout, default row major (SEE numpy.array)

        # FUTURE could use refactoring or a redesign, unfortunately it was rushed
        # the point of the _SummationArgs tuple was to have the same syntax to deal with
        # a) left / intra / right, and, b) superfamily / order
        # although there are so few intra calcs it might be easier to do that separately

        for te_ in zip(*sum_args.te_idx_name):  # for every superfam | order
            te_idx, te_name = te_  # the supefamily / order index and string
            for window in sum_args.windows:  # for every window
                w_idx = self._window_2_idx.get(window, None)
                for gene in overlap.gene_names:  # for every gene
                    g_idx = self._gene_2_idx[gene]
                    slice_out = sum_args.slice_out(
                        window_idx=w_idx, gene_idx=g_idx, group_idx=te_idx)
                    slice_in = sum_args.slice_in(w_idx, g_idx)
                    # find which genes match the superfam | order, and sum those
                    superfam_or_order_match = sum_args.where(te_name)
                    slice_in_only_te = (superfam_or_order_match, *slice_in[1:])
                    result = np.sum(sum_args.input[slice_in_only_te])
                    sum_args.output[slice_out] = result

    def _list_sum_args(self, overlap):

        superfam = self._list_sum_input_outputs(
            overlap,
            self.superfamily,
            self._config.transposons.superfamilies,
            self._config.transposons.superfamily_name_set,
            self._superfam_2_idx,
            overlap.windows
        )
        order = self._list_sum_input_outputs(
            overlap,
            self.order,
            self._config.transposons.orders,
            self._config.transposons.order_name_set,
            self._order_2_idx,
            overlap.windows
        )
        return superfam + order

    @classmethod
    def _list_sum_input_outputs(cls, overlap, density, te_group,
                                te_set, te_idx_map, windows):
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
            overlap.windows
        ]
        # we need to know the index of the subset (order|family) and the name
        te_set_names = list(te_set)
        te_set_idx = [te_idx_map[t] for t in te_set_names]
        te_set_idx_name = zip([te_set_idx]*3, [te_set_names]*3)
        # use this lambda to make the interface for all directions the same
        slice_in = [
            overlap.left_right_slice,
            lambda w, g : overlap.intra_slice(g),
            overlap.left_right_slice,
        ]
        slice_out = [
            cls.left_right_slice,
            cls.intra_slice,
            cls.left_right_slice,
        ]
        where_out = [*(partial(np.equal, te_group), )*3, ]  # compare ID to name
        summation_args = []
        for i, o, w, te, si, so, wo in \
            zip(arr_in, arr_out, win_idx_list, te_set_idx_name,
                slice_in, slice_out, where_out):
            s = _SummationArgs(
                input=i,
                output=o,
                windows=w,
                te_idx_name=te,
                slice_in=si,
                slice_out=so,
                where=wo
            )
            summation_args.append(s)
        return summation_args

    def _validate_chromosome(self, overlap):
        """ValueError if the chromosome ID of the overlap data does not match this.

        The OverlapData cannot be marged into a container created for another chromosome.
        This is because the genes may be different, and more importantly the density is
            not defined between chromosomes as they are different locations.
        """
        if self.chromosome_id is None:
            raise ValueError(f"MergeData.chromosome_id (self) is None")
        elif overlap.chromosome_id is None:
            raise ValueError(f"overlap.chromosome_id is None")
        elif self.chromosome_id != overlap.chromosome_id:
            raise ValueError(f"""The chromosome identifier for MergeData,
                             {self.chromosome_id} does not match with the
                             chromosome identifier of overlap:
                             {overlap.chromosome_id}""")

    def _validate_windows(self, overlap):
        """ValueError if the overlap windows are not in the destination."""
        if self.windows is None:
            raise ValueError(f"MergeData.windows (self) is None")
        elif overlap.windows is None:
            raise ValueError(f"overlap.windows is None")
        elif self.windows != overlap.windows:
            raise ValueError(f"""The windows for MergeData, {self.windows},
                             do not match with the windows of overlap:
                             {overlap.windows}""")

    def _validate_gene_names(self, overlap):
        """ValueError if the overlap genes are not in the destination."""
        if self.gene_names is None:
            raise ValueError(f"MergeData.gene_names (self) is None")
        elif overlap.gene_names is None:
            raise ValueError(f"overlap.gene_names is None")
        elif self.gene_names != overlap.gene_names:
            raise ValueError(f"""The gene_names of MergeData,
                             {self.gene_names}, do not match with the
                             gene_names of overlap: {overlap.gene_names}""")
