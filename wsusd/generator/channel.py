import functools
import os
import shutil

import numpy as np
import scipy
from casatasks.private import sdutil
from scipy.ndimage import convolve1d

from wsusd._logging import get_logger
from wsusd.generator.util import get_spw_dd_map, get_target_spws


logger = get_logger(__name__)


@functools.lru_cache
def gauss_normalized(n, sigma):
    gauss = scipy.signal.windows.gaussian(n, sigma)
    gauss /= gauss.sum()
    return gauss


def robust_stddev(data, clipthresh=3, clipniter=3):
    """Compute robust stddev using n-sigma clipping.

    Supports both 1D arrays (single spectrum) and 2D arrays (M rows x N channels).
    Vectorized across all rows for high performance.
    """
    arr = np.asarray(data)
    if arr.ndim == 1:
        mask = np.zeros(len(arr), dtype=bool)
        for _ in range(clipniter + 1):
            counts = max(int(np.sum(~mask)), 1)
            mean = np.sum(np.where(mask, 0.0, arr)) / counts
            var = np.sum(np.where(mask, 0.0, np.abs(arr - mean) ** 2)) / counts
            stddev = float(np.sqrt(var))
            thresh = stddev * clipthresh
            mask |= np.abs(arr) > thresh
        return stddev

    mask = np.zeros_like(arr, dtype=bool)
    for _ in range(clipniter + 1):
        counts = np.maximum(np.sum(~mask, axis=1, keepdims=True), 1)
        mean = np.sum(np.where(mask, 0.0, arr), axis=1, keepdims=True) / counts
        var = np.sum(np.where(mask, 0.0, np.abs(arr - mean) ** 2), axis=1, keepdims=True) / counts
        stddev = np.sqrt(var)
        thresh = stddev * clipthresh
        mask |= np.abs(arr) > thresh
    return stddev


def interpolate_data(data_in, cf_in, cf_out, sigma=2, robust=True):
    """Interpolates data onto a new frequency grid with added noise.

    Performs Gaussian smoothing, interpolation, and re-addition of noise based on
    residuals from the smoothing process.

    Args:
        data_in: Input data array with shape (npol, nchan, nrow)
        cf_in: Input channel frequencies.
        cf_out: Output channel frequencies.
        sigma: Standard deviation for the Gaussian kernel.
        robust: If True, use robust stddev estimation for noise addition (slower).

    Returns:
        Interpolated data array with shape (npol, nchan_out, nrow).
    """
    # data.shape should be (npol, nchan, nrow) in F-order in memory per casatools output convention.
    assert len(data_in.shape) == 3
    assert len(cf_in.shape) == 1
    assert len(cf_out.shape) == 1

    npol, nchan_in, nrow = data_in.shape
    nchan_out = len(cf_out)

    n = min(nchan_in, round(sigma * 10))
    # Ensure an odd number of channels for a symmetric zero-phase Gaussian kernel
    if n % 2 == 0 and n + 1 <= nchan_in:
        n += 1
    gauss = gauss_normalized(n, sigma)

    # Reshape for vectorized convolution: (npol * nrow, nchan)
    data_flat = data_in.transpose(0, 2, 1).reshape(-1, nchan_in)

    # Vectorized smoothing using scipy.ndimage
    # Set origin to align with symmetric filtering and eliminate any 1-sample shifts
    origin = -1 if n % 2 == 0 else 0
    smoothed_flat = convolve1d(
        data_flat, gauss, axis=1, mode='constant', cval=0.0, origin=origin
    )

    # Create interpolator for all rows
    interpolator = scipy.interpolate.interp1d(
        cf_in, smoothed_flat, axis=1, bounds_error=False, fill_value=(smoothed_flat[:, 0], smoothed_flat[:, -1])
    )
    interpolated_flat = interpolator(cf_out)

    # Add noise vectorized
    diff_flat = data_flat - smoothed_flat
    if robust:
        noise_std = robust_stddev(diff_flat)
        noise_std_native = np.std(diff_flat, axis=1, keepdims=True)
        logger.debug(
            'native std %s robust std %s (channel-wise noise - median over chunk)',
            np.median(noise_std_native),
            np.median(noise_std),
        )
    else:
        noise_std = np.std(diff_flat, axis=1, keepdims=True)
    noise = np.random.default_rng().normal(0, noise_std, interpolated_flat.shape)
    corrupted_flat = interpolated_flat + noise

    # Reshape back to (npol, nchan_out, nrow)
    return corrupted_flat.reshape(npol, nrow, nchan_out).transpose(0, 2, 1)


def interpolate_bool(data_in, cf_in, cf_out):
    """Nearest-neighbor channel interpolation for boolean flag arrays.

    Vectorized across all rows and polarizations simultaneously.
    """
    assert len(data_in.shape) == 3
    assert len(cf_in.shape) == 1
    assert len(cf_out.shape) == 1

    indices = np.searchsorted(cf_in, cf_out)
    indices = np.clip(indices, 0, len(cf_in) - 1)
    left_indices = np.maximum(indices - 1, 0)
    dist_right = np.abs(cf_in[indices] - cf_out)
    dist_left = np.abs(cf_in[left_indices] - cf_out)
    nearest_idx = np.where(dist_left <= dist_right, left_indices, indices)

    return data_in[:, nearest_idx, :]


class TableUpdater:
    @property
    def columns(self):
        return []

    def taql(self, spw: int) -> str:
        return ''

    def __init__(self, vis: str, target_spws: list[int], table_name, **kwargs):
        self.vis = vis
        self.target_spws = [int(spw) for spw in target_spws]
        self.table_name = table_name
        for attr in ['chan_factor', 'freq_in', 'freq_out']:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])

        self.data_in: dict = {}
        self.data_out: dict = {}

    def read(self):
        self.data_in = {}
        for spw in self.target_spws:
            taql = self.taql(spw)
            assert len(taql) > 0
            with sdutil.table_selector(
                self.table_name, taql, nomodify=True
            ) as tb:
                if tb.nrows() == 0:
                    continue

                existing_columns = tb.colnames()
                _data_in = dict(
                    (col, tb.getcol(col)) for col in self.columns
                    if col in existing_columns and tb.iscelldefined(col, 0)
                )

            self.data_in[spw] = _data_in

    def update(self):
        pass

    def flush(self):
        if not isinstance(self.data_out, dict):
            return

        for spw in self.target_spws:
            taql = self.taql(spw)
            assert len(taql) > 0
            with sdutil.table_selector(
                self.table_name, taql, nomodify=False
            ) as tb:
                _data_out = self.data_out[spw]
                for col, val in _data_out.items():
                    tb.putcol(col, val)


class SpectralWindowUpdater(TableUpdater):
    @property
    def columns(self):
        return [
            'CHAN_FREQ', 'CHAN_WIDTH', 'NUM_CHAN',
            'EFFECTIVE_BW', 'RESOLUTION'
        ]

    def taql(self, spw: int) -> str:
        # inside TaQL, ROWNUMBER returns 1-based
        # row number while spw is 0-based
        return f'ROWNUMBER() == {spw + 1}'

    def __init__(self, vis: str, target_spws: list[int], chan_factor: float):
        table_name = os.path.join(vis, 'SPECTRAL_WINDOW')
        super().__init__(vis, target_spws, table_name)
        self.chan_factor = chan_factor

    def get_chan_freq_in(self):
        if isinstance(self.data_in, dict):
            return dict(
                (spw, data['CHAN_FREQ'][:, 0]) for spw, data
                in self.data_in.items()
            )
        else:
            raise ValueError('data_in is not initialized properly.')

    def get_chan_freq_out(self):
        if isinstance(self.data_out, dict):
            return dict(
                (spw, data['CHAN_FREQ'][:, 0]) for spw, data
                in self.data_out.items()
            )
        else:
            raise ValueError('data_out is not initialized properly.')

    def _update_num_chan(self, spw: int, nchan: np.ndarray) -> np.ndarray:
        nchan_new = np.rint(nchan * self.chan_factor).astype(int)
        logger.info(
            f'spw {spw}: nchan(in) {nchan[0]}, nchan(out) {nchan_new[0]}'
        )
        return nchan_new

    def _update_chan_width(self, spw, chan_width):
        chan_freq = self.data_in[spw]['CHAN_FREQ']
        nchan_new = self.data_out[spw]['NUM_CHAN']

        start_chan = chan_freq[0][0] - chan_width[-1][0] / 2
        end_chan = chan_freq[-1][0] + chan_width[-1][0] / 2
        bandwidth = end_chan - start_chan

        cw_new = np.zeros(nchan_new, dtype=chan_width.dtype) \
            + (bandwidth / nchan_new)
        logger.debug(
            f'spw {spw}: chan_width(in) {chan_width[0][0]} '
            f'chan_width(out) {cw_new[0]}'
        )
        return cw_new[:, np.newaxis]

    def _update_chan_freq(self, spw, chan_freq):
        chan_width = self.data_in[spw]['CHAN_WIDTH']
        cw_new = self.data_out[spw]['CHAN_WIDTH']
        start_chan = chan_freq[0][0] - chan_width[-1][0] / 2
        end_chan = chan_freq[-1][0] + chan_width[-1][0] / 2
        _start = start_chan + cw_new[0][0] / 2
        _end = end_chan
        _step = cw_new[0][0]
        logger.debug(f'_start {_start} _end {_end} _step {_step}')

        cf_new = np.arange(_start, _end, _step, dtype=chan_freq.dtype)
        logger.debug(
            f'spw {spw}: chan_freq(in) {chan_freq[0][0]} '
            f'chan_freq(out) {cf_new[0]}')
        return cf_new[:, np.newaxis]

    def _scale_by_chan_width(self, spw, data):
        cw = self.data_in[spw]['CHAN_WIDTH']
        cw_new = self.data_out[spw]['CHAN_WIDTH']
        data_new = np.abs(cw_new * data[0][0] / cw[0][0])
        return data_new

    def _no_scale(self, spw, data):
        cw_new = self.data_out[spw]['CHAN_WIDTH']
        data_new = np.zeros_like(cw_new) + abs(data[0][0])
        return data_new

    def _update_col(self, column, update_func):
        for spw in self.target_spws:
            _data_in = self.data_in[spw][column]
            _data_out = update_func(spw, _data_in)
            self.data_out[spw][column] = _data_out

    def update(self):
        self.data_out = dict((spw, {}) for spw in self.target_spws)

        # NUM_CHAN
        self._update_col('NUM_CHAN', self._update_num_chan)

        # CHAN_WIDTH
        self._update_col('CHAN_WIDTH', self._update_chan_width)

        # CHAN_FREQ
        self._update_col('CHAN_FREQ', self._update_chan_freq)

        # EFFECTIVE_BW
        # effective noise bandwidth is proportional to channel width
        self._update_col('EFFECTIVE_BW', self._scale_by_chan_width)
        # self._update_col('EFFECTIVE_BW', self._no_scale)

        # RESOLUTION
        self._update_col('RESOLUTION', self._scale_by_chan_width)


class SyscalUpdator(TableUpdater):
    @property
    def columns(self):
        return [
            'TCAL_SPECTRUM', 'TRX_SPECTRUM', 'TSKY_SPECTRUM',
            'TSYS_SPECTRUM', 'TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM'
        ]

    def taql(self, spw):
        return f'SPECTRAL_WINDOW_ID == {spw}'

    def __init__(
            self, vis: str, spw: list[int], freq_in: dict, freq_out: dict
    ):
        table_name = os.path.join(vis, 'SYSCAL')
        super().__init__(vis, spw, table_name)
        self.freq_in = freq_in
        self.freq_out = freq_out

    def update(self):
        self.data_out = dict((spw, {}) for spw in self.target_spws)

        for spw, _data_in in self.data_in.items():
            _freq_in = self.freq_in[spw]
            _freq_out = self.freq_out[spw]
            for col, arr in _data_in.items():
                logger.info(f'spw {spw}: updating column {col}')
                # scale temperature data with 1 / sqrt(chan_width)
                # _cw_in = abs(_freq_in[1] - _freq_in[0])
                # _cw_out = abs(_freq_out[1] - _freq_out[0])
                # scale_factor = np.sqrt(_cw_in / _cw_out)
                # logger.info(
                #     f'spw {spw}: temperature scaling factor is '
                #     f'{scale_factor}'
                # )
                scale_factor = 1.0
                self.data_out[spw][col] = interpolate_data(
                    arr, _freq_in, _freq_out
                ) * scale_factor


def copy_table_structure(vis, outputvis):
    logger.debug(f'copying {vis} into {outputvis}')
    with sdutil.table_manager(vis) as tb:
        tout = tb.copy(outputvis, norows=True)
        tout.close()

    logger.debug('done copying')


def copy_subtable_rows(vis, outputvis):
    with sdutil.table_manager(vis) as tb:
        table_names = filter(
            lambda x: isinstance(x[1], str) and x[1].startswith('Table: '),
            tb.getkeywords().items()
        )

    for name, path in table_names:
        logger.debug(f'copying {name} rows')
        src = path[7:]
        dst = os.path.join(outputvis, name)
        with sdutil.table_manager(src) as tb:
            tb.copyrows(dst)

        logger.debug(f'done {name}')


def copy_main_columns(vis, outputvis, ignore):
    with sdutil.table_manager(outputvis, nomodify=False) as tb_out:
        nrow_out = tb_out.nrows()
        with sdutil.table_manager(vis) as tb_in:
            nrow = tb_in.nrows()
            tb_out.addrows(nrow - nrow_out)
            for col in tb_out.colnames():
                logger.debug(col)
                if col in ignore or not tb_in.iscelldefined(col, 0):
                    continue

                logger.debug(f'copy column {col}')
                data = tb_in.getcol(col)
                tb_out.putcol(col, data)
                logger.debug(f'done {col}')


def rename_table(src, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)

    os.rename(src, dst)


class ChunkInterpolator:
    def __init__(self, cf_in, cf_out):
        self.cf_in = cf_in
        self.cf_out = cf_out

        # Precompute nearest channel indices for FLAG interpolation
        indices = np.searchsorted(cf_in, cf_out)
        indices = np.clip(indices, 0, len(cf_in) - 1)
        left_indices = np.maximum(indices - 1, 0)
        dist_right = np.abs(cf_in[indices] - cf_out)
        dist_left = np.abs(cf_in[left_indices] - cf_out)
        self.nearest_flag_idx = np.where(dist_left <= dist_right, left_indices, indices)

    def __scale_data(self, data_in, factor):
        # shape of data_in should be (npol, nchan, nrow)
        assert len(data_in.shape) == 3
        assert data_in.shape[1] == len(self.cf_in)

        npol, _, nrow = data_in.shape
        nchan = len(self.cf_out)
        data_out = np.zeros((npol, nchan, nrow), dtype=data_in.dtype)
        data_out[:] = data_in[0, 0, 0] * factor
        return data_out

    def __update_weight_spectrum(self, data_in):
        # weight is proportional to channel width
        cw_in = abs(self.cf_in[1] - self.cf_in[0])
        cw_out = abs(self.cf_out[1] - self.cf_out[0])
        factor = cw_out / cw_in
        return self.__scale_data(data_in, factor)

    def __update_sigma_spectrum(self, data_in):
        # sigma is proportional to 1 / sqrt(channel width)
        cw_in = abs(self.cf_in[1] - self.cf_in[0])
        cw_out = abs(self.cf_out[1] - self.cf_out[0])
        factor = 1 / np.sqrt(cw_out / cw_in)
        return self.__scale_data(data_in, factor)

    def __call__(self, chunk):
        chunk_start, nrow_chunk, data_in = chunk

        data_out = {}

        # update data column (FLOAT_DATA/DATA)
        for column in ['FLOAT_DATA', 'DATA']:
            if column in data_in:
                data_out[column] = interpolate_data(
                    data_in[column], self.cf_in, self.cf_out
                )
                break

        # update FLAG
        column = 'FLAG'
        if column in data_in:
            data_out[column] = data_in[column][:, self.nearest_flag_idx, :]

        # update WEIGHT_SPECTRUM
        column = 'WEIGHT_SPECTRUM'
        if column in data_in:
            data_out[column] = self.__update_weight_spectrum(data_in[column])

        # update_SIGMA_SPECTRUM
        column = 'SIGMA_SPECTRUM'
        if column in data_in:
            data_out[column] = self.__update_sigma_spectrum(data_in[column])

        return chunk_start, nrow_chunk, data_out


DEFAULT_NROW_CHUNK = 10000


class MainUpdater(TableUpdater):
    @functools.lru_cache(1)
    def __columns(self):
        columns = []
        with sdutil.table_manager(self.vis) as tb:
            colnames = tb.colnames()
            colname = 'FLOAT_DATA'
            if colname in colnames:
                columns.append(colname)
            else:
                columns.append('DATA')

            columns.append('FLAG')

            for colname in ['SIGMA_SPECTRUM', 'WEIGHT_SPECTRUM']:
                if colname in colnames and tb.iscelldefined(colname, 0):
                    columns.append(colname)

        return columns

    @property
    def columns(self):
        return self.__columns()

    def taql(self, ddid):
        return f'DATA_DESC_ID == {ddid}'

    def __init__(
            self, vis: str, target_spws: list, freq_in: dict, freq_out: dict
    ):
        table_name = vis
        super().__init__(vis, target_spws, table_name)
        self.freq_in = freq_in
        self.freq_out = freq_out

        self.dd_entries = []
        self.ddid = {}
        spw_ddid_map = get_spw_dd_map(self.vis)
        for spw, ddids in sorted(spw_ddid_map.items()):
            if len(ddids) > 0 and spw not in self.ddid:
                self.ddid[spw] = ddids[0]
            for ddid in ddids:
                self.dd_entries.append((spw, ddid))
        self.all_spws = sorted(self.ddid.keys())

        vis_dir = os.path.dirname(vis)
        vis_base = os.path.basename(vis)
        self.tmp_vis = os.path.join(vis_dir, f'genwsusd.{vis_base}.tmp')
        self.backup_vis = os.path.join(vis_dir, f'genwsusd.{vis_base}.bak')

    def _read_main(self, spw, ddid):
        nrow_chunk_default = int(
            os.environ.get('WSUSD_NROW_CHUNK', DEFAULT_NROW_CHUNK)
        )

        taql = self.taql(ddid)
        with sdutil.table_selector(self.vis, taql) as tb:
            nrow = tb.nrows()
            nchunk = nrow // nrow_chunk_default
            nmod = nrow % nrow_chunk_default
            chunk_list = [nrow_chunk_default] * nchunk
            if nmod > 0:
                chunk_list.append(nmod)

            chunk_start = 0
            for i, nrow_chunk in enumerate(chunk_list):
                logger.info(
                    'spw %d (ddid %d): start reading chunk %d', spw, ddid, i
                )

                data_in = dict(
                    (name, tb.getcol(name, chunk_start, nrow_chunk))
                    for name in self.columns
                )

                yield chunk_start, nrow_chunk, data_in

                logger.info(
                    'spw %d (ddid %d): done reading chunk %d', spw, ddid, i
                )
                chunk_start += nrow_chunk

    def _get_chunk_updater(self, spw):
        if spw in self.target_spws:
            logger.info('spw %d: update chunk', spw)
            return ChunkInterpolator(self.freq_in[spw], self.freq_out[spw])
        else:
            logger.info('spw %d: leave input chunk as it is', spw)
            return lambda x: x

    def read(self):
        # here, copy input MS to temporary MS
        copy_table_structure(self.vis, self.tmp_vis)
        copy_subtable_rows(self.vis, self.tmp_vis)
        copy_main_columns(self.vis, self.tmp_vis, ignore=self.columns)

        # create generators for lazy read
        self.read_generators = [
            self._read_main(spw, ddid) for spw, ddid in self.dd_entries
        ]

    def update(self):
        self.update_generators = [
            map(self._get_chunk_updater(spw), chunk)
            for (spw, ddid), chunk in zip(self.dd_entries, self.read_generators)
        ]

    def flush(self):
        # flush to the disk
        try:
            for (spw, ddid), update_gen in zip(self.dd_entries, self.update_generators):
                logger.debug('spw %d (ddid %d): generator %s', spw, ddid, update_gen)
                taql = self.taql(ddid)
                with sdutil.table_selector(
                    self.tmp_vis, taql, nomodify=False
                ) as tb:
                    for i, (chunk_start, nrow_chunk, data_out) in enumerate(update_gen):
                        logger.info(
                            'spw %d (ddid %d): start writing chunk %d', spw, ddid, i
                        )
                        for column, chunk in data_out.items():
                            tb.putcol(column, chunk, chunk_start, nrow_chunk)

                        logger.info(
                            'spw %d (ddid %d): done writing chunk %d', spw, ddid, i
                        )

            # finalization
            rename_table(self.vis, self.backup_vis)
            try:
                rename_table(self.tmp_vis, self.vis)
            except Exception as e:
                # If renaming tmp_vis fails, restore original MS from backup_vis
                logger.error(
                    'Failed to rename %s to %s: %s. Attempting to restore original MS.',
                    self.tmp_vis, self.vis, e
                )
                if os.path.exists(self.backup_vis) and not os.path.exists(self.vis):
                    rename_table(self.backup_vis, self.vis)
                raise
        except Exception:
            if os.path.exists(self.tmp_vis):
                shutil.rmtree(self.tmp_vis)
            logger.error(
                'Error during processing or table swapping. Original MS preserved if backup restoration succeeded.'
            )
            raise
        finally:
            # Only remove backup if the destination MS exists safely on disk
            if os.path.exists(self.backup_vis) and os.path.exists(self.vis):
                shutil.rmtree(self.backup_vis)


class WSUChannelExpander:
    def __init__(self, vis: str, chan_factor: float):
        self.vis = vis
        self.chan_factor = chan_factor

        self.science_spws, self.atm_spws = get_target_spws(self.vis)
        self.target_spws = self.science_spws + self.atm_spws

    def expand(self, dry_run: bool = False):
        # process SPECTRAL_WINDOW table
        spw_updater = SpectralWindowUpdater(
            self.vis, self.target_spws, self.chan_factor
        )
        spw_updater.read()

        # dry run mode: just report updated NUM_CHAN values
        if dry_run:
            for spw, nchan in spw_updater.data_in.items():
                _ = spw_updater._update_num_chan(spw, nchan['NUM_CHAN'])
            return

        spw_updater.update()
        spw_updater.flush()

        cf_in = spw_updater.get_chan_freq_in()
        cf_out = spw_updater.get_chan_freq_out()

        del spw_updater

        # process SYSCAL table
        # - depends on SPECTRAL_WINDOW information
        syscal_updater = SyscalUpdator(self.vis, self.atm_spws, cf_in, cf_out)
        syscal_updater.read()
        syscal_updater.update()
        syscal_updater.flush()

        del syscal_updater

        # process MAIN table
        # - depends on SPECTRAL_WINDOW information
        main_updater = MainUpdater(self.vis, self.target_spws, cf_in, cf_out)
        main_updater.read()
        main_updater.update()
        main_updater.flush()
