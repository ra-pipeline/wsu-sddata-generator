import os

import numpy as np

from casatasks.private import sdutil

from wsusd._logging import get_logger
from wsusd.generator.util import get_spw_dd_map, get_target_spws

logger = get_logger(__name__)


class WSUSpwExpander:
    def __init__(self, vis: str, spw_factor: int):
        self.vis = vis
        self.spw_factor = spw_factor

        self.science_spws, self.atm_spws = get_target_spws(self.vis)
        self.target_spws = sorted(set(self.science_spws + self.atm_spws))
        logger.info(f'target spws: {self.target_spws}')

        self.spw_dd_map = get_spw_dd_map(self.vis)

    def expand(self, dry_run: bool = False):
        num_duplication = self.spw_factor - 1
        logger.info(
            f'number of output spws: {len(self.target_spws) * self.spw_factor}'
        )
        logger.info(f'duplicate spws {num_duplication} times')

        if dry_run:
            return

        for i in range(num_duplication):
            logger.info(f'start duplication cycle {i}')
            self.cycle_id = i
            self._duplicate()
            logger.info(f'done duplication cycle {i}')

    def __generate_spw_name(self, base_name, extra_digit):
        separator = '#'
        elements = base_name.split(separator)
        # append extra digit to the first element
        elements[0] += f'{extra_digit:0>2d}'
        return separator.join(elements)

    def __find_dst_spw_id(self, tb):
        taql = 'NAME == pattern("WVR#Antenna*")'
        tsel = tb.query(taql)
        rows = tsel.rownumbers()
        tsel.close()

        if len(rows) == 0:
            return -1
        else:
            return rows[0]

    def _expand_spectral_window(self):
        table_name = os.path.join(self.vis, 'SPECTRAL_WINDOW')
        self.extra_spw = []
        with sdutil.table_manager(table_name, nomodify=False) as tb:
            for base_spw in self.target_spws:
                startrow = self.__find_dst_spw_id(tb)
                if startrow < 0:
                    new_spw = int(tb.nrows())
                else:
                    new_spw = int(startrow)
                logger.info(
                    f'duplicating spw {base_spw}: spw {new_spw} will be added'
                )
                tb.copyrows(table_name, base_spw, startrowout=startrow, nrow=1)
                base_spw_name = tb.getcell('NAME', base_spw)
                new_spw_name = self.__generate_spw_name(
                    base_spw_name, self.cycle_id
                )
                logger.info(f'base spw name: "{base_spw_name}"')
                logger.info(f' new spw_name: "{new_spw_name}"')
                tb.putcell('NAME', new_spw, new_spw_name)
                self.extra_spw.append(new_spw)

    def __copy_selected_rows(self, tb, taql):
        table_name = tb.name()

        selected = tb.query(taql)
        try:
            if selected.nrows() == 0:
                return

            selected.copyrows(table_name)
        finally:
            selected.close()

    def __remove_preexisting_rows(self, tb, spw_id):
        taql = f'SPECTRAL_WINDOW_ID == {spw_id}'
        tsel = tb.query(taql)
        rows = tsel.rownumbers()
        tsel.close()

        if len(rows) > 0:
            tb.removerows(rows)

    def __expand_subtable(self, subtable_name):
        full_table_name = os.path.join(self.vis, subtable_name)
        with sdutil.table_manager(full_table_name, nomodify=False) as tb:
            nrow_before = tb.nrows()

            for base_spw, new_spw in zip(self.target_spws, self.extra_spw):
                logger.info(f'duplicating {subtable_name} rows for '
                            f'spw {base_spw}: '
                            f'spw {new_spw} will be assigned')

                self.__remove_preexisting_rows(tb, new_spw)

                startrow = tb.nrows()
                taql = f'SPECTRAL_WINDOW_ID == {base_spw}'
                self.__copy_selected_rows(tb, taql)
                nrow = tb.nrows() - startrow
                spwcol = np.zeros(nrow, dtype=int) + new_spw
                tb.putcol('SPECTRAL_WINDOW_ID', spwcol, startrow, nrow)

            nrow_after = tb.nrows()

        return nrow_before, nrow_after

    def _expand_data_description(self):
        full_table_name = os.path.join(self.vis, 'DATA_DESCRIPTION')
        self.dd_pairs: list[tuple[int, int, int, int]] = []
        with sdutil.table_manager(full_table_name, nomodify=False) as tb:
            for base_spw, new_spw in zip(self.target_spws, self.extra_spw):
                logger.info(
                    'duplicating DATA_DESCRIPTION rows for spw %d: spw %d will be assigned',
                    base_spw,
                    new_spw,
                )
                self.__remove_preexisting_rows(tb, new_spw)

                taql = f'SPECTRAL_WINDOW_ID == {base_spw}'
                selected = tb.query(taql)
                try:
                    base_dds = [int(r) for r in selected.rownumbers()]
                    if len(base_dds) == 0:
                        logger.warning('No DATA_DESCRIPTION row found for SPW %d', base_spw)
                        continue

                    startrow = tb.nrows()
                    selected.copyrows(tb.name())
                    nrow = tb.nrows() - startrow
                    spwcol = np.zeros(nrow, dtype=int) + new_spw
                    tb.putcol('SPECTRAL_WINDOW_ID', spwcol, startrow, nrow)

                    new_dds = list(range(startrow, startrow + nrow))
                    for base_dd, new_dd in zip(base_dds, new_dds):
                        self.dd_pairs.append((base_spw, new_spw, base_dd, new_dd))

                    self.spw_dd_map[new_spw] = new_dds
                finally:
                    selected.close()

        self.extra_dd = [new_dd for _, _, _, new_dd in self.dd_pairs]

    def _expand_syscal(self):
        self.__expand_subtable('SYSCAL')

    def _expand_feed(self):
        self.__expand_subtable('FEED')

    def _expand_source(self):
        self.__expand_subtable('SOURCE')

    def _expand_main(self):
        table_name = self.vis
        with sdutil.table_manager(table_name, nomodify=False) as tb:
            for base_spw, new_spw, base_dd, new_dd in self.dd_pairs:
                logger.info(
                    'duplicating MAIN rows for dd %d (spw %d): dd %d (spw %d) will be assigned',
                    base_dd,
                    base_spw,
                    new_dd,
                    new_spw,
                )

                startrow = tb.nrows()
                taql = f'DATA_DESC_ID == {base_dd}'
                self.__copy_selected_rows(tb, taql)
                nrow = tb.nrows() - startrow
                if nrow > 0:
                    ddcol = np.zeros(nrow, dtype=int) + new_dd
                    tb.putcol('DATA_DESC_ID', ddcol, startrow, nrow)

    def _duplicate(self):
        # duplicate science & atm spws
        # process SPECTRAL_WINDOW
        self._expand_spectral_window()

        # process DATA_DESCRIPTION
        self._expand_data_description()

        # process SYSCAL
        self._expand_syscal()

        # process FEED
        self._expand_feed()

        # process SOURCE
        self._expand_source()

        # process MAIN
        self._expand_main()
