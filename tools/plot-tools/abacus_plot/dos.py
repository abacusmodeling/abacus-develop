'''
Date: 2021-08-18 11:05:00
LastEditors: jiyuyang
LastEditTime: 2021-08-21 21:31:47
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
from os import PathLike
from typing import Dict, List, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import axes

from abacus_plot.utils import (energy_minus_efermi, get_angular_momentum_label,
                               get_angular_momentum_name, remove_empty)


class DosPlot:
    """Plot density of state(DOS)"""

    @classmethod
    def _val_to_zero(cls, string):
        if eval(string) <= 1e-20:
            string = "0"
        return string

    @classmethod
    def set_vcband(cls, energy: Sequence, dos: Sequence, prec=0.01) -> Tuple[namedtuple, namedtuple]:
        """Separate valence and conduct band

        :params energy: band energy after subtracting the Fermi level
        :params dos: density of state
        :params prec: dos below this value thought to be zero. Default: 0.01
        """

        Band = namedtuple('Band', ['band', 'value'])

        # valance band
        band_vbm_index = np.where(energy <= 0)[0]
        evbm = energy[np.where(dos[band_vbm_index] > prec)[0][-1]][0]
        band_vbm = energy[band_vbm_index]
        vb = Band(band_vbm-evbm, evbm-evbm)

        # conduct band
        band_cbm_index = np.where(energy > 0)[0]
        band_cbm = energy[band_cbm_index]
        ecbm = energy[np.where(dos[band_cbm_index] > prec)[
            0][0]+band_cbm_index[0]][0]
        cb = Band(band_cbm-evbm, ecbm-evbm)

        return vb, cb

    @classmethod
    def info(cls, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        gap = cls.bandgap(vb, cb)
        print("--------------------------Density of State--------------------------", flush=True)
        print(f"{'Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap

    @classmethod
    def read(cls, tdosfile: PathLike = '', pdosfile: PathLike = '') -> tuple:
        """Read DOS data file, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`

        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        """

        if tdosfile:
            dosdata = np.loadtxt(tdosfile)
            nsplit = dosdata.shape[1]
            return np.split(dosdata, nsplit, axis=1)

        elif pdosfile:
            def handle_data(data):
                data.remove('')

                def handle_elem(elem):
                    elist = elem.split(' ')
                    remove_empty(elist)  # `list` will be modified in function
                    return elist
                return list(map(handle_elem, data))

            from lxml import etree
            pdosdata = etree.parse(pdosfile)
            root = pdosdata.getroot()
            pdosdata = etree.parse(pdosfile)
            root = pdosdata.getroot()
            nspin = int(root.xpath('//nspin')[0].text.replace(' ', ''))
            norbitals = int(root.xpath('//norbitals')[0].text.replace(' ', ''))
            eunit = root.xpath('//energy_values/@units')[0].replace(' ', '')
            e_list = root.xpath(
                '//energy_values')[0].text.replace(' ', '').split('\n')
            remove_empty(e_list)
            orbitals = []
            for i in range(norbitals):
                orb = OrderedDict()
                orb['index'] = int(root.xpath(
                    '//orbital/@index')[i].replace(' ', ''))
                orb['atom_index'] = int(root.xpath(
                    '//orbital/@atom_index')[i].replace(' ', ''))
                orb['species'] = root.xpath(
                    '//orbital/@species')[i].replace(' ', '')
                orb['l'] = int(root.xpath('//orbital/@l')[i].replace(' ', ''))
                orb['m'] = int(root.xpath('//orbital/@m')[i].replace(' ', ''))
                orb['z'] = int(root.xpath('//orbital/@z')[i].replace(' ', ''))
                data = root.xpath('//data')[i].text.split('\n')
                data = handle_data(data)
                remove_empty(data)
                orb['data'] = np.asarray(data, dtype=float)
                orbitals.append(orb)

            return np.reshape(e_list, newshape=(-1, 1)).astype(float), orbitals

    @classmethod
    def _set_figure(cls, ax: axes.Axes, energy_range: Sequence, dos_range: Sequence):
        """set figure and axes for plotting

        :params ax: matplotlib.axes.Axes object
        :params dos_range: range of dos
        :params energy_range: range of energy
        """

        # y-axis
        if dos_range:
            ax.set_ylim(dos_range[0], dos_range[1])
        ax.set_ylabel("DOS")

        # x-axis
        if energy_range:
            ax.set_xlim(energy_range[0], energy_range[1])
        ax.set_xlabel(r"$E-E_{fermi}(eV)$")

        # others
        ax.axvline(0, linestyle="--", c='b', lw=1.0)
        ax.legend()

    @classmethod
    def _tplot(cls, res: tuple, efermi: float = 0, energy_range: Sequence[float] = [], dos_range: Sequence[float] = [], prec: float = 0.01):

        fig, ax = plt.subplots()

        nsplit = len(res)
        if nsplit == 2:
            energy, dos = res
            vb, cb = cls.set_vcband(
                energy_minus_efermi(energy, efermi), dos, prec)
            energy = np.concatenate((vb.band, cb.band))
            ax.plot(energy, dos, lw=0.8, c='gray', linestyle='-', label='TDOS')

        elif nsplit == 3:
            energy, dos_up, dos_dw = res
            vb, cb = cls.set_vcband(
                energy_minus_efermi(energy, efermi), dos_up, prec)
            energy = np.concatenate((vb.band, cb.band))
            dos_dw = -dos_dw
            ax.plot(energy, dos_up, lw=0.8, c='gray',
                    linestyle='-', label=r'$TDOS \uparrow$')
            ax.plot(energy, dos_up, lw=0.8, c='gray',
                    linestyle='--', label=r'$TDOS \downarrow$')

        return ax, energy_range, dos_range

    @classmethod
    def _all_sum(cls, orbitals: dict) -> Tuple[np.ndarray, int]:
        nsplit = orbitals[0]["data"].shape[1]
        res = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
        for orb in orbitals:
            res = res + orb['data']
        return res, nsplit

    @classmethod
    def plot(cls, tdosfile: PathLike = '', pdosfile: PathLike = '', efermi: float = 0, energy_range: Sequence[float] = [], dos_range: Sequence[float] = [], species: Union[Sequence[str], Dict[str, List[int]], Dict[str, Dict[str, List[int]]]] = [], tdosfig: PathLike = 'tdos.png', pdosfig:  PathLike = 'pdos.png', prec: float = 0.01):
        """Plot total dos or partial dos, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`

        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params dos_range: range of dos to plot, its length equals to two
        :params species: list of atomic species or dict of atomic species and its angular momentum list
        :params prec: dos below this value thought to be zero. Default: 0.01
        """

        res = cls.read(tdosfile, pdosfile)

        if tdosfile:
            ax, energy_range, dos_range = cls._tplot(
                res, efermi, energy_range, dos_range, prec)
            cls._set_figure(ax, energy_range, dos_range)
            plt.savefig(tdosfig)

        elif pdosfile and species:
            if isinstance(species, (list, tuple)):
                elements = species
            elif isinstance(species, dict):
                elements = list(species.keys())
                l = list(species.values())
            if not elements:
                raise TypeError(
                    "Only when `pdosfile` and `species` are both set, it will plot PDOS.")
            energy, orbitals = res

            # TDOS
            dos, nsplit = cls._all_sum(orbitals)
            vb, cb = cls.set_vcband(
                energy_minus_efermi(energy, efermi), dos, prec)
            cls.info(vb, cb)
            energy_f = np.concatenate((vb.band, cb.band))
            if nsplit == 1:
                ax, energy_range, dos_range = cls._tplot(
                    (energy, dos), efermi, energy_range, dos_range, prec)
            elif nsplit == 2:
                dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                ax, energy_range, dos_range = cls._tplot(
                    (energy, dos_up, dos_dw), efermi, energy_range, dos_range, prec)
            for elem in elements:
                dos = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
                for orb in orbitals:
                    if orb["species"] == elem:
                        dos += orb["data"]
                if nsplit == 1:
                    ax.plot(energy_f, dos, lw=0.8,
                            linestyle='-', label=f'{elem}')
                elif nsplit == 2:
                    dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                    ax.plot(energy_f, dos_up, lw=0.8, linestyle="-",
                            label=f"{elem}"+r"$\uparrow$")
                    dos_dw = -dos_dw
                    ax.plot(energy_f, dos_dw, lw=0.8, linestyle="--",
                            label=f"{elem}"+r"$\downarrow$")

            cls._set_figure(ax, energy_range, dos_range)
            plt.savefig(tdosfig)

            # PDOS
            if l:
                fig, ax = plt.subplots(
                    len(elements), 1, sharex=True, sharey=True)
                if len(elements) == 1:
                    ax = [ax]
                plt.subplots_adjust(hspace=0)
                for i, elem in enumerate(elements):
                    if isinstance(l[i], dict):
                        for ang, mag in l[i].items():
                            l_index = int(ang)
                            for m_index in mag:
                                dos = np.zeros_like(
                                    orbitals[0]["data"], dtype=np.float32)
                                for orb in orbitals:
                                    if orb["species"] == elem and orb["l"] == l_index and orb["m"] == m_index:
                                        dos += orb["data"]
                                if nsplit == 1:
                                    ax[i].plot(energy_f, dos, lw=0.8, linestyle='-',
                                               label=f'{elem}-{get_angular_momentum_name(l_index, m_index)}')
                                elif nsplit == 2:
                                    dos_up, dos_dw = np.split(
                                        dos, nsplit, axis=1)
                                    ax[i].plot(energy_f, dos_up, lw=0.8, linestyle="-",
                                               label=f"{elem}-{get_angular_momentum_name(l_index, m_index)}"+r"$\uparrow$")
                                    dos_dw = -dos_dw
                                    ax[i].plot(energy_f, dos_dw, lw=0.8, linestyle="--",
                                               label=f"{elem}-{get_angular_momentum_name(l_index, m_index)}"+r"$\downarrow$")
                    elif isinstance(l[i], list):
                        for l_index in l[i]:
                            dos = np.zeros_like(
                                orbitals[0]["data"], dtype=np.float32)
                            for orb in orbitals:
                                if orb["species"] == elem and orb["l"] == l_index:
                                    dos += orb["data"]
                            if nsplit == 1:
                                ax[i].plot(energy_f, dos, lw=0.8, linestyle='-',
                                           label=f'{elem}-{get_angular_momentum_label(l_index)}')
                            elif nsplit == 2:
                                dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                                ax[i].plot(energy_f, dos_up, lw=0.8, linestyle="-",
                                           label=f"{elem}-{get_angular_momentum_label(l_index)}"+r"$\uparrow$")
                                dos_dw = -dos_dw
                                ax[i].plot(energy_f, dos_dw, lw=0.8, linestyle="--",
                                           label=f"{elem}-{get_angular_momentum_label(l_index)}"+r"$\downarrow$")

                    cls._set_figure(ax[i], energy_range, dos_range)

                plt.savefig(pdosfig)
