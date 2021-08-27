'''
Date: 2021-05-08 11:47:09
LastEditors: jiyuyang
LastEditTime: 2021-08-26 12:07:22
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
from os import PathLike
from typing import Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import axes

from abacus_plot.utils import energy_minus_efermi, list_elem2str, read_kpt


class BandPlot:
    """Plot band structure"""

    @classmethod
    def set_vcband(cls, energy: Sequence) -> Tuple[namedtuple, namedtuple]:
        """Separate valence and conduct band

        :params energy: band energy after subtracting the Fermi level
        """

        e_T = energy.T
        num_gt_Ef = (e_T > 0).sum(axis=1)

        Band = namedtuple('Band', ['band_index', 'band', 'value', 'k_index'])

        # valance band
        band_vbm_index = np.where(num_gt_Ef == 0)[0]
        band_vbm = e_T[band_vbm_index]
        evbm = np.max(band_vbm)
        k_vbm_index = np.where(band_vbm == evbm)[1]
        vb = Band(band_vbm_index, band_vbm-evbm, evbm-evbm, k_vbm_index)

        # conduct band
        band_cbm_index = np.where(num_gt_Ef != 0)[0]
        band_cbm = e_T[band_cbm_index]
        ecbm = np.min(band_cbm)
        k_cbm_index = np.where(band_cbm == ecbm)[1]
        cb = Band(band_cbm_index, band_cbm-evbm, ecbm-evbm, k_cbm_index)

        return vb, cb

    @classmethod
    def read(cls, filename: PathLike) -> Tuple[np.ndarray, np.ndarray]:
        """Read band data file and return k-points and energy

        :params filename: string of band data file
        """

        data = np.loadtxt(filename)
        X, y = np.split(data, (1, ), axis=1)
        x = X.flatten()
        return x, y

    @classmethod
    def _set_figure(cls, ax: axes.Axes, index: dict, range: Sequence):
        """set figure and axes for plotting

        :params ax: matplotlib.axes.Axes object
        :params index: dict of label of points of x-axis and its index in data file. Range of x-axis based on index.value()
        :params range: range of y-axis
        """

        keys = []
        values = []
        for t in index:
            if isinstance(t, tuple):
                keys.append(t[0])
                values.append(t[1])
            elif isinstance(t, (int, float)):
                keys.append('')
                values.append(t)

        # x-axis
        ax.set_xticks(values)
        ax.set_xticklabels(keys)
        ax.set_xlim(values[0], values[-1])
        ax.set_xlabel("Wave Vector")

        # y-axis
        if range:
            ax.set_ylim(range[0], range[1])
        ax.set_ylabel(r"$E-E_{fermi}(eV)$")

        # others
        ax.grid(axis='x', lw=1.2)
        ax.axhline(0, linestyle="--", c='b', lw=1.0)
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    @classmethod
    def plot(cls, x: Sequence, y: Sequence, index: Sequence, efermi: float = 0, energy_range: Sequence[float] = [], label: str = None, color: str = None, outfile: PathLike = 'band.png'):
        """Plot band structure

        :params x, y: x-axis and y-axis coordinates
        :params index: special k-points label and its index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params label: band label. Default: ''
        :params color: band color. Default: 'black'
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()

        if not color:
            color = 'black'

        kpoints, energy = x, y
        energy = energy_minus_efermi(energy, efermi)

        ax.plot(kpoints, energy, lw=0.8, color=color, label=label)
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)

    @classmethod
    def singleplot(cls, datafile: PathLike, kptfile: str = [], efermi: float = 0, energy_range: Sequence[float] = [], shift: bool = False, label: str = None, color: str = None, outfile: PathLike = 'band.png'):
        """Plot band structure using data file

        :params datafile: string of band date file
        :params kptfile: k-point file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        :params label: band label. Default: ''
        :params color: band color. Default: 'black'
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()
        kpt = read_kpt(kptfile)

        if not color:
            color = 'black'

        kpoints, energy = cls.read(datafile)
        if shift:
            vb, cb = cls.set_vcband(energy_minus_efermi(energy, efermi))
            ax.plot(kpoints, np.vstack((vb.band, cb.band)).T,
                    lw=0.8, color=color, label=label)
            cls.info(kpt.full_kpath, vb, cb)
        else:
            ax.plot(kpoints, energy_minus_efermi(energy, efermi),
                    lw=0.8, color=color, label=label)
        index = kpt.label_special_k
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)

    @classmethod
    def multiplot(cls, datafile: Sequence[PathLike], kptfile: str = '', efermi: Sequence[float] = [], energy_range: Sequence[float] = [], shift: bool = True, label: Sequence[str] = None, color: Sequence[str] = None, outfile: PathLike = 'band.png'):
        """Plot more than two band structures using data file

        :params datafile: list of path of band date file 
        :params kptfile: k-point file
        :params efermi: list of Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        :params label: list of band labels, its length equals to `filename`
        :params color: list of band colors, its length equals to `filename`
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()
        kpt = read_kpt(kptfile)

        if not efermi:
            efermi = [0.0 for i in range(len(datafile))]
        if not label:
            label = ['' for i in range(len(datafile))]
        if not color:
            color = ['black' for i in range(len(datafile))]

        emin = -np.inf
        emax = np.inf
        for i, file in enumerate(datafile):
            kpoints, energy = cls.read(file)
            if shift:
                vb, cb = cls.set_vcband(energy_minus_efermi(energy, efermi[i]))
                energy_min = np.min(vb.band)
                energy_max = np.max(cb.band)
                if energy_min > emin:
                    emin = energy_min
                if energy_max < emax:
                    emax = energy_max

                ax.plot(kpoints, np.vstack((vb.band, cb.band)).T,
                        lw=0.8, color=color[i], label=label[i])
                cls.info(kpt.full_kpath, vb, cb)
            else:
                ax.plot(kpoints, energy_minus_efermi(energy, efermi[i]),
                        lw=0.8, color=color[i], label=label[i])

        index = kpt.label_special_k
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap

    @classmethod
    def info(cls, kpath: Sequence, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        def band_type(vbm_x, cbm_x):
            longone, shortone = (vbm_x, cbm_x) if len(
                vbm_x) >= len(cbm_x) else (cbm_x, vbm_x)
            for i in shortone:
                if i in longone:
                    btype = "Direct"
                else:
                    btype = "Indirect"
            return btype

        gap = cls.bandgap(vb, cb)
        print(
            "--------------------------Band Structure--------------------------", flush=True)
        print(
            f"{'Band character:'.ljust(30)}{band_type(vb.k_index, cb.k_index)}", flush=True)
        print(f"{'Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)
        print(f"{'Band index:'.ljust(30)}{'HOMO'.ljust(10)}{'LUMO'}", flush=True)
        print(
            f"{''.ljust(30)}{str(vb.band_index[-1]).ljust(10)}{str(cb.band_index[0])}", flush=True)
        print(f"{'Eigenvalue of VBM(eV):'.ljust(30)}{vb.value: .4f}", flush=True)
        print(f"{'Eigenvalue of CBM(eV):'.ljust(30)}{cb.value: .4f}", flush=True)
        vbm_k = np.unique(kpath[vb.k_index], axis=0)
        cbm_k = np.unique(kpath[cb.k_index], axis=0)
        print(
            f"{'Location of VBM'.ljust(30)}{' '.join(list_elem2str(vbm_k[0]))}", flush=True)
        for i, j in enumerate(vbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)
        print(
            f"{'Location of CBM'.ljust(30)}{' '.join(list_elem2str(cbm_k[0]))}", flush=True)
        for i, j in enumerate(cbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)
