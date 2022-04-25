'''
Date: 2021-12-29 10:27:01
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:06:14
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
import numpy as np
from os import PathLike
from typing import Sequence, Tuple, Union, Any, List, Dict
from matplotlib.figure import Figure
from matplotlib import axes

from abacus_plot.utils import energy_minus_efermi, list_elem2str, read_kpt, remove_empty, parse_projected_data, handle_data


class Band:
    """Parse Bands data"""

    def __init__(self, bandfile: Union[PathLike, Sequence[PathLike]] = None, kptfile: str = '') -> None:
        self.bandfile = bandfile
        self.read()
        self.kptfile = kptfile
        self.kpt = None
        if self.kptfile:
            self.kpt = read_kpt(kptfile)
        self.nspin = None
        self.norbitals = None
        self.orbitals = []

    def read(self):
        """Read band data file and return k-points and energy

        :params filename: string of band data file
        """

        data = np.loadtxt(self.bandfile)
        X, self.k_index = np.split(data, (1, ), axis=1)
        self.nkpoints = len(self.k_index)
        self.energy = X.flatten()
        self.nbands = self.energy.shape[-1]
        self.eunit = 'eV'

    @classmethod
    def direct_bandgap(cls, vb: namedtuple, cb: namedtuple, klength: int):
        """Calculate direct band gap"""

        gap_list = []
        i_index = []
        for i in range(klength):
            gap_list.append(np.min(cb.band[:, i])-np.max(vb.band[:, i]))
            i_index.append(i)
        dgap = np.min(gap_list)

        return dgap, i_index[np.argmin(gap_list)]

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap

    @classmethod
    def band_type(cls, vb: namedtuple, cb: namedtuple):
        vbm_x, cbm_x = vb.k_index, cb.k_index
        longone, shortone = (vbm_x, cbm_x) if len(
            vbm_x) >= len(cbm_x) else (cbm_x, vbm_x)
        for i in shortone:
            if i in longone:
                btype = "Direct"
            else:
                btype = "Indirect"

        return btype

    @classmethod
    def info(cls, kpath: Sequence, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        gap = cls.bandgap(vb, cb)
        dgap, d_i = cls.direct_bandgap(vb, cb, len(kpath))
        btype = cls.band_type(vb, cb)
        print(
            "--------------------------Band Structure--------------------------", flush=True)
        print(
            f"{'Band character:'.ljust(30)}{btype}", flush=True)
        if btype == "Indirect":
            print(f"{'Direct Band gap(eV):'.ljust(30)}{dgap: .4f}", flush=True)
            print(f"{'Indirect Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)
        elif btype == "Direct":
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

    def plot(self, x: Sequence, y: Sequence, index: Sequence, efermi: float = 0, energy_range: Sequence[float] = []):
        """Plot band structure

        :params x, y: x-axis and y-axis coordinates
        :params index: special k-points label and its index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        """

        if not self._color:
            self._color = 'black'

        kpoints, energy = x, y
        energy = energy_minus_efermi(energy, efermi)

        self.ax.plot(kpoints, energy, lw=self._lw, color=self._color,
                     label=self._label, linestyle=self._linestyle)
        self._set_figure(index, energy_range)

    def singleplot(self, efermi: float = 0, energy_range: Sequence[float] = [], shift: bool = False):
        """Plot band structure using data file

        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        """

        if not self._color:
            self._color = 'black'

        energy = energy_minus_efermi(self.energy, efermi)
        if shift:
            vb, cb = self.set_vcband(energy)
            self.ax.plot(self.k_index, np.vstack((vb.band, cb.band)).T,
                         lw=self._lw, color=self._color, label=self._label, linestyle=self._linestyle)
            self.info(self.kpt.full_kpath, vb, cb)
        else:
            self.ax.plot(self.k_index, energy,
                         lw=self._lw, color=self._color, label=self._label, linestyle=self._linestyle)
        index = self.kpt.label_special_k
        self._set_figure(index, energy_range)

    def multiplot(self, efermi: Sequence[float] = [], energy_range: Sequence[float] = [], shift: bool = True):
        """Plot more than two band structures using data file

        :params efermi: list of Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        """

        if not efermi:
            efermi = [0.0 for i in range(len(self.bandfile))]
        if not self._label:
            self._label = ['' for i in range(len(self.bandfile))]
        if not self._color:
            self._color = ['black' for i in range(len(self.bandfile))]
        if not self._linestyle:
            self._linestyle = ['solid' for i in range(len(self.bandfile))]

        emin = -np.inf
        emax = np.inf
        for i, file in enumerate(self.bandfile):
            if shift:
                vb, cb = self.set_vcband(
                    energy_minus_efermi(self.energy, efermi[i]))
                energy_min = np.min(vb.band)
                energy_max = np.max(cb.band)
                if energy_min > emin:
                    emin = energy_min
                if energy_max < emax:
                    emax = energy_max

                self.ax.plot(self.k_index, np.vstack((vb.band, cb.band)).T,
                             lw=self._lw, color=self._color[i], label=self._label[i], linestyle=self._linestyle[i])
                self.info(self.kpt.full_kpath, vb, cb)
            else:
                self.ax.plot(self.k_index, energy_minus_efermi(self.energy, efermi[i]),
                             lw=self._lw, color=self._color[i], label=self._label[i], linestyle=self._linestyle[i])

        index = self.kpt.label_special_k
        self._set_figure(index, energy_range)


class BandPlot:
    """Plot band structure"""

    def __init__(self, fig: Figure, ax: axes.Axes, **kwargs) -> None:
        self.fig = fig
        self.ax = ax
        self._lw = kwargs.pop('lw', 2)
        self._bwidth = kwargs.pop('bwdith', 3)
        self._label = kwargs.pop('label', None)
        self._color = kwargs.pop('color', None)
        self._linestyle = kwargs.pop('linestyle', 'solid')
        self.plot_params = kwargs

    def _set_figure(self, index: dict, range: Sequence):
        """set figure and axes for plotting

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
        self.ax.set_xticks(values)
        self.ax.set_xticklabels(keys)
        self.ax.set_xlim(values[0], values[-1])
        if "xlabel_params" in self.plot_params.keys():
            self.ax.set_xlabel(
                "Wave Vector", **self.plot_params["xlabel_params"])
        else:
            self.ax.set_xlabel("Wave Vector", size=25)

        # y-axis
        if range:
            self.ax.set_ylim(range[0], range[1])
        if "ylabel_params" in self.plot_params.keys():
            self.ax.set_ylabel(
                "Energy(eV)", **self.plot_params["ylabel_params"])
        else:
            self.ax.set_ylabel("Energy(eV)", size=25)

        # notes
        if "notes" in self.plot_params.keys():
            from matplotlib.offsetbox import AnchoredText
            if "s" in self.plot_params["notes"].keys() and len(self.plot_params["notes"].keys()) == 1:
                self.ax.add_artist(AnchoredText(self.plot_params["notes"]["s"], loc='upper left', prop=dict(size=25),
                                                borderpad=0.2, frameon=False))
            else:
                self.ax.add_artist(AnchoredText(**self.plot_params["notes"]))

        # ticks
        if "tick_params" in self.plot_params.keys():
            self.ax.tick_params(**self.plot_params["tick_params"])
        else:
            self.ax.tick_params(labelsize=25)

        # frame
        bwidth = self._bwidth
        self.ax.spines['top'].set_linewidth(bwidth)
        self.ax.spines['right'].set_linewidth(bwidth)
        self.ax.spines['left'].set_linewidth(bwidth)
        self.ax.spines['bottom'].set_linewidth(bwidth)

        # guides
        if "grid_params" in self.plot_params.keys():
            self.ax.grid(axis='x', **self.plot_params["grid_params"])
        else:
            self.ax.grid(axis='x', lw=1.2)
        if "hline_params" in self.plot_params.keys():
            self.ax.axhline(0, **self.plot_params["hline_params"])
        else:
            self.ax.axhline(0, linestyle="--", c='b', lw=1.0)

        if self._label:
            handles, labels = self.ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if "legend_prop" in self.plot_params.keys():
                self.ax.legend(by_label.values(), by_label.keys(),
                               prop=self.plot_params["legend_prop"])
            else:
                self.ax.legend(by_label.values(),
                               by_label.keys(), prop={'size': 15})


class PBand(Band):
    def __init__(self, bandfile: PathLike = None) -> None:
        self.bandfile = bandfile
        self.read()

    def _check_energy(self):
        assert self.energy.shape[0] == self.nkpoints, "The dimension of band structure dismatches with the number of k-points."
        assert self.energy.shape[1] == self.nbands, "The dimension of band structure dismatches with the number of bands."

    def _check_weights(self, weights:np.ndarray, prec=1e-5):
        assert weights.shape[0] == self.norbitals, "The dimension of weights dismatches with the number of orbitals."
        assert weights.shape[1] == self.nkpoints, "The dimension of weights dismatches with the number of k-points."
        assert weights.shape[2] == self.nbands, "The dimension of weights dismatches with the number of bands."
        one_mat = np.ones((self.nkpoints, self.nbands))
        assert (np.abs(weights.sum(axis=0)-one_mat) < prec).all(), f"np.abs(weights.sum(axis=0)-np.ones(({self.nkpoints}, {self.nbands}))) < {prec}"

    @property
    def weights(self):
        data = np.empty((self.norbitals, self.nkpoints, self.nbands))
        for i, orb in enumerate(self.orbitals):
            data[i] = orb['data']
        self._check_weights(data)
        return data

    def read(self):
        """Read projected band data file and return k-points, energy and Mulliken weights

        :params filename: string of projected band data file
        """

        from lxml import etree
        pbanddata = etree.parse(self.bandfile)
        root = pbanddata.getroot()
        self.nspin = int(root.xpath('//nspin')[0].text.replace(' ', ''))
        self.norbitals = int(root.xpath('//norbitals')
                             [0].text.replace(' ', ''))
        self.eunit = root.xpath('//band_structure/@units')[0].replace(' ', '')
        self.nbands = int(root.xpath('//band_structure/@nbands')
                          [0].replace(' ', ''))
        self.nkpoints = int(root.xpath('//band_structure/@nkpoints')
                            [0].replace(' ', ''))
        self.k_index = np.arange(self.nkpoints)
        self.energy = root.xpath('//band_structure')[0].text.split('\n')
        self.energy = handle_data(self.energy)
        remove_empty(self.energy)
        self.energy = np.asarray(self.energy, dtype=float)
        self._check_energy()

        self.orbitals = []
        for i in range(self.norbitals):
            orb = OrderedDict()
            o_index_str = root.xpath(
                '//orbital/@index')[i]
            orb['index'] = int(o_index_str.replace(' ', ''))
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
            self.orbitals.append(orb)

    def _write(self, species: Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], keyname='', outdir: PathLike = './'):
        """Write parsed partial dos data to files

        Args:
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            keyname (str): the keyword that extracts the PDOS. Allowed values: 'index', 'atom_index', 'species'
        """

        return 


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pathlib import Path
    parent = Path(r"C:\Users\YY.Ji\Desktop")
    name = "PBANDS_1"
    path = parent/name
    # notes = {'s': '(b)'}
    # datafile = [path/"soc.dat", path/"non-soc.dat"]
    # kptfile = path/"KPT"
    #fig, ax = plt.subplots(figsize=(12, 12))
    # label = ["with SOC", "without SOC"]
    # color = ["r", "g"]
    # linestyle = ["solid", "dashed"]
    # band = BandPlot(fig, ax, notes=notes, label=label,
    #                 color=color, linestyle=linestyle)
    # energy_range = [-5, 6]
    # efermi = [4.417301755850272, 4.920435541999894]
    # shift = True
    # band.multiplot(datafile, kptfile, efermi, energy_range, shift)
    # fig.savefig("band.png")
    pband = PBand(str(path))
