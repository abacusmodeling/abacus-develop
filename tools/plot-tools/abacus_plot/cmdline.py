'''
Date: 2021-08-21 11:46:01
LastEditors: jiyuyang
LastEditTime: 2021-08-21 20:42:44
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from os import PathLike
from typing import Dict, List, Sequence, Union

from abacus_plot.utils import read_json


class Show:
    """Show auto-test information"""

    @classmethod
    def show_bandinfo(cls, datafile: Union[PathLike, Sequence[PathLike]], kptfile: PathLike, efermi: Union[float, Sequence[float]] = None, energy_range: Sequence[float] = [], blabel: Union[str, Sequence[str]] = None, color: Union[str, Sequence[str]] = None, outfile: PathLike = "band.png"):
        """Show band structure information

        :params datafile: path of band date file 
        :params kptfile: path of k-points file
        :params efermi: Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params blabel: band labels, its length equals to `filename`.
        :params color: band colors, its length equals to `filename`.
        :params outfile: band picture file name. Default: 'band.png'
        """

        from abacus_plot.band import BandPlot

        if isinstance(datafile, (str, PathLike)):
            BandPlot.singleplot(datafile, kptfile, efermi,
                                energy_range, blabel, color, outfile)
        elif isinstance(datafile, (list, tuple)):
            BandPlot.multiplot(datafile, kptfile, efermi,
                               energy_range, blabel, color, outfile)

    @classmethod
    def show_dosinfo(cls, tdosfile: PathLike = '', pdosfile: PathLike = '', efermi: float = 0, energy_range: Sequence[float] = [], dos_range: Sequence[float] = [], species: Union[Sequence[str], Dict[str, List[int]]] = [], tdosfig: PathLike = 'tdos.png', pdosfig: PathLike = 'pdos.png', prec: float = 0.01):
        """Plot total dos or partial dos, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`

        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params dos_range: range of dos to plot, its length equals to two
        :params species: list of atomic species or dict of atomic species and its angular momentum list
        :params prec: dos below this value thought to be zero. Default: 0.01
        """

        from abacus_plot.dos import DosPlot

        DosPlot().plot(tdosfile, pdosfile, efermi, energy_range,
                       dos_range, species, tdosfig, pdosfig, prec)

    @classmethod
    def show_cmdline(cls, args):
        if args.band:
            text = read_json(args.band)
            filename = text["filename"]
            kptfile = text["kptfile"]
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            blabel = text.pop("blabel", None)
            color = text.pop("color", None)
            outfile = text.pop("outfile", "band.png")
            cls.show_bandinfo(filename, kptfile, efermi,
                              energy_range, blabel, color, outfile)

        if args.dos:
            text = read_json(args.dos)
            tdosfile = text.pop("tdosfile", '')
            pdosfile = text.pop("pdosfile", '')
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            dos_range = text.pop("dos_range", [])
            species = text.pop("species", [])
            tdosfig = text.pop("tdosfig", "tdos.png")
            pdosfig = text.pop("pdosfig", "pdos.png")
            prec = text.pop("prec", 0.01)
            cls.show_dosinfo(tdosfile, pdosfile, efermi,
                             energy_range, dos_range, species, tdosfig, pdosfig, prec)
