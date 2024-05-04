#!/usr/bin/env python3

from xbout import open_boutdataset

ds = open_boutdataset("data/BOUT.dmp.*.nc", info=False).squeeze()

#ds.bout.animate_list([["n", "V_i", "pperp_i", "par_i", "pperp_e", "ppar_e"]], vmin=-0.2, vmax=0.4, show=True)
ds.bout.animate_list([["n", "V_i", "pperp_i", "ppar_i", "pperp_e", "ppar_e"]], show=True)
