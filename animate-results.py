#!/usr/bin/env python3

from xbout import open_boutdataset

ds = open_boutdataset("data/BOUT.dmp.*.nc", info=False).squeeze()

#ds.bout.animate_list([["n", "V_i", "pperp_i", "ppar_i", "pperp_e", "ppar_e", "Epar"]], vmin=-0.1,
#                     vmax=1.4, show=True)
ds.bout.animate_list([["n", "pperp_i", "ppar_i", "pperp_e", "ppar_e"]], show=True)

ds.bout.animate_list([["V_i", "Epar"]], show=True)
