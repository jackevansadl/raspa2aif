# -*- coding: utf-8 -*-
"""Plot AIF from command line"""
import sys
import os
from gemmi import cif  # pylint: disable-msg=no-name-in-module
import matplotlib.pyplot as plt
import numpy as np

FILENAME = 'test.aif'

aif = cif.read(FILENAME)
block = aif.sole_block()
ads_press = np.array(block.find_loop('_adsorp_pressure'), dtype=float)
ads_amount_absolute = np.array(block.find_loop('_adsorp_amount_absolute'), dtype=float)
ads_amount_excess = np.array(block.find_loop('_adsorp_amount_excess'), dtype=float)

#material_id = block.find_pair('_sample_material_id')[-1]

plt.plot(ads_press, ads_amount_absolute, 'o-', color='C0',label='absolute')
plt.plot(ads_press, ads_amount_excess, 'o-', color='C1', label='excess')
plt.legend()
plt.ylabel('quantity adsorbed / ' + block.find_pair('_units_loading')[-1])
plt.xlabel('pressure / ' + block.find_pair('_units_pressure')[-1])

plt.show()
