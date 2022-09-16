#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# RP-1 Density, pressure, temperature characteristics
#
# This loads data from "Table 2" on document https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=832303
# into a two dimensional array. Then we use interp2d to allow to find densities for various pressures/temperatures.
#
# In some cases, the pressures are of by 0.01, but we ignore this and use pressures as defined in pressures lis
#
# 

# Initialize and import packages
import numpy as np
from scipy.interpolate import interp2d


def loadDensityData():
    PRHO_file = "rp1DensityData.txt"
    data = open(PRHO_file)
    content = data.readlines()
    data.close()

    temps = [270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470]
    pressures = [40, 35, 30., 25.0, 20, 15, 10, 6, 5, 3, 2, 1, 0.5, 0.083]
    density_array = np.empty((14, 11), float)

    lnum = 0
    for line in content:
        pair = 0
        idx = 0
        for record in line.split(" "):
            pair = pair + 1
            if pair % 2 == 0:
                density = record
                density_array[lnum][idx] = density
                idx = idx + 1
            else:
                pressure = record
        lnum = lnum + 1

    f = interp2d(temps, pressures, density_array, kind='linear', fill_value='-1')
    # print(density_array)
    return f


# Load RP-1 Density data
f = loadDensityData()


t1 = 360
p1 = 5
print("Density of RP-1 at ", t1, "degrees Kelvin, and pressure ", p1, " MPa, has a density of ", f(t1, p1), " kg/m^3")
