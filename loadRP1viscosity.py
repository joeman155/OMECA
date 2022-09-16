#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# RP-1 Viscosity, pressure, temperature characteristics
#


# Initialize and import packages
import numpy as np
from scipy.interpolate import interp2d
from scipy import interpolate



def loadViscosityData():
    file = "rp1ViscosityData.txt"
    data = open(file)
    content = data.readlines()
    data.close()

    pressures = np.array([0.101, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], float)
    temps = np.array([], float)
    viscosities = np.empty((11, 7), float)

    lnum = 0
    for line in content:
        # print(line)
        idx = 0
        col = 0
        for record in line.split(" "):
            if idx == 0:
                temp = record
                temps = np.append(temps, float(temp))
                # print("Temp: ", temp)
            else:
                # print("Appending: ", float(record), " to ", lnum, ", ", col)
                viscosities[lnum][col] = float(record)
                col = col + 1
            idx = idx + 1

        lnum = lnum + 1

    # print("Pressures")
    # print(pressures)

    # print("Temps")
    # print(temps)

    f = interp2d(pressures, temps, viscosities, kind='linear', fill_value='-1')

    return f

f = loadViscosityData()


t1 = 360
p1 = 5
print("Viscosity of RP-1 at ", t1, "degrees Kelvin, and pressure ", p1, " MPa, has a density of ", f(p1, t1), " mPa s")


