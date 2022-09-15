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
    file = "rp1ViscosityATM.txt"
    data = open(file)
    contentATM = data.readlines()
    data.close()

    file = "rp1ViscosityHP.txt"
    data = open(file)
    contentHP = data.readlines()
    data.close()

    pressures = np.array([0.1, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], float)
    temps = np.array([], float)
    v = []

    lnum = 0
    for line in contentHP:
        # print(line)
        idx = 0
        for record in line.split(" "):
            if idx == 0:
                temp = record
                temps = np.append(temps, float(temp))
                # print("Temp: ", temp)
                viscosities = np.array([])
            else:
                viscosities = np.append(viscosities, float(record))
                # print("Appending: ", float(record))
                # print(viscosities)
            idx = idx + 1

        v.append(viscosities)

    print(pressures)
    print(temps)
    print(v)

loadViscosityData()

