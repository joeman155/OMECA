#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# RP-1 Thermal conductivity
#
# This loads data from "Table 10" on page 53 of https://www.govinfo.gov/content/pkg/GOVPUB-C13-18e6f8b04c5baea8df583f63c83ce583/pdf/GOVPUB-C13-18e6f8b04c5baea8df583f63c83ce583.pdf
# into a two dimensional array. Then we use interp2d to allow to find thermal conductivites for various pressures/temperatures.
#
# the data in file is space separated values:-
# "Point ID   T0(K)   Te(K)   Pe(MPa)   ρe(g·cm-3)   λe(W·m-1K-1)   λc(W·m-1K-1)"
#  numeric
#  id
#
# We are interested in Te, Pe and λe
#
# 

# Initialize and import packages
import numpy as np
from scipy.interpolate import interp2d
from scipy import interpolate


def interCond(k, temps, pressure, temperature):
    max_temp = np.max(temps)
    min_temp = np.min(temps)
    if temperature > max_temp or temperature < min_temp:
        print("Temperature is outside bounds of data")
        return
    lower_temp = 0
    upper_temp = 0
    lower_idx = 0
    upper_idx = 0
    idx = 0
    for temp in temps:
        # print("Comparing ", temperature, " with ", temp)
        if temperature < float(temp) and lower_temp > 0:
            upper_temp = previous_temp
            upper_idx = idx - 1
            break
        if temperature < float(temp) and lower_temp == 0:
            lower_temp = float(previous_temp)
            lower_idx = idx - 1
        previous_temp = float(temp)
        idx = idx + 1

    if lower_temp > 0 and upper_temp > 0:
        # Lower Temp side
        pressures = list(k[lower_idx].keys())
        conductivities = list(k[lower_idx].values())

        f = interpolate.interp1d(pressures, conductivities)
        lower_kval = f(pressure)

        # Upper Temp side
        pressures = list(k[upper_idx].keys())
        conductivities = list(k[upper_idx].values())
        f = interpolate.interp1d(pressures, conductivities)
        upper_kval = f(pressure)

        # Interpolate between the two temperatures, to give us the final result
        f = interpolate.interp1d([lower_temp, upper_temp], [lower_kval, upper_kval])
        kval = f(temperature)

        return kval

    else:
        print("Unable to get temperature bounded")
        return


def loadConductivityData():
    file = "rp1heatk.txt"
    data = open(file)
    content = data.readlines()
    data.close()

    temps = np.array([], float)
    k = []

    lnum = 0
    idx = 0
    for line in content:
        # print(line)
        values = line.split(" ")
        if lnum > 0:
            if values[0][0] != index_prefix:
                idx = idx + 1
                temp = float(values[2])  # Adopt first temp as "index" (close enough)
                temps = np.append(temps, temp)
                k.append(data)
                data = {}
        else:
            temp = values[2]  # Adopt first temp as "index" (close enough)
            temps = np.append(temps, float(temp))
            data = {}
        # print("For Index: ", values[0], " Temp: ", values[2], ", Pressure: ", values[3], " we have conductivity of: ", values[5])
        # print(values[0][0], ", Index: ", idx, ", index temp: ", temp)

        # idx  = row  (TEMP)
        # lnum = col (PRESSURE)
        data[float(values[3])] = float(values[5])
        # print(data)
        lnum = lnum + 1
        index_prefix = values[0][0]

    return k, temps


k, temps = loadConductivityData()
print(temps)

pressure = 4.0
temperature = 460
print("Finding thermal conductivity for Pressure: ", pressure, " MPa, Temperature: ", temperature, " Kelvin")
kval = interCond(k, temps, pressure, temperature)
print("Final thermal conductivity: ", kval)
