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


def interCond(pressure, temperature):
	print("Finding thermal conductivity for Pressure: ", pressure, " MPa, Temperature: ", temperature, " Kelvin")







file = "rp1heatk.txt"
data = open(file)
content = data.readlines()
data.close()


temps = np.array([])
# k = np.empty((14,11), float)   # rows, cols
k = []




lnum = 0
idx = 0
for line in content:
    # print(line)
    values = line.split(" ")
    if lnum > 0:
       if values[0][0] != index_prefix:
           idx = idx + 1
           temp = values[2]  # Adopt first temp as "index" (close enough)
           temps = np.append(temps, temp)
           k.append(data)
           data = {}
    else:
       temp = values[2]      # Adopt first temp as "index" (close enough)
       temps = np.append(temps, temp)
       data = {}
    # print("For Index: ", values[0], " Temp: ", values[2], ", Pressure: ", values[3], " we have conductivity of: ", values[5])
    # print(values[0][0], ", Index: ", idx, ", index temp: ", temp)

    # idx  = row  (TEMP)
    # lnum = col (PRESSURE)
    data[values[3]] = values[5]
    # print(data)
    lnum = lnum + 1 
    index_prefix = values[0][0]


print(k)

for t in temps:
    print(t)


interCond(5, 330)


