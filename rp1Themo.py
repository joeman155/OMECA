#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# RP-1 Thermodynamic properties
#


# Initialize and import packages
import numpy as np
from scipy.interpolate import interp2d
from scipy import interpolate



# Define rp-1 constants
rp1_M = 123     # Grams per mole
# Define
# other
# variables
#


#
# Load specific heat capacity data data from file into memory structures
#
# Example usage
#
#
# Range of data
# Temperatures from ??? K to ??? K
# Pressures    from ??? MPa to ??? MPa
#



#
# Load density data from file into memory structures
#
# Example usage
#
#
# Range of data
# Temperatures from 270 to 470K
# Pressures    from 0.083MPa to 40 MPa
#
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
    return f


#
# Load density data from file into memory structures
#
# Example usage
#
#
# Range of data
# Temperatures from 270 to 470K
# Pressures    from 0.083MPa to 40 MPa
`
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

    f = interp2d(pressures, temps, viscosities, kind='linear', fill_value='-1')
    return f



#
# Load thermal conductivity data from file into memory structures
#
# Example usage
#
#     k, temps = loadConductivityData()
#
#
#  See interCond function below for example of how to interogate this data to get thermal conducitivies at other temperatures/pressures
#
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

        # idx  = row  (TEMP)
        # lnum = col (PRESSURE)
        data[float(values[3])] = float(values[5])
        # print(data)
        lnum = lnum + 1
        index_prefix = values[0][0]

    return k, temps


#
# Routine used to interpolate thermal conducvitiy data structures to get value
# Couldn't use normal scipy routines unfortunately.
#
# Example usage
#
#     k, temps = loadConductivityData()
#     pressure = 4.0
#     temperature = 460
#     print("Finding thermal conductivity for Pressure: ", pressure, " MPa, Temperature: ", temperature, " Kelvin")
#     kval = interCond(k, temps, pressure, temperature)
#     print("Final thermal conductivity: ", kval)
#
# Range of data
# Temperatures from 300 to 645K
# Pressures    from 0.17MPa to 63 MPa
#
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