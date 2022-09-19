#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# RP-1 Thermodynamic properties
#


# Initialize and import packages
import numpy as np
from scipy.interpolate import interp2d
from scipy import interpolate


def bartz(T0, Tw, p0, M, Dt, area, visc, cp, Pr, gamma, cstar):
    s = (0.5 * Tw / T0 * (1 + (gamma - 1) / 2 * M ** 2) + 0.5) ** (-0.68) * (1 + (gamma - 1) / 2 * M ** 2) ** (-0.12)
    h = 0.026 / Dt ** 0.2 * visc ** 0.2 * cp / Pr ** 0.6 * (p0 / cstar) ** 0.8 * (1.0 / area) ** 0.9 * s
    return h


def bartzFilm(T0, Tw, p0, M, Dt, area, muFilm, cpFilm, PrFilm, gamma, cstar):
    s = (0.5 * Tw / T0 * (1 + (gamma - 1) / 2 * M ** 2) + 0.5) ** (-0.68)
    h = 0.026 / Dt ** 0.2 * muFilm ** 0.2 * cpFilm / PrFilm ** 0.6 * (p0 / cstar) ** 0.8 * (1.0 / area) ** 0.9 * s
    return h


def adiabatic_wall(T_free, gamma, M, Pr):
    r = Pr ** (1.0 / 3.0)
    return T_free * (1 + r * (gamma - 1) / 2 * M ** 2)


def interpol(Tv, param, Tf):
    intpFunction = interp1d(Tv[::-1], param[::-1])
    return intpFunction(Tf)


#
# See https://www.thermal-engineering.org/what-is-dittus-boelter-equation-definition/  for some details
#
#
def dittusBoelter(Re, Pr):
    if Re < 0.6 or Re > 160:
        print("Warning: Reynolds number (", Re, ") is outside the valid range of 0.6 ... 160. Best to assess suitability of dittusBoelter correlation for this job.")
    if Pr < 10000:
        print("Warning: Prandtl number (", Pr, ") is less than 10000. Best to assess suitability of dittusBoelter correlation for this job.")
    return 0.023 * Re ** 0.8 * Pr ** 0.4


def Taylor(Re, Pr, T, Tw, Dh, x):
    return 0.023 * Re ** 0.8 * Pr ** 0.4 * (T / Tw) ** (0.57 - 1.59 * Dh / x)


def Ruan(Re, Pr, rho, rhow, Dh, x):
    return 0.0069 * Re ** 0.9 * Pr ** 0.66 * (rhow / rho) ** 0.43 * (1 + 2.4 * Dh / x)


def Dh_shell(Do, Di):
    return Do - Di


def Dh_rect(w, h):
    return 2 * w * h / (w + h)


def colebrook(Dh, roughness, Re, f):
    return 1 / (-2 * np.log10(roughness / (3.7 * Dh) + 2.51 / (Re * np.sqrt(f)))) ** 2


def frictionFactor(Dh, roughness, Re):
    guess = 1e-5
    for i in range(5):
        guess = colebrook(Dh, roughness, Re, guess)
    return guess


class rp1thermo:
    # Define rp-1 constants
    rp1_M = 175  # Grams per mole
    cp = 2001    # kJ/Kg/K

    def __init__(self):
        print("instantiating class.")
        self.loadViscosityData()
        self.loadDensityData()
        self.loadConductivityData()


    #
    # Return the cp
    #
    def getCp(self):
        return self.cp

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
    def loadDensityData(self):
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

        self.densityfn = interp2d(temps, pressures, density_array, kind='linear', fill_value='-1')


    #
    # getDensity - Get density for given temperature/pressure
    #
    def getDensity(self, pressure, temperature):
        return self.densityfn(temperature, pressure)

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
    def loadViscosityData(self):
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

        self.viscosityfn = interp2d(pressures, temps, viscosities, kind='linear', fill_value='-1')


    #
    # getViscosity - Get viscosity for given temperature/pressure
    #
    def getViscoity(self, pressure, temperature):
        return self.viscosityfn(pressure, temperature)


    #
    # Load thermal conductivity data from file into memory structures
    #
    # Example usage
    #
    #  Instantiating the class loads this data
    #
    #
    #  See interCond function below for example of how to interogate this data to get thermal conducitivies at other temperatures/pressures
    #
    def loadConductivityData(self):
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

        self.cond_k = k
        self.cond_temps = temps





    #
    # Routine used to interpolate thermal conducvitiy data structures to get value
    # Couldn't use normal scipy routines unfortunately.
    #
    # Example usage
    #
    #     import rp1themoClass as rp1
    #
    #     rp = rp1.rp1thermo()
    #     pressure = 4.0
    #     temperature = 460
    #     print("Finding thermal conductivity for Pressure: ", pressure, " MPa, Temperature: ", temperature, " Kelvin")
    #     kval = rp.interCond(pressure, temperature)
    #     print("Final thermal conductivity: ", kval)
    #
    # Range of data
    # Temperatures from 300 to 645K
    # Pressures    from 0.17MPa to 63 MPa
    #
    def getConductivity(self, pressure, temperature):
        max_temp = np.max(self.cond_temps)
        min_temp = np.min(self.cond_temps)
        if temperature > max_temp or temperature < min_temp:
            print("Temperature (", temperature, ") is outside bounds of data")
            return
        lower_temp = 0
        upper_temp = 0
        lower_idx = 0
        upper_idx = 0
        idx = 0
        for temp in self.cond_temps:
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
            pressures = list(self.cond_k[lower_idx].keys())
            conductivities = list(self.cond_k[lower_idx].values())

            f = interpolate.interp1d(pressures, conductivities)
            lower_kval = f(pressure)

            # Upper Temp side
            pressures = list(self.cond_k[upper_idx].keys())
            conductivities = list(self.cond_k[upper_idx].values())
            f = interpolate.interp1d(pressures, conductivities)
            upper_kval = f(pressure)

            # Interpolate between the two temperatures, to give us the final result
            f = interpolate.interp1d([lower_temp, upper_temp], [lower_kval, upper_kval])
            kval = f(temperature)

            return kval

        else:
            print("Unable to get temperature bounded")
            return