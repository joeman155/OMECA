#!/usr/bin/python3
# -*- coding: latin-1 -*-

#
# Example of how to get thermodynamic properties for RP-1
#


import rp1themoClass as rp1


rp = rp1.rp1thermo()

print("Viscosity is: ", rp.getViscoity(5,360))

print("Density is: ", rp.getDensity(5,360))

print("Conductivity is: ", rp.getConductivity(5, 360))