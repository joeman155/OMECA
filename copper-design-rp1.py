#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:59:24 2015

@author: luka
"""

# Initialize and import packages
import numpy as np
import CEAbly as CEA
import rp1themoClass as th
import scipy.interpolate as intp
import matplotlib.pyplot as plt

# Calculate fuel flow from engine parameters
mflow   = 0.89229  # kg/s
OFratio = 2.00
fFlow   = 1 / (1 + OFratio) * mflow


rhoratio = 1.00
kratio   = 1

# Determined that is we multiple these thermal conductivities by 0.7 , each, we get close to what RPA predicts.
hgratio  = 1
hcratio  = 1


# Define nozzle material and thickness
tChamber = 3.0e-3  # wall thickness
kChamber = 295  # W/(m2 K)
rhoChamber = 9134
mu = 0.34  # Poisson's ratio
E = 85e9  # Modulus of elasticity
s_yield = 120932000  # Yield strength

# Define channel geometry
# Not sure where this is defined...throat? Exit? Combustion chamber?
# I am sure this geometry changes.


# We probably want to set the tRib ...from a manufacturing point of view.
# tRib = Rnozzle * 2 * np.pi / NChannels - channelWidth = 
# NChannels = (Rnozzle + tChamber) * 2 * np.pi / (tRib + channelWidth) = (10 + 3) * 2 * 3.1415 / (1.5 + 1.5) = 27.2263 - Close enough to 27.
NChannels = 26  # Number of channels
tRib      = 1.5e-3  # Thickness of Rib
channelHeight = 3e-3  # In RPA, the channel height varies between 2.5 and 3mm  # THIS SETTING IS NOT USED IN THIS CODE. IT IS OVERWRITTEN
roughness = 6e-6  # Not sure what this is, but assume it is right.

# Initialize coolant pressure and temperature
p = pin = 50e5  # Pressure in Pascals
T = Tin = 302  # Temperature of coolant at nozzle exit in Kelvin

# Read nozzle coordinates
cont = np.genfromtxt("nozzleContourbly.csv", delimiter=",")


# Define function for radius of curvature based on coordinates of 3 points
def radiusCurvature(x1, y1, x2, y2, x3, y3):
    num = np.sqrt(
        ((x2 - x1) ** 2 + (y2 - y1) ** 2) * ((x2 - x3) ** 2 + (y2 - y3) ** 2) * ((x3 - x1) ** 2 + (y3 - y1) ** 2))
    den = 2 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2)
    if den == 0:
        return 0
    return num / den


xVals = cont[0, ::]
yVals = cont[1, ::]
# Define engine size (throat radius and area)
rt = 0.01  # radius of throat of nozzle in meters. Later on, we will calculate this from cont
At = rt ** 2 * np.pi
aRatioMinm = min(yVals ** 2 / rt ** 2)


def interpol(x, y, xNew, how="linear"):
    f = intp.interp1d(x, y, kind=how)
    return f(xNew)


# CHANNEL HEIGHT DIMENSIONS AT VARIOUS STATIONS
# It only allows for varying heights, not varying widths. I think we should allow for varying widths.
xHeight = np.array([0, 9, 11, 13, 15, 16, 18, 20, 30]) * 1e-2
# Height = np.array([1.5, 1.5, 2.0, 2.3, 3.0, 4.0, 3, 2.5, 2.5]) * 1e-3
# Height = np.array([1.5, 2.0, 2.2, 2.3, 2.3, 2.0, 1.8, 1.6, 1.8]) * 1e-3
# Height = np.array([2.5, 2.0, 2.2, 2.3, 2.3, 2.0, 1.8, 1.6, 1.8]) * 1e-3
# Height = np.array([2.5, 2.0, 2.2, 2.3, 2.3, 2.0, 1.8, 1.6, 1.2]) * 1e-3

# Level = height .... the same
Height = np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]) * 1e-3
Height = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]) * 1e-3

# big peak
# Height = np.array([2.5, 2.5, 2.5, 2.5, 3.5, 2.5, 2.5, 2.5, 2.5]) * 1e-3


# Height = np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.5]) * 1e-3
# Height = np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.2, 1.5]) * 1e-3
# Height = np.array([3.5, 3.5, 4.0, 4.3, 5.0, 6.0, 5, 4.5, 4.5]) * 1e-3
# Height   = np.array([1,   1,   1.3, 2.0, 2.6, 3.0, 3, 2.0, 1.5]) * 1e-3

# Check for inward buckling (due to coolant pressure)
l = max(xVals)
r = max(yVals)
t = tChamber

gamma = 1
pcrit = 0.855 * E * np.sqrt(gamma) / ((1 - mu ** 2) ** (3. / 4.) * (r / t) ** (5. / 2.) * (l / r))
if pcrit > pin:
    print("Buckling pressure okay: Pin(", round(pin / 1e5), ") is less than critical pressure: ", round(pcrit / 1e5),
          "bar")
else:
    print("**DANGER**: Buckling pressure exceeded. Pin(", round(pin / 1e5), ") excceds ", round(pcrit / 1e5), "bar")

# Check for hoop stress (due to chamber pressure)
s_h = pin * r / t
if s_h < s_yield:
    print("Hoop stress okay. Max Stress experienced (", round(s_h / 1e6), ") is less than yield stress: ",
          round(s_yield / 1e6), "MPa")
else:
    print("**DANGER**: Hoop stress exceeded. Max Stress experienced: ", round(s_h / 1e6), " exceeds yield stress: ",
          round(s_yield / 1e6), "MPa")

# Read CEA file to find adiabatic wall temperature and convective coefficient
CEAfile = "dataCea/rp-1lox.out"
AreaCEA, pCEA, TCEA, rhoCEA, MCEA, muCEA, cpCEA, PrCEA, gCEA, pH2O, pCO2, cstar, velCEA = CEA.read(CEAfile)
T0 = TCEA[0]
p0 = pCEA[0]
# Create class with methane thermophysical model
rp1 = th.rp1thermo()

# Initialize empty lists for parameter storage
pvals = []
Tvals = []
p0vals = []
T0vals = []
V0vals = []
rhovals = []
Tgvals = []
Twvals = []
Twchvals = []
hcvals = []
hgvals = []
wvals = []
Revals = []
Nuvals = []
qvals = []
qradvals = []
Tawvals = []
Civals = []
Dhvals = []
Prvals = []
muvals = []
kapvals = []
cpvals = []
channelHeightvals = []
# Set channel wall temperature to coolant inlet temperature for 1st channel
TwChannel = T
# Pointer to indicate what the current CEA station is, start at nozzle end
CEAval_curr = len(AreaCEA) - 1
# Start channel length at nonzero value to make sure Taylor equation does not crash
x = 0.01
Q = 0
Atot = 0
V = 0
rho = rp1.getDensity(p, T)
rho = rho / rhoratio
cp = rp1.getCp()
Tw = 400
mTot = 0
# Start calculation loop from end of nozzle towards combustion chamber
for i in range(1, len(xVals)):

    # Calculate lenght of channel part and geometry of chamber
    l = np.sqrt((xVals[-i - 1] - xVals[-i]) ** 2 + (yVals[-i - 1] - yVals[-i]) ** 2)
    Rnozzle = yVals[-i]
    aRatio = yVals[-i] ** 2 / rt ** 2

    channelHeight = interpol(xHeight, Height, xVals[-i])
    channelHeightvals.append(channelHeight)
    # Calculate channel cross-sectional dimensions at this nozzle station
    if NChannels == 1:
        A = np.pi * ((Rnozzle + channelHeight) ** 2 - Rnozzle ** 2)
        Dh = th.Dh_shell(Rnozzle + channelHeight, Rnozzle)
    else:
        channelWidth = (tChamber + Rnozzle) * 2 * np.pi / NChannels - tRib
        if channelWidth < 0:
            print("Error: channel width smaller than 0")
        A = NChannels * channelWidth * channelHeight
        Dh = th.Dh_rect(channelWidth, channelHeight)


    # COOLANT: Calculate dynamic pressure and temperature at previous station
    dynPres1 = 0.5 * rho * V ** 2
    dynTemp1 = 0.5 * V ** 2 / cp

    # COOLANT: Calculate density and flow velocity
    rho = rp1.getDensity(p, T) 
    rho = rho / rhoratio
    V = fFlow / (A * rho)

    # COOLANT: Calculate/update static pressure and temperature
    dynPres2 = 0.5 * rho * V ** 2
    p = p - (dynPres2 - dynPres1)
    dynTemp2 = 0.5 * V ** 2 / cp



    T = T - (dynTemp2 - dynTemp1)
    # COOLANT: Calculate thermodynamic properties of methane at current (rho,T)

    mu = rp1.getViscoity(p, T) * rho / 1000000
    cp = rp1.getCp()
    kap = kratio * rp1.getConductivity(p, T)

    # COOLANT:  Calculate bulk flow properties of coolant
    Re = V * rho * Dh / mu
    Pr = mu * cp / kap

    # COOLANT RELATED: Correct for curvature of channel alongside nozzle
    if i > 1 and i < len(xVals):
        (x1, y1) = (xVals[-i - 1], yVals[-i - 1])
        (x2, y2) = (xVals[-i], yVals[-i])
        (x3, y3) = (xVals[-i + 1], yVals[-i + 1])
        Rc = radiusCurvature(x1, y1, x2, y2, x3, y3)
        # Use Niino's formula
        Ci = (Re * (Dh / 4 / abs(Rc)) ** 2) ** (np.sign(Rc) * 0.05)
        # If radius is too high, set correction to 1 (no correction)
        if abs(Rc) > 1:
            Ci = 1
            Rc = 1e9
    else:
        Ci = 1
        Rc = 1e9
    ksi = th.frictionFactor(Dh, roughness, Re) / th.frictionFactor(Dh, 0, Re)
    Cksi = (1 + 1.5 * Pr ** (-1. / 6.) * Re ** (-1. / 8.) * (Pr - 1)) * ksi / (
                1 + 1.5 * Pr ** (-1. / 6.) * Re ** (-1. / 8.) * (Pr * ksi - 1))

    # Check if CEA station should be shifted, depending on current area ratio
    while (not ((aRatio >= AreaCEA[CEAval_curr - 1] and aRatio <= AreaCEA[CEAval_curr]) or (aRatio <= AreaCEA[CEAval_curr - 1] and aRatio >= AreaCEA[CEAval_curr]))):
       CEAval_curr = CEAval_curr - 1




    # HOT GASES: Calculate hot gas parameters depending on CEA values
    pWater = CEA.interpol(aRatio, AreaCEA, CEAval_curr, pH2O)
    pCarbDiox = CEA.interpol(aRatio, AreaCEA, CEAval_curr, pCO2)
    Tg = CEA.interpol(aRatio, AreaCEA, CEAval_curr, TCEA)
    Mg = CEA.interpol(aRatio, AreaCEA, CEAval_curr, MCEA)
    gg = CEA.interpol(aRatio, AreaCEA, CEAval_curr, gCEA)
    Prg = CEA.interpol(aRatio, AreaCEA, CEAval_curr, PrCEA)
    cpg = CEA.interpol(aRatio, AreaCEA, CEAval_curr, cpCEA)
    mug = CEA.interpol(aRatio, AreaCEA, CEAval_curr, muCEA)
    rhog = CEA.interpol(aRatio, AreaCEA, CEAval_curr, rhoCEA)
    Pg   = CEA.interpol(aRatio, AreaCEA, CEAval_curr, pCEA)
    Vg   = CEA.interpol(aRatio, AreaCEA, CEAval_curr, velCEA) 
    Taw = th.adiabatic_wall(Tg, gg, Mg, Prg)



    # HOT GASES: Increase TwNew to avoid missing loop
    TwNew = Tw + 10
    TwChannelNew = TwChannel + 10

    while (abs(TwNew - Tw) > 0.1) and (abs(TwChannel - TwChannelNew) > 0.1):
        Tw = TwNew
        TwChannel = TwChannelNew
        # HOT GASES
        # Calculate convective coefficient using Bartz
        hg = th.bartz(T0, Tw, p0, Mg, rt * 2, aRatio, mug, cpg, Prg, gg, cstar)
        # hgg = th.bartz2(T0, Tw, p0, Mg, rt * 2, aRatio, mug, cpg, Prg, gg, cstar, Rnozzle)
        # print("Compare hg: ", hg , " with hgg: ", hgg)
        # hg = hg / 0.026 * 0.0195
        hg = hg * hgratio
        

        # COOLANT
        # Calculate Nusselt number
        # Nu = th.Taylor(Re, Pr, T, TwChannel, Dh, x)
        Nu = th.dittusBoelter(Re,Pr)
        f = th.frictionFactor(Dh, roughness, Re)
        # Nug = th.gnielinski(Re, Pr, Dh, f)


        # rhow = methane.eqState(p,TwChannel)
        # Nu = th.Ruan(Re,Pr,rho,rhow,Dh,x)
        # Apply correction to Nusselt number
        # Nu = Nu * Ci * Cksi
        # Calculate coolant convective coefficient
        hc = Nu * kap / Dh
        hc = hc * hcratio



        # joe
        # EXPERIMENTAL CODE
        s = 0.001           # Gap size of 1 mm
        Vcool = Vg          # Velocity of liquid film = velocity of hot gas
        Pcool = Pg          # Pressure is the same as hot gases
        # rhoc = rp1.getDensity(Pcool, Tcool) 
        # Recool = Vcool * rhocool * distance / mucool
        # Prcool = mucool * cpcool / kapcool
        # St = Nu / Re / Pr

        # Location of orifice
        # xloc = 0.02
        # x = xVals[-i - 1] - xloc   # How far away we are from the orifice

        # rho = density of coolant in channel, not in the film... but we can fix this up later
        # 
        # F = rho * vcool / (rhog * v)

        # eff = th.coolingEffiency(St, x, F, s, Recool, Prcool, v, vcool)
        # eff = th.adbiaticCoolingEff(3.09, x, vc, v, Recool)
        # print("eff = ", eff)






        # Incorporate heat transfer fin effectiveness into hc (Heat Transfer coefficient for hot gases)
        m = np.sqrt(2 * hc * tRib / kChamber)

        # Calculate finEffectiveness, a number between 0 and 1 We want as close as possible.
        finEffectiveness = np.tanh(m / tRib * channelHeight) / (m / tRib * channelHeight)

        # Calculate new hc.
        hc = hc * (channelWidth + finEffectiveness * 2 * channelHeight) / (channelWidth + tRib)
        hcboost = (channelWidth + finEffectiveness * 2 * channelHeight) / (channelWidth + tRib)

        # Calculate radiative heat transfer (Suspect this is for the MAIN products of combustion)
        qW = 5.74 * (pWater / 1e5 * Rnozzle) ** 0.3 * (Taw / 100) ** 3.5  # Water
        qC = 4 * (pCarbDiox / 1e5 * Rnozzle) ** 0.3 * (Taw / 100) ** 3.5  # Carbon Dioxide
        qRad = qW + qC
        # qRad = 0

        # Calculate heat flux
        q = (Taw - T + qRad / hg) / (1 / hg + tChamber / kChamber + 1 / hc)

        effk = 1 / (1 / hg + tChamber / kChamber + 1 / hc)
        # print("effective heat coeff = ", effk, ", q  = ", q, ", Taw = ", Taw, ", hg = ", hg, ", kChamber/tChamber = ", kChamber/tChamber, ", hc = ", hc)

        # Calculate hot gas wall temperature and channel wall temperature
        TwNew = Taw - (q - qRad) / hg
        TwChannelNew = T + q / hc
        # These get fed back into top at beginning of While loop.
        # print("Comparing TwNew: ", TwNew , " with Tw: ", Tw)
        # print("Comparing TwChannelNew ", TwChannelNew , " with TwChannel: ", TwChannel)

    Tw = TwNew
    TwChannel = TwChannelNew


    xx = np.round(xVals[-i] * 100, 2)
    print("")
    print("")
    print("")
    print("X POSITION: ", xx, " mm")
    print("---------------------------------")
    print("channelWidth = ", np.round(1000 * channelWidth, 2), ", ChannelHeight = ", np.round(1000 * channelHeight, 2), ", Area = ", np.round(1000000 * A, 2))
    print("V,rho,Dh,mu = ", np.round(V, 2), np.round(rho, 0), np.round(Dh, 6), np.round(mu, 6))
    print("Coolant Pressure, Temperature = ", np.round(p, 2), np.round(T, 1))
    print("Heat Coefficients: hg = ", np.round(hg, 2), " - hc = ", np.round(hc, 2), ", Nu = ", np.round(Nu, 2), ", RE = ", np.round(Re,3), ", Pr = ", np.round(Pr, 3), ", kap = ", np.round(kap, 2), ", Dh = ", np.round(1000 * Dh, 2), ", FinEffectiveness = ", finEffectiveness, ", hcboost = ", hcboost)
    print("Temperatures:  Adbiatic Gas Temp: ", Taw, " Wall(gas) Temperature: ", Tw, ", Channel Wall Temp: ", TwChannel)
    print("Pressure of hot gases: ", Pg)
    Vgspeed = Vg * Mg
    print("Velocity of hot gases: ", Vgspeed, " ms-1, MACH: ", Mg)
    print("Heat Flux: ", q)

    # Calculate change in temperature and pressure
    A_heat = 2 * np.pi * Rnozzle * l  # Area of station being considered
    deltaT = q * A_heat / (fFlow * cp)  # Change in temperature of coolant
    # deltap = th.frictionFactor(Dh, roughness, Re) * l / Dh * rho * V ** 2 / 2.0  # Change in pressure of coolant

    # We expect flow to be mostly laminar, so we use frictionFactorLaminar to calculate the frictionfactor coefficient.
    # deltap = th.frictionFactorLaminar(Re) * l / Dh * rho * V ** 2 / 2.0  # Change in pressure of coolant
    deltap = f * l / Dh * rho * V ** 2 / 2.0  # Change in pressure of coolant


    Q = Q + q * A_heat  # Change in enthalpy of coolant
    Atot = Atot + A_heat  # NOT SURE WHAT THIS IS...
    mCur = (2 * np.pi * Rnozzle * l * tChamber + l * tRib * channelHeight * NChannels) * rhoChamber  # ?? change in ??
    mTot = mTot + mCur  # ??

    # Update pressure, temperature and channel length
    p = p - deltap  # Update pressure of coolant
    T = T + deltaT  # Update temperature of coolant
    x = x + l  # new X position.

    p0vals.append(p + 0.5 * rho * V ** 2)
    T0vals.append(T + 0.5 * V ** 2 / cp)
    V0vals.append(V)

    # Store parameters in lists
    pvals.append(p)
    Tvals.append(T)
    rhovals.append(rho)
    Tgvals.append(Tg)
    Twvals.append(Tw)
    Twchvals.append(TwChannel)
    Tawvals.append(Taw)
    hcvals.append(hc)
    hgvals.append(hg)
    wvals.append(channelWidth)
    Revals.append(Re)
    Nuvals.append(Nu)
    qvals.append(q)
    Civals.append(Ci)
    Dhvals.append(Dh)
    Prvals.append(Pr)
    muvals.append(mu)
    kapvals.append(kap)
    cpvals.append(cp)

# Print output for user
print(min(wvals) * 1e3, "mm minimum channel width")
print(max(Twvals), "K maximum wall temperature")
print((p0vals[0] - p0vals[-1]) / 1e5, "bar pressure loss.")
print(T - Tvals[0], "K temperature rise")
print(Q, "Total heat input")
print(mTot, "kg chamber mass")




# Plot results


# Create figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(36 / 2.54, 15 / 2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = list(range(7))

# Wall temperature
lins[0] = ax.plot(xVals[1:] * 100, Twvals[::-1], 'g--', lw=2, label=r'$T_w$')
lins[1] = ax.plot(xVals[1:] * 100, Twchvals[::-1], 'y--', lw=2, label=r'$T_cchan$')
lins[2] = ax.plot(xVals[1:] * 100, Tawvals[::-1], 'p--', lw=2, label=r'$T_gadbiatic$')
lins[3] = ax.plot(xVals[1:] * 100, Tgvals[::-1], 'm--', lw=2, label=r'$T_ggas$')
ax.set_ylim([0, np.round(max(Twvals) + 100, 2)])

# Heat flux
ax2 = ax.twinx()
lins[4] = ax2.plot(xVals[1:] * 100, np.array(qvals[::-1]) / 1e7, 'r-.', lw=2, label=r'$q$')

# Geometry
heights = interpol(xHeight, Height, xVals)
lins[5] = ax2.plot(xVals * 100, heights * 1e3, 'b:', lw=2, label=r'$d_c$')
lins[6] = ax2.plot(xVals * 100, yVals * 10, 'k-', label=r'Contour')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines, labs, loc=6, labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Temperature [K]")
ax2.set_ylabel(r"$d_c$ [mm]; $q$ [$\mathrm{10^7 W/m^2}$]; Radius [10 cm]")
ax.set_ylim([400, 3400])
ax2.set_ylim([0, 4])
ax.grid()
plt.show()




# Create figure
fig = plt.figure(2)
fig.clf()
fig.set_size_inches(36 / 2.54, 15 / 2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = list(range(3))

# Reynolds number
lins[0] = ax.plot(xVals[1:] * 100, np.array(Revals[::-1]) / 1e4, 'r--', lw=2, label=r'Re')
# ax.set_ylim([0,round(max(Revals)+100,2)])

# Nusselt number
ax2 = ax.twinx()
lins[1] = ax2.plot(xVals[1:] * 100, Nuvals[::-1], 'b-.', lw=2, label=r'Nu')

# Nozzle contour
lins[2] = ax.plot(xVals * 100, yVals * 100, 'k-', label=r'Contour')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines, labs, loc=0, labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Radius [cm]; Re [$\mathrm{10^4}$]")
ax2.set_ylabel(r"Nu [-]")
ax.set_ylim([0, 80])
ax2.set_ylim([0, 200])
ax.grid()
plt.show()

# Create figure
fig = plt.figure(3)
fig.clf()
fig.set_size_inches(36 / 2.54, 15 / 2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = list(range(3))

# Pressure
lins[0] = ax.plot(xVals[1:] * 100, np.array(p0vals[::-1]) / 1e5, 'b--', lw=2, label=r'$p_0$')
lins[1] = ax.plot(xVals[1:] * 100, V0vals[::-1], 'b:', lw=2, label=r'$V_0$')

# Temperature
ax2 = ax.twinx()
lins[2] = ax2.plot(xVals[1:] * 100, T0vals[::-1], 'r-.', lw=2, label=r'$T_0$')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines, labs, loc=7, labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"$p_{0,c,b}$ [bar]; Velocity [m^-1]")
ax2.set_ylabel(r"$T_{0,c,b}$ [K]")
ax.set_ylim([00, 75])
ax2.set_ylim([000, 1000])
ax.grid()
plt.show()

# Create figure
fig = plt.figure(4)
fig.clf()
fig.set_size_inches(36 / 2.54, 15 / 2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = list(range(1))

# Pressure
lins[0] = ax.plot(xVals * 100, yVals * 100, 'k-', label=r'Contour')
ax.plot(xVals * 100, -yVals * 100, 'k-')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines, labs, loc=7, labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Radius [cm]")
ax.grid()
plt.show()
