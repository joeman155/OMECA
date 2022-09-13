import CEA

# Read CEA file to find adiabatic wall temperature and convective coefficient
CEAfile = "dataCea/methalox.out"
AreaCEA,pCEA,TCEA,rhoCEA,MCEA,muCEA,cpCEA,PrCEA,gCEA,pH2O,pCO2,cstar = CEA.read(CEAfile)



CEAval_curr = len(AreaCEA)-1

print("Length AreaCEA: ", len(AreaCEA))
print("Length pCEA: ", len(pCEA))
print("Length TCEA: ", len(TCEA))
print("Length rhoCEA: ", len(rhoCEA))
print("Length MCEA: ", len(MCEA))
print("Length muCEA: ", len(muCEA))
print("Length cpCEA: ", len(cpCEA))
print("Length PrCEA: ", len(PrCEA))
print("Length gCEA: ", len(gCEA))
print("Length pH2O: ", len(pH2O))
print("Length pCO2: ", len(pCO2))
print("cstar: ", cstar)

print(str("AreaRatio").ljust(10), "   ", str("Press").ljust(8), "   ", str("Temp").ljust(8), "   ", str("Density").ljust(8), "   ", str("Mach").ljust(8), "   ", str("Viscosity").ljust(8), "   ", str("Heat Cap").ljust(8), "   ", str("Prandtl").ljust(8), "   ", str("gamma").ljust(8), "   ", str("PR H2O").ljust(8), "   ", str("PR CO2").ljust(8))
for i in range(1, CEAval_curr):
#    print("i: ", i)
    print(str(AreaCEA[i]).ljust(10), "   ", str(round(pCEA[i])).ljust(8), "   ", str(TCEA[i]).ljust(8), "   ", str(rhoCEA[i]).ljust(8), "   ", str(MCEA[i]).ljust(8), "   ", str(round(muCEA[i], 6)).ljust(8), "   ", str(round(cpCEA[i], 1)).ljust(8), "   ", str(PrCEA[i]).ljust(8), "   ", str(gCEA[i]).ljust(8), "   ", str(round(pH2O[i])).ljust(8), "   ", str(round(pCO2[i])).ljust(8))

