# Initialize and import packages
import numpy as np
import scipy.interpolate as intp

# I BELIEVE this is channel dimensions.
xHeight = np.array([0,  9,  11, 13, 15, 16, 18, 20, 30])*1e-2
Height = np.array([ 0.8,0.8,0.6,1.0,3.0,1.0,0.4,1.1,2])*1e-3


# Read nozzle coordinates
cont = np.genfromtxt("nozzleContour.csv",delimiter=",")


xVals = cont[0,::]
yVals = cont[1,::]

def interpol(x,y,xNew,how="linear"):
    f = intp.interp1d(x,y,kind=how)
    return f(xNew)


print("Starting")


# Start calculation loop from end of nozzle towards combustion chamber
for i in range(1,len(xVals)):
    # print("i = ", i)
    channelHeight = interpol(xHeight,Height,xVals[-i])   
    print("At X: ", xVals[-i]," ChannelHeight: ", channelHeight)
