from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math

fig, axes = plt.subplots()

subcycleNumber = 7680

operatorMethods = ["wachspress","pwl","weak"]

for operatorMethod in operatorMethods:

    # data in
    filenameIn = "./output_hex_%s_%i/output.2000.nc" %(operatorMethod,subcycleNumber)
    filein = Dataset(filenameIn, "r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])
    nTimes = len(filein.dimensions["Time"])

    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex -= 1

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    uVelocity = filein.variables["uVelocity"][-1,:]
    vVelocity = filein.variables["vVelocity"][-1,:]

    uVelocities = filein.variables["uVelocity"][:,:]

    filein.close()

    xmin = np.amin(xVertex)
    xmax = np.amax(xVertex)
    ymin = np.amin(yVertex)
    ymax = np.amax(yVertex)

    us = []
    for iTime in range(0,nTimes):

        x = []
        u = []
        for iVertex in range(0,nVertices):
            if (math.fabs(yVertex[iVertex] - 508068.236886871) < 1e-8):
                x.append(xVertex[iVertex])
                u.append(uVelocities[iTime,iVertex])
        x = np.array(x)
        u = np.array(u)

        sortedIdxs = x.argsort()

        x = x[sortedIdxs]
        u = u[sortedIdxs]

        us.append(math.sqrt(np.sum(np.power(u,2))))

        if (iTime == nTimes-1):
            axes.plot(x, u, label=operatorMethod)

#axes.plot(x, np.zeros(x.shape[0]), zorder=1, c='k')

uAir = 1.0
rhoair = 1.3
rhow = 1026.0
cocn = 0.00536
cair = 0.0012
Pstar = 2.75e4
Cstar = 20.0
e = 2
alpha = math.sqrt(1.0 + math.pow(1.0 / e, 2))
Lx = 1280000


uu = []
for xx in x:

    a = xx / Lx
    v = 2.0 * a
    dadx = (1.0 / Lx)
    dvdx = 2.0 * dadx

    oceanStressCoeff = rhow * cocn * a
        
    airStress = rhoair * uAir * uAir * a * cair

    P = Pstar * v * math.exp(-Cstar * (1-a))

    dPdx = Pstar * math.exp(-Cstar * (1-a)) * (dvdx + v * Cstar * dadx)

    print(xx, a, -Cstar * (1-a), P, dPdx)
        
    u = max((airStress - 0.5*(alpha + 1.0) * dPdx) / oceanStressCoeff, 0.0)
    uu.append(u)

axes.plot(x, uu, zorder=2, c='r')

axes.set_xlabel("time")
axes.set_ylabel("uVelocity")
axes.legend()

plt.savefig("1D_velocity_operator.png",dpi=300)
