import plotting
import processing

"""

The inputs are:
x,
y,
z,
Initial Ball Speed,
vertical release angle,
horizontal release angle,
Spin Rate,
Tilt in Hours,
Tilt in minutes,
Spin Efficiency,
Y seam orientation angle,
Z seam orientation angle,




Primay inputs are: initial position, x0, y0, and z0 with origin at the
point of home plate, x to the rright of the catcher, y from the catcher
towards the pitcher, and z straight up. Initial velocities
u0, v0, and w0 which are the speeds of the ball in x, y, and z
respectivley. And spin rates

UMBA1.0: This code uses a constant Cd and rod cross's model for CL
Predictions. Seam Orientation is not accounted for. Air Density is
considered only at sea level at 60% relative humidity. but can be easily
altered

UMBA2.0 Adding seam positions and attempting to model CL from seams.
"""
#def runAllCode(
pX = []
pY = []
pZ = []
IX = []
IY = []
IZ = []
DX = []
DY = []
DZ = []
FX = []
FY = []
FZ = []
VariableS = []
#num = []

variable = [-38, 0, 38]
iters = len(variable)
for i in range(iters):
    print(i, 'out of ',iters)

    """

    When running a single case you can imput the variables you want to use
    below in the function argument and put a zero in the "Variable" array above

    """
    # all the variables below are at release
    x = 0 #initial x ft
    y = 5.5 #initial y location  ft
    z = 6 #initial z location ft
    v = 75 #initial velocity  mph
    vRelAng = 15 #initial vertical release angle deg
    hRelAng = 0 #initial horizaontal release angle deg
    rpm = 1200 #initial spin rate
    tHrs = 3 #initial Tilt in hours
    tMin = 0 #initial tilt mins
    eff = 100 #efficiency as a percentage
    Yang = 0 #initial seam orientation y angle
    Zang = variable[i] #initial seam orientation z angle
    LorR = 'r' #if efficiency is less than 100% which poll should be forward l or r

    positions = processing.umba(x,y,z,v,vRelAng,hRelAng,rpm,tHrs,tMin,eff,Yang,Zang,LorR, i)

    pX.append(positions[0])
    pY.append(positions[1])
    pZ.append(positions[2])
    IX.append(positions[3])
    IY.append(positions[4])
    IZ.append(positions[5])
    DX.append(positions[6])
    DY.append(positions[7])
    DZ.append(positions[8])
    FX.append(positions[9])
    FY.append(positions[10])
    FZ.append(positions[11])
    VariableS.append(variable[i])
#    num.append(i)

plotting.plotSFinalIter(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,iters,variable)

