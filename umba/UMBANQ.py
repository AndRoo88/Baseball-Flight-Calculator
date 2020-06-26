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
TF = []
aX = []
aY = []
aZ = []
VariableS = []
#num = []

variable = [100]
iters = len(variable)
for i in range(iters):
    print(i, 'out of ',iters)

    """

    When running a single case you can imput the variables you want to use
    below in the function argument and put a zero in the "Variable" array above

    """
    # all the variables below are at release
    x = 0 #initial x location ft
    y = 6.5 #initial y location  ft
    z = 50 #initial z location ft
    v = 10 #initial velocity  mph
    vRelAng = -0 #initial vertical release angle deg
    hRelAng = 0.0 #initial horizaontal release angle deg
    rpm = 1300 #initial spin rate
    tHrs = 3 #initial Tilt in hours
    tMin = 0 #initial tilt mins
    eff = variable[i] #efficiency as a percentage
    Yang = 0 #initial seam orientation y angle
    Zang = -20 #initial seam orientation z angle
    LorR = 'r' #if efficiency is less than 100% which poll should be forward l or r
    seamsOn = True
    FullRot = False

    positions = processing.umba(x,y,z,v,vRelAng,hRelAng,rpm,tHrs,tMin,eff,Yang,Zang,LorR, i, seamsOn, FullRot)

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
    TF.append(positions[12])
    aX.append(positions[13])
    aY.append(positions[14])
    aZ.append(positions[15])
    VariableS.append(variable[i])


#plotting.AccelPlots(aX,aY,aZ,TF,iters,variable)
if len(variable) > 1:
    plotting.plotSFinalIter(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,iters,variable)
#plotting.AccelPlots(pX,pY,pZ,TF,j,variable)
