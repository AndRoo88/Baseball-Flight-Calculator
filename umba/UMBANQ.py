import numpy as np
import plotting
import processing


def umba(x,y,z,Vtot,Theta, Psi, SpinRate, TiltH, Tiltm, SpinE, Yang, Zang, LorR):
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

    """

    """
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

    Yang = (Yang) * np.pi/180
    Zang = -Zang * np.pi/180
    i = 0

    seamsOn = True
    frameRate = 0.002



    Tilt = processing.TimeToTilt(TiltH, Tiltm)
    if LorR == 'l':
        Gyro = np.arcsin(SpinE/100)
    elif LorR =='r':
        Gyro = np.pi - np.arcsin(SpinE/100)
    else:
        while LorR != 'l' or LorR != 'r':
            if LorR == 'l':
                Gyro = np.arcsin(SpinE/100)
            elif LorR =='r':
                Gyro = np.pi - np.arcsin(SpinE/100)
            else:
                LorR = input('please type in an "l" or an "r" for which pole goes forward')

    positions = (processing.PitchedBallTraj(x, y, z, Vtot, Theta, Psi, SpinRate, Tilt, Gyro, Yang, Zang, i, frameRate, seamsOn))
    plotting.Plotting(positions)

    pX = (positions[0])
    pY = (positions[1])
    pZ = (positions[2])
    IX = (positions[3])
    IY = (positions[4])
    IZ = (positions[5])
    DX = (positions[6])
    DY = (positions[7])
    DZ = (positions[8])
    FX = (positions[9])
    FY = (positions[10])
    FZ = (positions[11])

    return(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ)


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

variable = [100, 90]
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
    v = 90 #initial velocity  mph
    vRelAng = 0 #initial vertical release angle deg
    hRelAng = 0 #initial horizaontal release angle deg
    rpm = 1200 #initial spin rate
    tHrs = 3 #initial Tilt in hours
    tMin = 0 #initial tilt mins
    eff = variable[i] #efficiency as a percentage
    Yang = 0 #initial seam orientation y angle
    Zang = 38 #initial seam orientation z angle
    LorR = 'l' #if efficiency is less than 100% which poll should be forward l or r

    positions = umba(x,y,z,v,vRelAng,hRelAng,rpm,tHrs,tMin,eff,Yang,Zang,LorR)

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

