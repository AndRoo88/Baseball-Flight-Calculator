import numpy as np
from scipy.spatial.transform import Rotation as R
import math
import plotting
#import time


def normalPlane(VelVec,SpinVecT,t):
    """
    finds an acceptable range of seam locations where they can have an affect
    on the aerodynamics.

    Plans: the code finds the acceptable range wherein the seams can produce a
    force on the ball normal to the direction of flight. the force the ball
    will also be normal to the surface of the ball at that location.
    However, htis portion of code is only to find the range.
    """

    SpinVecTNew = np.zeros([3])
    SpinVecTMag = np.linalg.norm(SpinVecT)
    SpinVecTUnit = SpinVecT/SpinVecTMag
    dt = 0.001
    SpinShiftFactor = 1.25 #this number effects how much the separation location
                            #Will change based on the spin rate. Bigger, Move shift
    forwardBackV = 0. # allows for the moving the effectiveness of the seams
    # forwards or backwards.

    AngleOfActivation = 5

    diameter = (2. + 15/16) #in
    acceptableRange = diameter*1.1
    acceptableThickness = diameter/2*np.sin(AngleOfActivation*2*np.pi/180)
#
    #produces a 3d box
    xmin0 = -acceptableRange*.5
    xmax0 =  acceptableRange*.5
    zmin0 = -acceptableRange*.5
    zmax0 =  acceptableRange*.5
    ymin0 = -acceptableThickness + forwardBackV
    ymax0 =  acceptableThickness + forwardBackV

    node10 = [xmin0,ymin0,zmin0]
    node20 = [xmin0,ymax0,zmin0]
    node30 = [xmax0,ymax0,zmin0]
    node40 = [xmax0,ymin0,zmin0]
    node50 = [xmin0,ymin0,zmax0]
    node60 = [xmin0,ymax0,zmax0]
    node70 = [xmax0,ymax0,zmax0]
    node80 = [xmax0,ymin0,zmax0]
    nodes0 = np.asarray([node10, node20, node30, node40, node50, node60, node70, node80])

    VRotVec = findRotVec(0, -1., 0, VelVec[0], VelVec[1], VelVec[2])
    SpinVecTNew[0] = SpinVecT[0] * SpinShiftFactor * dt
    SpinVecTNew[1] = SpinVecT[1] * SpinShiftFactor * dt
    SpinVecTNew[2] = SpinVecT[2] * SpinShiftFactor * dt

    r = R.from_rotvec(VRotVec)
    rr = R.from_rotvec(SpinVecTNew)
    nodes1 = r.apply(nodes0)
    nodes2 = rr.apply(nodes1)
#    print(nodes1-nodes2)

#    x0check = nodes0[0,0]
#    y0check = nodes0[0,1]
#    v1 = [x0check,y0check]
#
#    x2check = nodes2[0,0]
#    y2check = nodes2[0,1]
#    v2 = [x2check,y2check]

#    cosang = np.dot(v1, v2)
#    sinang = np.linalg.norm(np.cross(v1, v2))
#    angle = (np.arctan2(sinang, cosang))
#    print(angle*180/np.pi)
#    print(x0check, y0check, '\n', x2check, y2check)

#    print(np.arctan2(nodes2[1],nodes2[0])*180/np.pi)
    return (nodes2)

###############################################################################

def derivs(t, BallState, BallConsts, activeSeams):
    """
    This is where the magic happens all models are input here
    Ball State:

    1, x
    2, y
    3, z
    4, u
    5, v
    6, w
    7, spinx
    8, spiny
    9, spinz
    """

    dy = np.zeros(len(BallState))

    u = BallState[3]
    v = BallState[4]
    w = BallState[5]
    Spinx = BallState[6]
    Spiny = BallState[7]
    Spinz = BallState[8]#rad/sec
    VelTot = np.sqrt(u**2 + v**2 + w**2)

    SpinRate = np.sqrt(Spinx**2 + Spiny**2 + Spinz**2)
    diameter = BallConsts[1] #ft
    c0 = BallConsts[4]

    rw = (diameter/2)*SpinRate
    S = (rw/VelTot)*np.exp(-t/10000) #the "np.exp(-t/NUM) is for spin decay
    #for no spin decay NUM should be large. When better data is available on
    #spin decay will account for it here likely

    Cl = ClCross(S)
#    Cl = ClKensrud(S) This is not right yet
    CdConst = 0.33

    # The coefficient of Seams "Cseams" is the essentially the Lift coeficient
    # per seam per length away from the origin.
    Cseams = .02 #per active seam

    aDragx = -c0*CdConst*VelTot*u
    aDragy = -c0*CdConst*VelTot*v
    aDragz = -c0*CdConst*VelTot*w

    aSpinx = c0*(Cl/SpinRate)*VelTot*(Spiny*w - Spinz*v)
    aSpiny = c0*(Cl/SpinRate)*VelTot*(Spinz*u - Spinx*w)
    aSpinz = c0*(Cl/SpinRate)*VelTot*(Spinx*v - Spiny*u)

    SeamXLength = 0
    SeamYLength = 0
    SeamZLength = 0

    if len(activeSeams) > 4:
        for i in range(len(activeSeams)):
            SeamXLength = SeamXLength + activeSeams[i,0]
            SeamYLength = SeamYLength + activeSeams[i,1]
            SeamZLength = SeamZLength + activeSeams[i,2]

    aSeamsx = -c0*Cseams*(VelTot**2)*SeamXLength
    aSeamsy = -c0*Cseams*(VelTot**2)*SeamYLength
    aSeamsz = -c0*Cseams*(VelTot**2)*SeamZLength


#    print(SeamXLength,SeamYLength, SeamZLength)

    ax = aDragx + aSpinx + aSeamsx
    ay = aDragy + aSpiny + aSeamsy
    az = aDragz + aSpinz + aSeamsz - 32.2

    dSpinx = 0
    dSpiny = 0
    dSpinz = 0

    dy[0] = u
    dy[1] = v
    dy[2] = w
    dy[3] = ax
    dy[4] = ay
    dy[5] = az
    dy[6] = dSpinx
    dy[7] = dSpiny
    dy[8] = dSpinz

    return dy

###############################################################################


def PitchedBallTraj(x,y,z,Vtot, Theta, Psi, SpinRate, Tilt, Gyro, Yangle, Zangle,i, frameRate, seamsOn):


    #this is wehre the work needs to happen now

    FullState = anglesTOCart(Vtot, Theta, Psi, SpinRate, Tilt, Gyro, Yangle, Zangle)
    print(FullState)

    x0 = x
    y0 = 60.5 - y
    z0 = z
    u0 = FullState[0]
    v0 = FullState[1]
    w0 = FullState[2]
    Spinx0 = FullState[3]
    Spiny0 = FullState[4]
    Spinz0 = FullState[5]
    Yangle = FullState[6] #angle 1 is the angle from
    Zangle = -FullState[7]

    Spinx0 = Spinx0 * .104719754 #converts rps to rad/s
    Spiny0 = Spiny0 * -.104719754
    Spinz0 = Spinz0 * .104719754

    xSeam, ySeam, zSeam = initializeSeam() #initialized seams to a 90 deg x rotations
    # see https://www.baseballaero.com/2020/03/09/describing-ball-orientation-post-51/
    # for further info
    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, -np.pi/2, 0, 0, 1)
    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, 0, np.pi/2, 0, 1)
    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, 0, Yangle, 0, 1)
    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, 0, 0, Zangle, 1)
#    xSeam3, ySeam3, zSeam3 = rotateSeam(1, 0, 0, Spinx0, Spiny0, Spinz0, 1)
#    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, 0, Yangle, 0, 1)
    RotVec = findRotVec(1,0,0,Spinx0,Spiny0,Spinz0)

    xSeam, ySeam, zSeam = rotateSeam(xSeam, ySeam, zSeam, RotVec[0],RotVec[1],RotVec[2], 1)
    xSeam0,ySeam0,zSeam0 = xSeam, ySeam, zSeam
#    xSeam0, ySeam0, zSeam0 = rotateSeam(xSeam, ySeam, zSeam, 0, Yangle, 0, 1)
#    RotVec = findRotVec(Yangle, Zangle, 0, Spinx0, Spiny0, Spinz0)
#    xSeam0, ySeam0, zSeam0 = rotateSeam(xSeam2, ySeam2, zSeam2, RotVec[0],RotVec[1],RotVec[2],1)

#    xSeam0, ySeam0, zSeam0 = xSeam, ySeam, zSeam

    # extablished the (0,0) initial condition as a 2-seam spin

    # All air properties are the approximate averages for sea level over the season
#    rhoDRY = 0.0765 #lb/ft^3
#    relHum = 0.73
#    Temp = 60 #deg fahrenheit
    rho = 0.074  #lb/ft^3, with humidity at sea level

    circ = 9.125/12 #ft
    diameter = (2. + 15/16)/12 #ft
    Area = .25*np.pi*diameter**2 #ft^2
    mass = 0.3203125 #lbm
    c0 = 0.5*rho*Area/mass
    BallConsts = [circ,diameter,Area,mass,c0]

    t0 = 0.0
    t = t0
    dt = 0.001

    u0 = u0 *1.467#ft/sec
    v0 = -v0 *1.467#ft/sec
    w0 = w0 *1.467#ft/sec

    decisionPoint = 0.2 #sec #time before ball arrives when batter has
                        #to decide to hit or not.

    SpinVec = [Spinx0,Spiny0,Spinz0]
    Vel = [u0,v0,w0]
    VelTot = np.sqrt(u0**2 + v0**2 + w0**2)
    SpinRate0 = np.sqrt(Spinx0**2 + Spiny0**2 + Spinz0**2)
    SpinEfficiency0 = 1-abs(np.dot(Vel, SpinVec)/(SpinRate0*VelTot))
    #assumes that the efficiency is non-linear and that it follows the sin of the
    #angle between the ball direction and the spin.


    BallState0 = [x0,y0,z0,u0,v0,w0,Spinx0,Spiny0,Spinz0]

    fileBT = open(str(i) + "BallTrajectoryNEW.txt","w+")
    fileBT.write("time        x        y         z          u            v       w      Spin x     Spin y    Spin z\n")
    fileBT.write("==================================================================================================\n")
    fileBT.write("{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}\n"
                 .format(t,x0,y0,z0,u0,v0,w0,Spinx0,Spiny0,Spinz0))

    xP = []
    yP = []
    zP = []
    uP = []
    vP = []
    wP = []
    xD = BallState0[0]
    yD = BallState0[1]
    zD = BallState0[2]
    uD = BallState0[3]
    vD = BallState0[4]
    wD = BallState0[5]
    while BallState0[1] > 0. and BallState0[2] > 0. and t < 20:

        SpinVec = [Spinx0,Spiny0,Spinz0]
        #need to input a non-magnus ball path indicator.
        if seamsOn == True:
            xSeam1, ySeam1, zSeam1 = rotateSeam(xSeam0, ySeam0, zSeam0, BallState0[6],BallState0[7],BallState0[8],dt)
            seamPoints = np.asarray([xSeam1, ySeam1, zSeam1])

            VelVec = np.asarray([BallState0[3], BallState0[4], BallState0[5]])
            activeSeams, inactiveSeams, nodes = findSSWseams(VelVec,seamPoints,SpinVec,t)
#            plotting.plotPointsTest(activeSeams, inactiveSeams, nodes,t)

            xSeam0 = xSeam1
            ySeam0 = ySeam1
            zSeam0 = zSeam1

#            if t == 0:
#                seamPoints0 = np.asarray([xSeam0, ySeam0, zSeam0])
#                activeSeams0, inactiveSeams0, nodes0 = findSSWseams(VelVec,seamPoints0,SpinVec,t)
#                plotting.plotSeams(activeSeams0, inactiveSeams0, Spinx0, Spiny0, Spinz0, 0, VelVec,nodes)
#                time.sleep(10)
            if t % frameRate > -0.0000001 and t % frameRate < 0.0000001 and (SpinRate0*t) < np.pi*.2 and t <= .0501:# and SpinRate0 > 100:
                plotting.plotSeams(activeSeams, inactiveSeams, Spinx0, Spiny0, Spinz0, t, VelVec, nodes)
        else:
            activeSeams = [0,0,0]
#         # This section is for showing the spin behaviour of the ball and where
#         # the seams are moving


        BallState1 = RK4(t, BallState0, dt, BallConsts, activeSeams)

#            plotting.plotPointsTest(activeSeams, inactiveSeams, nodes,t)

        fileBT.write("{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}\n"
                 .format(t,BallState1[0],BallState1[1],BallState1[2],BallState1[3],BallState1[4],BallState1[5],BallState1[6],BallState1[7],BallState1[8]))

        BallState0 = BallState1

        xP.append(BallState1[0]*12)
        yP.append(BallState1[1])
        zP.append(BallState1[2]*12)
        uP.append(BallState1[3])
        vP.append(BallState1[4])
        wP.append(BallState1[5])
        t = t + dt

    DecisionPointStep = int(t - (.2/dt))
    if t < decisionPoint:
        print("WOW! no batter has enough skill to hit a ball thrown that fast")

        xD = -10
        yD = -10
        zD = -10
        uD = 0
        vD = 0
        wD = 0
    else:
        xD = xP[DecisionPointStep]/12
        yD = yP[DecisionPointStep]
        zD = zP[DecisionPointStep]/12
        uD = uP[DecisionPointStep]
        vD = vP[DecisionPointStep]
        wD = wP[DecisionPointStep]

    BallStateF = BallState1
    xF, yF, zF = BallStateF[0], BallStateF[1], BallStateF[2]
    fileBT.close()

    dzNoSpin = w0*t - (32.2/2)*t*t
    zfg = z0 + dzNoSpin
    vBreak = BallStateF[2] - zfg

    dxNoSpin = u0*t
    xfg = x0 + dxNoSpin
    hBreak = BallStateF[0] - xfg

    SpinVecF = [BallStateF[6],BallStateF[7],BallStateF[8]]
    VelF = [BallStateF[3],BallStateF[4],BallStateF[5]]
    VelTotF = np.sqrt(BallStateF[3]**2 + BallStateF[4]**2 + BallStateF[5]**2)
    SpinRateF = np.sqrt(BallStateF[6]**2 + BallStateF[7]**2 + BallStateF[8]**2)
    SpinEfficiencyF = 1-abs(np.dot(VelF, SpinVecF)/(SpinRateF*VelTotF))

    totalRotations = SpinRateF/(2*np.pi) #assumes no spin decay

    finalApproachAngleyz = np.arctan2(abs(BallStateF[5]), abs(BallStateF[4]))
    finalApproachAnglexy = np.arctan2(abs(BallStateF[3]), abs(BallStateF[4]))

    Hrs, mins = TiltToTime(Tilt)
#    Tiltdegs = TimeToTilt(Hrs,mins)


    print('initial conditions:')
    print('x0 (ft)------------------------------- ', to_precision(x0,4))
    print('y0 (ft)------------------------------- ', to_precision(y0,4))
    print('z0 (ft)------------------------------- ', to_precision(z0,4))
    print('u0 (mph)------------------------------ ', to_precision(u0/1.467,4))
    print('v0 (mph)------------------------------ ', to_precision(v0/1.467,4))
    print('w0 (mph)------------------------------ ', to_precision(w0/1.467,4))
    print('Total Velocity (mph)------------------ ', to_precision(VelTot/1.467,4))
    print('Spinx0 (rpm)-------------------------- ', to_precision(Spinx0/0.104719754,4))
    print('Spiny0 (rpm)-------------------------- ', to_precision(Spiny0/-0.104719754,4))
    print('Spinz0 (rpm)-------------------------- ', to_precision(Spinz0/0.104719754,4))
    print('Total Spin Rate (rpm)----------------- ', to_precision(SpinRate0/0.104719754,4))
    print('Tilt (clock face)----------------------', Hrs,':',mins)
#    print('Tilt (deg) --------------------------- ', to_precision(Tiltdegs,4))
    if SpinRate0 == 0:
        print('Initial Efficiency (%)---------------- NA')
    else:
        print('Initial Efficiency (%)---------------- ', to_precision(SpinEfficiency0*100,4))

    print('\n\nconditions at decision point:')
    print('x at decision point (ft)------------- ', to_precision(xD,4))
    print('y at decision point (ft)------------- ', to_precision(yD,4))
    print('z at decision point (ft)--------------', to_precision(zD,4))
    print('u at decision point (ft)--------------', to_precision(uD,4))
    print('v at decision point (ft)--------------', to_precision(vD,4))
    print('w at decision point (ft)--------------', to_precision(wD,4))

    print('\n\nconditions across the plate:')
    print('xf (ft)-------------------------------', to_precision(BallStateF[0],4))
    print('yf (ft)-------------------------------', to_precision(BallStateF[1],4)) # actually just the last point data was taken
    print('zf (ft)-------------------------------', to_precision(BallStateF[2],4))
    print('uf (mph)------------------------------', to_precision(BallStateF[3]/1.467,4))
    print('vf (mph)------------------------------', to_precision(-BallStateF[4]/1.467,4))
    print('wf (mph)------------------------------', to_precision(BallStateF[5]/1.467,4))
    print('Total Velocity (mph)------------------', to_precision(VelTotF/1.467,4))
    print('Spinxf (rpm)--------------------------', to_precision(BallStateF[6]/0.104719754,4))
    print('Spinyf (rpm)--------------------------', to_precision( BallStateF[7]/0.104719754,4))
    print('Spinzf (rpm)--------------------------', to_precision(BallStateF[8]/0.104719754,4))
    print('Total Spin Rate (rpm)-----------------', to_precision(SpinRateF/0.104719754,4))
    print('Approach Angle (yz, deg)--------------', to_precision(finalApproachAngleyz*180/np.pi,4))
    print('Approach Angle (xy, deg)--------------', to_precision(finalApproachAnglexy*180/np.pi,4))
    print('Final Efficiency (%)------------------', to_precision(SpinEfficiencyF*100,4))
    print('dx after decision point (ft)----------', to_precision((BallStateF[0] - xD)/12,4))
    print('dy after decision point (ft)----------', to_precision((BallStateF[1] - yD)/12,4))
    print('dz after decision point (ft)----------', to_precision((BallStateF[2] - zD)/12,4))

    print('\n\nTotals:')
    print('flight time (t)-----------------------', to_precision(t,4))
    print('Vertical break (in)-------------------', to_precision(vBreak*12,4))
    print('Horizontal break (in)-----------------', to_precision(hBreak*12,4))
    print('Number of Revolutions-----------------', to_precision(totalRotations*t,4))

    fileBART = open("y2traj.txt","w+")
    fileBART.write("     y         Traj\n")
    fileBART.write("===================\n")
#    for i in range(len(xP)):
#        traj = 180/np.pi*np.arctan2(uP[i],np.sqrt(vP[i]**2 + wP[i]**2))
##        print(yP[i],traj)
#        fileBART.write("{:<10.3f}{:<10.3f}\n"
#                 .format(yP[i], traj))
    positions = [xP,yP,zP,x0,y0,z0,xD,yD,zD,xF,yF,zF]

    return positions

###############################################################################

def TiltToTime(Tilt):

    TiltTime = (((Tilt)%360)/360)*12
    Hrs = int(TiltTime)
    if Hrs == 0:
        Hrs = 12
    mins = int(TiltTime*60)%60
    return(Hrs,mins)

###############################################################################

def TimeToTilt(Hrs, mins):
    """
    Take the tilt in hrs and mins and turns it into radians
    """
    radHrs = ((Hrs-3)*np.pi/6)
    radmins = (mins*np.pi/360)
    return(radHrs + radmins)

###############################################################################

def anglesTOCart(Vtot, Theta, Psi, SpinRate, Tilt, Gyro, Yangle, Zangle):
    """
    This function is designed merely to generate the balls initial conditions
    It will take various options and output x0,y0,z0,u0,v0,w0,Spinx0,\
    Spiny0,Spinz0,Yangle,Zangle angle 1 and angle 2 are for seam effects
    """

    Theta = Theta*np.pi/180
    Psi = Psi*np.pi/180

    uvmag = Vtot*np.cos(Theta)
    w0 = Vtot*np.sin(Theta)

    u0 = -uvmag*np.sin(Psi)
    v0 = uvmag*np.cos(Psi)

    Tilt = (Tilt) # rad tilt
    Gyro = (Gyro) # rad gyro

    #this is where the changes need to occur to fix the problems with the gyro

    Spinx0 = SpinRate*np.sin(Gyro)*np.sin(Tilt)
    Spiny0 = SpinRate*np.cos(Gyro)
    Spinz0 = -SpinRate*np.sin(Gyro)*np.cos(Tilt)
#    Yangle = 0
#    Zangle = 0

#    print('\nu:',u0,'\nv:',v0,'\nw:',w0)
#    print('\nSpinx0:',Spinx0,'\nSpiny0:',Spiny0,'\nSpinz0:',Spinz0)

    FullState = [u0,v0,w0,Spinx0,Spiny0,Spinz0,Yangle,Zangle]
    return FullState

###############################################################################

def findSSWseams(VelVec,seamPoints,SpinVec, t):
    """
    VelVec is the an array of the velocities u, v, w
    seamPoints is an array of the seam locations with location (0,0,0) being
    the ball's center.

    This function calculates which of the seams are inside a volume whose
    corners are called nodes. The nodes are determined by the normalPlane
    funtion

    Function outputs:
    innerPoints which are the points inside the volume,
    outerPoints which are the points outside the volume,(mostly for plotting)
    and the node points (mostly for plotting)
    """
    nodes = normalPlane(VelVec,SpinVec,t)
    seamPoints = np.transpose(seamPoints)
    innerPointsIndex = (inside_test(seamPoints,nodes))

    innerPoints = []
    outerPoints = []
    for i in range(len(seamPoints)):
        if i in innerPointsIndex:
            outerPoints.append(seamPoints[i])
        else:
            active = activeTest(VelVec, seamPoints, i)
            if active == True:
                innerPoints.append(seamPoints[i])
            else:
                outerPoints.append(seamPoints[i])
    innerPoints = np.asarray(innerPoints)
    outerPoints = np.asarray(outerPoints)

#    inlineTest(VelVec,innerPoints)
    #just to check
#    plotting.plotPointsTest(innerPoints,outerPoints,nodes)
    return innerPoints,outerPoints,nodes

###############################################################################

def inside_test(points , cube3d):
    """
    cube3d  =  numpy array of the shape (8,3) with coordinates in the clockwise order. first the bottom plane is considered then the top one.
    points = array of points with shape (N, 3).

    Returns the indices of the points array which are outside the cube3d
    """
    b1 = cube3d[0]
    b2 = cube3d[1]
#    b3 = cube3d[2]
    b4 = cube3d[3]
    t1 = cube3d[4]
#    t2 = cube3d[5]
    t3 = cube3d[6]
#    t4 = cube3d[7]

    dir1 = (t1-b1)
    size1 = np.linalg.norm(dir1)
    dir1 = dir1 / size1

    dir2 = (b2-b1)
    size2 = np.linalg.norm(dir2)
    dir2 = dir2 / size2

    dir3 = (b4-b1)
    size3 = np.linalg.norm(dir3)
    dir3 = dir3 / size3

    cube3d_center = (b1 + t3)/2.0

    dir_vec = points - cube3d_center

    res1 = np.where( (np.absolute(np.dot(dir_vec, dir1)) * 2) > size1 )[0]
    res2 = np.where( (np.absolute(np.dot(dir_vec, dir2)) * 2) > size2 )[0]
    res3 = np.where( (np.absolute(np.dot(dir_vec, dir3)) * 2) > size3 )[0]

    return list( set().union(res1, res2, res3) )

###############################################################################

def activeTest(VelVec,seamPoints, i):
    """
    since seams in the activation region cannot cause a separated flow to
    become separated again this function will eliminate any inline seams
    """

    AACUrrent = seamPoints[i]
    unit_vector_V = VelVec / np.linalg.norm(VelVec)

    if i  > 0:
        seamLineVecD = seamPoints[i] - seamPoints[i-1]
        unit_vector_Sd = seamLineVecD / np.linalg.norm(seamLineVecD)
        dot_productd = np.dot(unit_vector_Sd, unit_vector_V)
        angled = np.arccos(dot_productd)*180/np.pi
    else:
        seamLineVecD = seamPoints[i] - seamPoints[107]
        unit_vector_Sd = seamLineVecD / np.linalg.norm(seamLineVecD)
        dot_productd = np.dot(unit_vector_Sd, unit_vector_V)
        angled = np.arccos(dot_productd)*180/np.pi

    if i < (107):
        seamLineVecU = seamPoints[i+1] - seamPoints[i]
        unit_vector_Su = seamLineVecD / np.linalg.norm(seamLineVecU)
        dot_productu = np.dot(unit_vector_Su, unit_vector_V)
        angleu = np.arccos(dot_productu)*180/np.pi
    else:
        seamLineVecU = seamPoints[0] - seamPoints[i]
        unit_vector_Su = seamLineVecD / np.linalg.norm(seamLineVecU)
        dot_productu = np.dot(unit_vector_Su, unit_vector_V)
        angleu = np.arccos(dot_productu)*180/np.pi

    if angled > 90:
        angled = angled - 180
    if angleu > 90:
        angleu = angleu - 180

    if abs(angleu) < 35 or abs(angleu) < 35: # or (np.linalg.norm(seamLineVecNU)) > 0.15
        return False
    else:
        return True

###############################################################################

def initializeSeam():
    """
    This function defines the seams of a baseball. It is
    based, in large extant, on the work from
    http://www.darenscotwilson.com/spec/bbseam/bbseam.html
    """
    n = 108 #number of points were calculating on the seam line
    alpha = np.linspace(0,np.pi*2,n)
    x = np.zeros(len(alpha))
    y = np.zeros(len(alpha))
    z = np.zeros(len(alpha))
    R = (2 + 15/16.)/2
    for i in range(len(alpha)):

        x[i] = ((1/13)*R*((9*np.cos(alpha[i]) - 4*np.cos(3*alpha[i]))))
        y[i] = ((1/13)*R*((9*np.sin(alpha[i]) + 4*np.sin(3*alpha[i]))))
        z[i] = ((12/13)*R*np.cos(2*alpha[i]))

    return x,y,z

###############################################################################

def rotateSeam(x, y, z, Spinx,Spiny,Spinz,dt):
    """
    takes an initial seam orientation calculates new seam positions based on
    a cartesian spin rate vector. A rotation vecotr is calculated based on the
    spin rate vector and the time step.
    """
    xn = np.zeros(len(x))
    yn = np.zeros(len(y))
    zn = np.zeros(len(z))
    RotVec = [Spinx*dt,Spiny*dt,Spinz*dt]
    r = R.from_rotvec(RotVec)
    for i in range(len(x)):
        vec = [x[i],y[i],z[i]]
        vecN = r.apply(vec)
        xn[i] = vecN[0]
        yn[i] = vecN[1]
        zn[i] = vecN[2]

    return(xn,yn,zn)

###############################################################################

def findRotVec(sx0,sy0,sz0, sx1,sy1,sz1):
    SpinVecMag0 = np.sqrt(sx0**2 + sy0**2 + sz0**2)
    if SpinVecMag0 == 0:
        return(0,0,0)
    s = (3,3)
    RM = np.zeros(s)

    nx0 = sx0/(SpinVecMag0)
    ny0 = sy0/(SpinVecMag0)
    nz0 = sz0/(SpinVecMag0)

    nvec0 = [nx0,ny0,nz0]


    SpinVecMag1 = np.sqrt(sx1**2 + sy1**2 + sz1**2)

    nx1 = sx1/(SpinVecMag1)
    ny1 = sy1/(SpinVecMag1)
    nz1 = sz1/(SpinVecMag1)

    nvec1 = [nx1,ny1,nz1]

    axis = np.cross(nvec0, nvec1)
    axisLength = np.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    if axisLength != 0:
        axis = axis/axisLength

    x = axis[0]
    y = axis[1]
    z = axis[2]

    angle = np.arccos(np.dot(nvec0,nvec1))

    ca = np.cos(angle)
    sa = np.sin(angle)


    RM[0,0] = 1.0 + (1.0 - ca)*(x**2 - 1.0)
    RM[0,1] = -z*sa + (1.0 - ca)*x*y
    RM[0,2] = y*sa + (1.0 - ca)*x*z
    RM[1,0] = z*sa+(1.0 - ca)*x*y
    RM[1,1] = 1.0 + (1.0 - ca)*(y**2 - 1.0)
    RM[1,2] = -x*sa+(1.0 - ca)*y*z
    RM[2,0] = -y*sa+(1.0 - ca)*x*z
    RM[2,1] = x*sa+(1.0 - ca)*y*z
    RM[2,2] = 1.0 + (1.0 - ca)*(z**2 - 1.0)

    r = R.from_dcm(RM)
    V = R.as_rotvec(r)
#    i(V)
    return(V)

###############################################################################

def ClKensrud(S):
    """ S is the spin factor calulated above, Not in UMBA1.0 Some changes
    need ot be made before this will work
    """
    return (1.1968*np.log(abs(S) + 4.7096))

###############################################################################

def ClCross(S):

    return (1/(2.42 + (0.4/S)))

###############################################################################

def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p
    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

###############################################################################

def RK4(t0,y0,dt,BallConsts, activeSeams):

    n = len(y0)

    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)

    ym = np.zeros(n)
    ye = np.zeros(n)
    y = np.zeros(n)
    slope = np.zeros(n)

    k1 = derivs(t0,y0, BallConsts, activeSeams)
    ym = y0 + (k1*dt*0.5)

    k2 = derivs(t0+dt*0.5, ym, BallConsts, activeSeams)
    ym = y0 + k2*dt*0.5

    k3 = derivs(t0+dt*0.5,ym, BallConsts, activeSeams)
    ye = y0 + k3*dt


    k4 = derivs(t0+dt, ye, BallConsts, activeSeams)

    slope = (k1 + 2*(k2+k3) + k4)/6.0

    y = y0 + slope*dt

    return y
