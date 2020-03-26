import numpy as np
import math
import matplotlib.pyplot as plt


def main():
        
    print("This baseball trajectory calculator models \
still air at sea level with about 60% humidity.\nThe ambient\
 conditions are fixed can be adjusted in the code \
if needed or made into variable at a later time.\n\
\nAll initial ball state variables have default values that approximate\
 a 90 mph ball with no spin.\nThe release point is only asked once\
and remains the same for subsequent pitches.\nAll other variables\
can be changed for comparison. To retain the current value, just press return")
    
    
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
    
    #x y and z here are typical of an approximately 6' tall rhp
    # these are the defualt initial values and can be changed here or as the code runs
    # as the code
    x = -1
    y = 5.5
    z = 6
    Vtot = 90
    Theta = 0
    Psi = 0
    SpinRate = 0.001
    TiltH = 0
    Tiltm = 0
    SpinE = 100
            
    print('\n\n\ncurent release distance left to right is', x)
    Qx = (input('distance left from center of rubber (right haded pitchers should have negative numbers) (ft): '))
    if Qx == "":
        x = x
    else:
        x = float(x)
        
    print('\n\n\ncurent release distance from rubber', y)
    Qy = (input('release distance from rubber (should be approximatly the stride lenght)(ft): '))
    if Qy == "":
        y = y
    else:
        y = float(Qy)
        
    print('\n\n\ncurent release height is', z)
    Qz = (input('height of ball from field at release (ft): '))
    if Qz == "":
        z = z
    else:
        z = float(Qz)
        
    i = 0
    repeat = True
    while repeat == True:
            
        if i == 0:
            print('\n\nIf you want to keep the current value simply hit return. Otherwise enter a new value.\n\n')
        print("Current initial speed set to ",Vtot)
        QVtot = (input('what is the ball\'s total initial speed (mph): '))
        if QVtot == "":
            Vtot = Vtot
        else:
            Vtot = float(QVtot)
            
        print("Current vertical release angle set to ",Theta)
        QTheta = (input('what is the ball\'s vertical release angle (deg): '))
        if QTheta == "":
            Theta = Theta
        else:
            Theta = float(QTheta)
            
        print("Current horizontal release angle set to ", Psi)
        QPsi = (input('what is the ball\'s horizontal release angle(deg): '))
        if QPsi == "":
            Psi = Psi
        else:
            Psi = float(QPsi)
            
        print("Current initial Spin Rate set to ", SpinRate)
        QSpinRate = (input('what is the ball\'s initial spin rate (rpm): '))
        if QSpinRate == "":
            SpinRate = SpinRate
        else:
            SpinRate = float(QSpinRate)
            
        print("Current initial tilt hours set to ", TiltH)
        QTiltH = (input('what is the ball\'s initial hours tilt (hrs): '))
        if QTiltH == "":
            TiltH = TiltH
        else:
            TiltH = float(QTiltH)
            
        print("Current initial tilt minutes set to ", Tiltm)
        QTiltm = (input('what is the ball\'s initial minutes tilt (mins): '))
        if QTiltm == "":
            Tiltm = Tiltm
        else:
            Tiltm = float(QTiltm)
            
        print("Current initial spin efficiency set to ", SpinE)
        QSpinE = (input('what is the ball\'s initial spin efficiency (%): '))
        if QSpinE == "":
            SpinE = SpinE
        else:
            SpinE = float(QSpinE)
        if SpinE == 100:
            Gyro = np.arccos(1 - (SpinE/100))
        else:
            TiltHnewUp = TiltH + 3
            if TiltHnewUp > 12:
                TiltHnewUp = int(TiltHnewUp - 12)
            else:
                TiltHnewUp = int(TiltHnewUp)
            TiltHnewDn = TiltH - 3
            if TiltHnewDn < 1:
                TiltHnewDn = int(TiltHnewDn + 12)
            else:
                TiltHnewDn = int(TiltHnewDn)
                
            print('if', TiltHnewUp,':',int(Tiltm),'is forward enter " r "')
            print('if', TiltHnewDn,':',int(Tiltm),'is forward enter " l "')
            
            leftRightGyro = ''
            while leftRightGyro  != 'l' and leftRightGyro != 'r':
                leftRightGyro = input("l/r")
                if leftRightGyro == 'l':
                    Gyro = np.arccos(1 - (SpinE/100))
                elif leftRightGyro == 'r':
                    Gyro = np.pi - np.arccos(1- (SpinE/100))
        print(Gyro)
        Tiltr = TimeToTilt(TiltH, Tiltm)
        
        positions = (PitchedBallTraj(x,y,z,Vtot,Theta,Psi,SpinRate,Tiltr,Gyro,0,0,i))
        Plotting(positions)
        
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
        
        leave = False
        k = 0
        
        while leave == False:
            k = k+1
            if k == 1:
                Again = input("Would you like to look at another pitch?\n")
            elif k > 6:
                Again = 'n'
            elif k > 5:
                Again = input("Last Chance.Would you like to look at another pitch (y/n)?\n")
            else:
                Again = input("Would you like to look at another pitch (y/n)?\n")
             
                
            if Again == 'y' or Again == 'yes' or Again == 'Y' or Again == 'YES' or Again == 'Yes':
                repeat = True
                leave = True
            elif Again == 'n' or Again == 'no' or Again == 'N' or Again == 'NO' or Again == 'No':
                repeat = False
                leave = True
            else:
                leave = False
            
        i = i + 1
    
    plotSFinal(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,i)
    

def anglesTOCart(Vtot, Theta, Psi, SpinRate, Tilt, Gyro, angle1, angle2):
    """
    This function is designed merely to generate the balls initial conditions
    It will take various options and output x0,y0,z0,u0,v0,w0,Spinx0,\
    Spiny0,Spinz0,angle1,angle2 angle 1 and angle 2 are for seam effects
    """
    
    Theta = Theta*np.pi/180
    Psi = Psi*np.pi/180
    
    uvmag = Vtot*np.cos(Theta)
    w0 = Vtot*np.sin(Theta)

    u0 = -uvmag*np.sin(Psi)
    v0 = uvmag*np.cos(Psi)
    
    Tilt = (Tilt) # rad tilt
    Gyro = (Gyro) # rad gyro
    
    Spinx0 = SpinRate*np.sin(Gyro)*np.sin(Tilt)
    Spiny0 = SpinRate*np.cos(Gyro)
    Spinz0 = -SpinRate*np.sin(Gyro)*np.cos(Tilt)
    angle1 = 0
    angle2 = 0
    
    print('\nu:',u0,'\nv:',v0,'\nw:',w0)
    print('\nSpinx0:',Spinx0,'\nSpiny0:',Spiny0,'\nSpinz0:',Spinz0)
    
    
    
    FullState = [u0,v0,w0,Spinx0,Spiny0,Spinz0,angle1,angle2]
    
    
    
    return FullState

def PitchedBallTraj(x,y,z,Vtot, Theta, Psi, SpinRate, Tilt, Gyro, angle1, angle2,i):
    
    """
        Primay inputs are: initial position, x0, y0, and z0 with origin at the
        point of home plate, x to the rright of the catcher, y from the catcher 
        towards the pitcher, and z straight up. Initial velocities
        u0, v0, and w0 which are the speeds of the ball in x, y, and z 
        respectivley. And spin rates
        
    
       GUMBA1.0: This code uses a constant Cd and rod cross's model for CL
       Predictions. Seam Orientation is not accounted for. Air Density is
       considered only at sea level at 60% relative humidity. but can be easily
       altered
    """
    
###############################################################################  

    FullState = anglesTOCart(Vtot, Theta, Psi, SpinRate, Tilt, Gyro, 0,0)
#    print(FullState)
    
    x0 = x
    y0 = 60.5 - y
    z0 = z
    u0 = FullState[0]
    v0 = FullState[1]
    w0 = FullState[2]
    Spinx0 = FullState[3]
    Spiny0 = FullState[4]
    Spinz0 = FullState[5]
    angle1 = FullState[6]
    angle2 = FullState[7]
###############################################################################
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
    
    Spinx0 = Spinx0 * .104719754 #rad/s
    Spiny0 = Spiny0 * -.104719754
    Spinz0 = Spinz0 * .104719754
    
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
    while BallState0[1] > 0. and BallState0[2] > 0. and t < 10:
        #need to input a non-magnus ball path indicator.
        t = t + dt
        BallState1 = RK4(t, BallState0, dt, BallConsts)
        
        fileBT.write("{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}\n"
                 .format(t,BallState1[0],BallState1[1],BallState1[2],BallState1[3],BallState1[4],BallState1[5],BallState1[6],BallState1[7],BallState1[8]))
        
        BallState0 = BallState1
            
        xP.append(BallState1[0]*12)
        yP.append(BallState1[1])
        zP.append(BallState1[2]*12)
        uP.append(BallState1[3])
        vP.append(BallState1[4])
        wP.append(BallState1[5])
        
    DecisionPointStep = int(.2/dt)
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
    
    positions = [xP,yP,zP,x0,y0,z0,xD,yD,zD,xF,yF,zF]
    
    return positions

def Plotting(positions):
    
    """
    xP, yP, and zP are the arrays that contain all the velocity positions. x0,
    xD, and xF are the beginning, decision point, and last x poisitions of the 
    ball. Same for Y and Z
    """
    xP = positions[0]
    yP = positions[1]
    zP = positions[2]
    x0 = positions[3]
    y0 = positions[4]
    z0 = positions[5]
    xD = positions[6]
    yD = positions[7]
    zD = positions[8]
    xF = positions[9]
    yF = positions[10]
    zF = positions[11]
    
    plt.figure(1,figsize=(3,10))
    plt.xlabel('x (in)')
    plt.ylabel('y (ft)')
    plt.title('Bird\'s Eye View')
    plt.ylim(0,max(yP) +2.5)
    plt.plot(xP,yP)
    plt.scatter(x0*12,y0, s=100, c = 'g')
    plt.scatter(xD*12,yD, s=100, c = 'y')
    plt.scatter(xF*12,yF, s=100, c = 'r')
#    hold on
    plt.show()
    
    plt.figure(2, figsize=(3,6))
    plt.xlabel('x (in)')
    plt.ylabel('z (in)')
    plt.title('Catcher\'s Perspective')
    plt.plot(x0,z0,'g')
#    plt.xlim(-17., 17.)
    plt.ylim(0,max(zP) + 4)
    plt.plot(xP,zP)
    plt.scatter(x0*12,z0*12, s=100, c = 'g')
    plt.scatter(xD*12,zD*12, s=100, c = 'y')
    plt.scatter(xF*12,zF*12, s=100, c = 'r')
    plt.show()
    
    plt.figure(3,figsize=(10,3))
    plt.xlabel('y (ft)')
    plt.ylabel('z (in)')
    plt.title('Side View')
    plt.xlim(0,62.5)
    plt.ylim(0,max(zP) + 9)
    plt.plot(yP,zP)
    plt.scatter(y0,z0*12, s=100, c = 'g')
    plt.scatter(yD,zD*12, s=100, c = 'y')
    plt.scatter(yF,zF*12, s=100, c = 'r')
    plt.show()
    
def plotSFinal(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,j):
    
    plt.figure(4,figsize=(3,10))
    plt.xlabel('x (in)')
    plt.ylabel('y (ft)')
    plt.title('Bird\'s Eye View')
#    plt.ylim(0,max(pY) + 2.5)
    for i in range(j):
        plt.plot(pX[i],pY[i], label = i)
        plt.legend()
        plt.scatter(IX[i]*12,IY[i], s=100, c = 'g')
        plt.scatter(DX[i]*12,DY[i], s=100, c = 'y')
        plt.scatter(FX[i]*12,FY[i], s=100, c = 'r')
    plt.show()
    plt.savefig("BirdsEye.jpg")
    
    plt.figure(5, figsize=(3,6))
    plt.xlabel('x (in)')
    plt.ylabel('z (in)')
    plt.title('Catcher\'s Perspective')
#    plt.ylim(0,max(pZ) + 4)
    for i in range(j):
        plt.plot(pX[i],pZ[i], label=i)
        plt.legend()
        plt.scatter(IX[i]*12,IZ[i]*12, s=100, c = 'g')
        plt.scatter(DX[i]*12,DZ[i]*12, s=100, c = 'y')
        plt.scatter(FX[i]*12,FZ[i]*12, s=100, c = 'r')
    plt.show()
    plt.savefig("Catcher.jpg")
    
    plt.figure(6,figsize=(10,3))
    plt.xlabel('y (ft)')
    plt.ylabel('z (in)')
    plt.title('Side View')
    plt.xlim(0,62.5)
#    plt.ylim(0,max(pZ) + 9)
    for i in range(j):
        plt.plot(pY[i],pZ[i], label=i)
        plt.legend()
        plt.scatter(IY[i],IZ[i]*12, s=100, c = 'g')
        plt.scatter(DY[i],DZ[i]*12, s=100, c = 'y')
        plt.scatter(FY[i],FZ[i]*12, s=100, c = 'r')
    plt.show()
    plt.savefig("Side.jpg")
    
def TiltToTime(Tilt):
    
    TiltTime = (((Tilt)%360)/360)*12    
    Hrs = int(TiltTime)
    if Hrs == 0:
        Hrs = 12
    mins = int(TiltTime*60)%60
    return(Hrs,mins)
    
def TimeToTilt(Hrs, mins):
    """
    Take the tilt in hrs and mins and turns it into deg
    """
    radHrs = ((Hrs-3)*np.pi/6)
    radmins = (mins*np.pi/360)
    return(radHrs + radmins)
    
def derivs(t, BallState, BallConsts):
    
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
    
    Cl = 1/(2.32 + (0.4/S))
    CdConst = 0.33
    
    aDragx = -c0*CdConst*VelTot*u
    aDragy = -c0*CdConst*VelTot*v
    aDragz = -c0*CdConst*VelTot*w
    
    aSpinx = c0*(Cl/SpinRate)*VelTot*(Spiny*w - Spinz*v)
    aSpiny = c0*(Cl/SpinRate)*VelTot*(Spinz*u - Spinx*w)
    aSpinz = c0*(Cl/SpinRate)*VelTot*(Spinx*v - Spiny*u)
    
    ax = aDragx + aSpinx
    ay = aDragy + aSpiny
    az = aDragz + aSpinz - 32.2
    
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

def ClKensrud(r,Spin,u):
    """ r is the radius of the ball, Spin is the spin rate and u is the velocity
    It still needs to be worked out. As it is it only works for one dimensional
    movement. Not in GUMBA1.0
    """
    if u == 0. or Spin == 0:
        CL = 0.
    elif Spin < 0: 
        CL = -1.1968*np.log(abs(r*Spin/u)) - 4.7096
    else:
        CL = 1.1968*np.log(abs(r*Spin/u)) + 4.7096
    return CL

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

def RK4(t0,y0,dt,BallConsts):
    
    n = len(y0)
    
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)
    
    ym = np.zeros(n)
    ye = np.zeros(n)
    y = np.zeros(n)
    slope = np.zeros(n)
    
    k1 = derivs(t0,y0, BallConsts)
    ym = y0 + (k1*dt*0.5)
    
    k2 = derivs(t0+dt*0.5, ym, BallConsts)
    ym = y0 + k2*dt*0.5
    
    k3 = derivs(t0+dt*0.5,ym, BallConsts)
    ye = y0 + k3*dt
    
    
    k4 = derivs(t0+dt, ye, BallConsts)
    
    slope = (k1 + 2*(k2+k3) + k4)/6.0
    
    y = y0 + slope*dt
    
    return y


main()
