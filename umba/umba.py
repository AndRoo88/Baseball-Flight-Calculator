import numpy as np
import plotting
import processing


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
        Tiltr = processing.TimeToTilt(TiltH, Tiltm)
        
        positions = (processing.PitchedBallTraj(x,y,z,Vtot,Theta,Psi,SpinRate,Tiltr,Gyro,0,0,i))
        plotting.Plotting(positions)
        
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
    
    plotting.plotSFinal(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,i)


    

main()
