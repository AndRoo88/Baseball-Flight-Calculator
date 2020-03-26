import matplotlib.pyplot as plt

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
        plt.savefig("BirdsEye.jpg")
    plt.show()
    
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
        plt.savefig("Catcher.jpg")
    plt.show()
    
    plt.figure(6,figsize=(10,3))
    plt.xlabel('y (ft)')
    plt.ylabel('z (in)')
    plt.title('Side View')
#    plt.xlim(0,62.5)
#    plt.ylim(0,max(pZ) + 9)
    for i in range(j):
        plt.plot(pY[i],pZ[i], label=i)
        plt.legend()
        plt.scatter(IY[i],IZ[i]*12, s=100, c = 'g')
        plt.scatter(DY[i],DZ[i]*12, s=100, c = 'y')
        plt.scatter(FY[i],FZ[i]*12, s=100, c = 'r')
        plt.savefig("Side.jpg")
    plt.show()