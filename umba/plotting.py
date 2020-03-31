import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits import mplot3d

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
    
###############################################################################    
    
def plotSeams(xSeam,ySeam,zSeam,sx,sy,sz,num):
    """
    plots in 3 dimensions 3 arrays
    """
    radius = (2. + 15/16)/2 #in
    
    SpinVecMag = np.sqrt(sx**2 + sy**2 + sz**2)
    
    nx = sx/(SpinVecMag)
    ny = sy/(SpinVecMag)
    nz = sz/(SpinVecMag)
    
    nvec = [nx,ny,nz] #normalized Vector
    SpinAxis = [i*radius for i in nvec]
    SpinAxisneg = [-i*radius for i in nvec]
    
    xax = [SpinAxis[0],SpinAxisneg[0]]
    yax = [SpinAxis[1],SpinAxisneg[1]]
    zax = [SpinAxis[2],SpinAxisneg[2]]
    
    
    x = xSeam
    y = ySeam
    z = zSeam

    fig = plt.figure(num)
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(x, y, z, color = 'r')
    ax.scatter(x[0],y[0],z[0], color = 'b', s = 80, marker = 's')
    ax.plot(xax,yax,zax)
#    ax.scatter(SpinAxisneg[0], SpinAxisneg[1], SpinAxisneg[2], s = 60, marker = 'v')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
#    t = num*.0021
    plt.title('time: %1.3f' %num)

    max_range = 1.3554770667211864
    mid_x = -0.00044881302015253866
    mid_y = 0.0
    mid_z = 0.00029216404804432994
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    plt.show()
    
###############################################################################    
        
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
    plt.xlim(0,62.5)
    plt.ylim(0,89.)
    for i in range(j):
        plt.plot(pY[i],pZ[i], label=i)
        plt.legend()
        plt.scatter(IY[i],IZ[i]*12, s=100, c = 'g')
        plt.scatter(DY[i],DZ[i]*12, s=100, c = 'y')
        plt.scatter(FY[i],FZ[i]*12, s=100, c = 'r')
        plt.savefig("Side.jpg")
    plt.show()
