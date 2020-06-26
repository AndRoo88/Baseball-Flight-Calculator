import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

def AccelPlots(aX,aY,aZ,TF,j,variable):

    plt.figure(7)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in x')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aX[i])))
        plt.plot(t, aX[i], label = variable[i])
        plt.legend()
        plt.savefig('AX.jpg')
    plt.show()

    plt.figure(8)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in y')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aY[i])))
        plt.plot(t, aY[i], label = variable[i])
        plt.legend()
        plt.savefig('AY.jpg')
    plt.show()

    plt.figure(9)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in z')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aZ[i])))
#        t = np.linspace(0,TF,num = len(aZ[i]))
        plt.plot(t, aZ[i], label = variable[i])
        plt.legend()
        plt.savefig('AZ.jpg')
    plt.show()

###############################################################################

def AccelPlots2(aX,aY,aZ,TF,j):

    plt.figure(7)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in x')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aX[i])))
        plt.plot(t[i], aX[i], label = i)
        plt.legend()
        plt.savefig('AX.jpg')
    plt.show()

    plt.figure(8)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in y')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aY[i])))
        plt.plot(t[i], aY[i], label = i)
        plt.legend()
        plt.savefig('AY.jpg')
    plt.show()

    plt.figure(9)
    plt.xlabel('time (sec)')
    plt.ylabel('acceleration (ft/s^2)')
    plt.title('acceleration in z')
    for i in range(j):
        t=(np.linspace(0,TF[i],num = len(aZ[i])))
        plt.plot(t[i], aZ[i], label = i)
        plt.legend()
        plt.savefig('AZ.jpg')
    plt.show()

###############################################################################


def plotSFinalIter(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,j,variable):
    """

    works with UMBANQ.py and plots all the pitches that were tested.
    the p variables (pX, pY, pZ) are all the position arrays. the I variables
    indicate the initial positions, the D variables indicate the ball's
    position at the decision point, and the F variables show the final point of
    the ball as it crosses the plate. j is the iteration and variable is the
    variable that changes if multiple pitches are tested.

    """

    plt.figure(4,figsize=(3,10))
    plt.xlabel('x (in)')
    plt.ylabel('y (ft)')
    plt.title('Bird\'s Eye View')
#    plt.ylim(0,max(pY) + 2.5)
    for i in range(j):
        plt.plot(pX[i],pY[i], label = variable[i])
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
    plt.axis('equal')
    for i in range(j):
        plt.plot(pX[i],pZ[i], label=variable[i])
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
        plt.plot(pY[i],pZ[i], label=variable[i])
        plt.legend()
        plt.scatter(IY[i],IZ[i]*12, s=100, c = 'g')
        plt.scatter(DY[i],DZ[i]*12, s=100, c = 'y')
        plt.scatter(FY[i],FZ[i]*12, s=100, c = 'r')
        plt.savefig("Side.jpg")
    plt.show()

###############################################################################

def plotSFinal(pX,pY,pZ,IX,IY,IZ,DX,DY,DZ,FX,FY,FZ,j):
    """

    works with umba.py and plots all the pitches that were tested.
    the p variables (pX, pY, pZ) are all the position arrays. the I variables
    indicate the initial positions, the D variables indicate the ball's
    position at the decision point, and the F variables show the final point of
    the ball as it crosses the plate

    """

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
    plt.axis('equal')
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

###############################################################################

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

def plotSeams(innerPoints,outerPoints,sx,sy,sz,num, VelVec, nodes):
    """

    plots in 3 dimensions 3 arrays and includes the activates seams, spin axis,
    velocity vector.

    """
#    ims = []
    VelVecN = np.linalg.norm(VelVec)
    VelVec = VelVec/VelVecN
    if len(innerPoints) > 0:
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


        xo = outerPoints[:,0]
        yo = outerPoints[:,1]
        zo = outerPoints[:,2]

        xi = innerPoints[:,0]
        yi = innerPoints[:,1]
        zi = innerPoints[:,2]

        xn = nodes[:,0]
        yn = nodes[:,1]
        zn = nodes[:,2]

        fig = plt.figure(num)
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(xo, yo, zo, color = 'r')
        ax.scatter(xi, yi, zi, color = 'g', s = 90)
        ax.scatter(0,0,0, color = 'k', s = 80)
        ax.scatter(xn,yn,zn, color = 'b', s = 70)
    #    ax.scatter(x[0],y[0],z[0], color = 'b', s = 80, marker = 's')
        ax.plot(xax,yax,zax)
        ax.quiver(0,0,0,VelVec[0],VelVec[1],VelVec[2],color = 'k')
    #    ax.scatter(SpinAxisneg[0], SpinAxisneg[1], SpinAxisneg[2], s = 60, marker = 'v')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    #    t = num*.0021
        plt.title('time: %1.3f' %num)


        max_range = 1.3554770667211864
        mid_x = -0.00
        mid_y = 0.0
        mid_z = 0.00
        adjustment = 0.7
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range*adjustment, mid_z + max_range*adjustment)
        plt.show()

#        plt.figure(3)
#        plt.scatter(xo,zo)
#        plt.scatter(xi,zi)
#        plt.show()

    else:
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


        xo = outerPoints[:,0]
        yo = outerPoints[:,1]
        zo = outerPoints[:,2]

#        xi = innerPoints[:,0]
#        yi = innerPoints[:,1]
#        zi = innerPoints[:,2]

        xn = nodes[:,0]
        yn = nodes[:,1]
        zn = nodes[:,2]

        fig = plt.figure(num)
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(xo, yo, zo, color = 'r')
#        ax.scatter(xi, yi, zi, color = 'g')
    #    ax.scatter(x[0],y[0],z[0], color = 'b', s = 80, marker = 's')
        ax.scatter(0,0,0, color = 'k', s = 80)
        ax.scatter(xn,yn,zn, color = 'b', s = 70)
        ax.plot(xax,yax,zax)
        ax.quiver(0,0,0,VelVec[0],VelVec[1],VelVec[2],color = 'k')
    #    ax.scatter(SpinAxisneg[0], SpinAxisneg[1], SpinAxisneg[2], s = 60, marker = 'v')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    #    t = num*.0021
        plt.title('time: %1.3f' %num)

#this next section of code makes the axis all be to the same scale

        max_range = 1.3554770667211864
        mid_x = -0.00
        mid_y = 0.0
        mid_z = 0.00
        adjustment = 0.7
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range*adjustment, mid_z + max_range*adjustment)
        plt.show()

#        plt.figure(3)
#        plt.scatter(xo,zo)
#    #    plt.scatter(xi,zi)
#        plt.show()


###############################################################################

def plotPointsTest(outerPoints,innerPoints,nodes,num):
    """

    Can be useful for bug finding and fixing

    """
    fig = plt.figure(num)
    ax = fig.add_subplot(111,projection='3d')

    ax.plot(nodes[:,0],nodes[:,1],nodes[:,2],color = 'r')

    ax.scatter(outerPoints[:,0],outerPoints[:,1],outerPoints[:,2], color = 'r')
    ax.scatter(innerPoints[:,0],innerPoints[:,1],innerPoints[:,2], color = 'g')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show
