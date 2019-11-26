from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def trans(tx,ty,tz):
    Tx = np.matrix('{}; {}; {}'.format(tx,ty,tz))
    print Tx
    return Tx;

def rotx(x):
    degx = np.radians(x)
    c,s = np.cos(degx),np.sin(degx)
    Rx = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(1, 0, 0, 0, c, -s, 0, s, c))
    print "The rotation matrix about x axis is: "; print Rx
    return Rx;

def roty(y):
    degy = np.radians(y)
    c,s = np.cos(degy),np.sin(degy)
    Ry = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(c, 0, s, 0, 1, 0, -s, 0, c))
    print "The rotation matrix about y axis is: "; print Ry
    return Ry;

def rotz(z):
    degz = np.radians(z)
    c,s = np.cos(degz),np.sin(degz)
    Rz = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(c, -s, 0, s, c, 0, 0, 0, 1))
    print "The rotation matrix about z axis is: "; print Rz
    return Rz;

def vect(mat):
    VecStart_x = [mat[0,0],mat[0,1],mat[0,2],mat[0,3],mat[0,0],mat[0,4],mat[0,5],mat[0,6],mat[0,7],mat[0,1],mat[0,2],mat[0,3]]
    VecStart_y = [mat[1,0],mat[1,1],mat[1,2],mat[1,3],mat[1,0],mat[1,4],mat[1,5],mat[1,6],mat[1,7],mat[1,1],mat[1,2],mat[1,3]]
    VecStart_z = [mat[2,0],mat[2,1],mat[2,2],mat[2,3],mat[2,0],mat[2,4],mat[2,5],mat[2,6],mat[2,7],mat[2,1],mat[2,2],mat[2,3]]
    VecEnd_x = [mat[0,1],mat[0,2],mat[0,3],mat[0,0],mat[0,4],mat[0,5],mat[0,6],mat[0,7],mat[0,4],mat[0,5],mat[0,6],mat[0,7]]
    VecEnd_y = [mat[1,1],mat[1,2],mat[1,3],mat[1,0],mat[1,4],mat[1,5],mat[1,6],mat[1,7],mat[1,4],mat[1,5],mat[1,6],mat[1,7]]
    VecEnd_z  =[mat[2,1],mat[2,2],mat[2,3],mat[2,0],mat[2,4],mat[2,5],mat[2,6],mat[2,7],mat[2,4],mat[2,5],mat[2,6],mat[2,7]]
    return VecStart_x,VecStart_y,VecStart_z,VecEnd_x,VecEnd_y,VecEnd_z;

tx = input("Enter the amount of translation in x axis: ")
ty = input("Enter the amount of translation in y axis: ")
tz = input("Enter the amount of translation in z axis: ")
print "The translation matrix is: " ;tr = trans(tx,ty,tz)

b = raw_input("Enter whether the rotation is about x,y or z axis: ")

if(b == 'x' or b == 'X'):
    a = input("Enter the degree of x rotation: ")
    rx = rotx(a)
    r = np.hstack((rx,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
   
elif(b == 'y' or b == 'Y'):
    a = input("Enter the degree of y rotation: ")
    ry = roty(a)
    r = np.hstack((ry,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
    
elif(b == 'z' or b == 'Z' ):
    a = input("Enter the degree of z rotation: ")
    rz = rotz(a)
    r = np.hstack((rz,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
if(b == 'xy' or b == 'XY'):
    a = input("Enter the degree of x rotation: ")
    xy = input("Enter the degree of y rotation: ")
    rx = rotx(a)
    ry = roty(xy)
    rxy = np.matmul(rx,ry)
    r = np.hstack((rxy,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
   
elif(b == 'yz' or b == 'YZ'):
    a = input("Enter the degree of y rotation: ")
    yz = input("Enter the degree of z rotation: ")
    ry = roty(a)
    rz = rotz(yz)
    ryz = np.matmul(ry,rz)
    r = np.hstack((ryz,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
    
elif(b == 'zx' or b == 'ZX' ):
    a = input("Enter the degree of x rotation: ")
    zx = input("Enter the degree of z rotation: ")
    rz = rotz(zx)
    rx = rotx(a)
    rzx = np.matmul(rz,rx)
    r = np.hstack((rzx,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res

elif(b == 'xyz' or b == 'XYZ' ):
    a = input("Enter the degree of x rotation: ")
    xyz = input("Enter the degree og y rotation: ")
    zx = input("Enter the degree of z rotation: ")
    rz = rotz(zx)
    ry = roty(xyz)
    rx = rotx(a)
    rxy = np.matmul(rx,ry)
    rxyz = np.matmul(rxy,rz)
    r = np.hstack((rxyz,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
    
else:
    print("Enter a valid data")
    
xs=([0,10,10,0,0,10,10,0])
ys=([0,0,0,0,10,10,10,10])
zs=([0,0,10,10,0,0,10,10])
rs=([1,1,1,1,1,1,1,1])
mat = np.matrix((xs,ys,zs,rs))

rs = np.matmul(res,mat)
print "The multiplied Matrix is: "
print rs

for c, m in [('r', 'o')]:
    ax.scatter(xs, ys, zs, c=c, marker=m)
    
pa=np.asarray(rs[0,])
pb=np.asarray(rs[1,])
pc=np.asarray(rs[2,])   

for c, m in [('b', '^')]:
    ax.scatter(pa, pb, pc, c=c, marker=m)
    
vsx,vsy,vsz,vex,vey,vez=vect(mat)
for i in range(12):
    ax.plot([vsx[i], vex[i]], [vsy[i],vey[i]],[vsz[i],vez[i]],color='r')

asx,asy,asz,aex,aey,aez=vect(rs)
for i in range(12):
    ax.plot([asx[i], aex[i]], [asy[i],aey[i]],[asz[i],aez[i]],color='b')
        
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
plt.close()
