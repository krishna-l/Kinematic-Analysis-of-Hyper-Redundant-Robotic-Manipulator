import numpy as np
import math

def trans(tx,ty,tz):
    Tx = np.matrix('{}; {}; {}'.format(tx,ty,tz))
    print Tx
    return Tx;

def rotx(x):
    degx = np.radians(x)
    c,s = np.cos(degx),np.sin(degx)
    Rx = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(c, -s, 0, s, c, 0, 0, 0, 1))
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
    Rz = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(1, 0, 0, 0, c, -s, 0, s, c))
    print "The rotation matrix about z axis is: "; print Rz
    return Rz;

tx = input("Enter the amount of translation in x axis: ")
ty = input("Enter the amount of translation in y axis: ")
tz = input("Enter the amount of translation in z axis: ")
print "The translation matrix is: " ;tr = trans(tx,ty,tz)

a = input("Enter the degree of rotation: ")
b = raw_input("Enter whether the rotation is about x,y or z axis: ")

if(b == 'x' or b == 'X'):
    rx = rotx(a)
    r = np.hstack((rx,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
   
elif(b == 'y' or b == 'Y'):
    ry = roty(a)
    r = np.hstack((ry,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
    
elif(b == 'z' or b == 'Z' ):
    rz = rotz(a)
    r = np.hstack((rz,tr))
    res = np.vstack((r, np.array((0,0,0,1))))
    print " Resultant Transformation Matrix is: "; print res
else:
    print("Enter a valid data")
