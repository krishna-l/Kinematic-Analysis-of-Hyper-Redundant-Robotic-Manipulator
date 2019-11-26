import numpy as np
from numpy.linalg import inv
import math
import matplotlib.pyplot as plt
import time
beta = 0.01
tol = 0.1
lamb = 40
iterations = 10000
mg = []
ii = []

def fk(l1,l2,d1,d2):
    dr1,dr2 = np.radians(d1), np.radians(d2)
    x = float(((l1*np.cos(dr1))+(l2*np.cos(dr1+dr2))))
    y = float(((l1*np.sin(dr1))+(l2*np.sin(dr1+dr2))))
    return x,y;

def jacobian(l1,l2,d1,d2): 
    t1,t2=np.radians(d1),np.radians(d2)
    x1 = float(-l1*np.sin(t1)-l2*np.sin(t1+t2)); x2 = float(-l2*np.sin(t1+t2))
    y1 = float(l1*np.cos(t1)+l2*np.cos(t1+t2)); y2 = float(l2*np.cos(t1+t2))
    J = np.matrix([[x1,x2],[y1,y2]])
    Jt = J.transpose()
    I = np.identity(2)
    JI =np.matmul(J,Jt); lambI= lamb*lamb*I
    JJI = np.add(JI,lambI)
    JJJI = inv(JJI)
    IJ = np.matmul(Jt,JJJI)
    return IJ;

def error(x,y,ex,ey):
    errorx = float(ex - x); errory = float(ey - y)
    magnitude = float(math.sqrt(math.pow(errorx, 2) + math.pow(errory, 2)))
    return errorx , errory, magnitude;
   
d1 = float(input("Enter the angle of the first joint: "))
d2 = float(input("Enter the angle of the second joint: "))
l1 = input("Enter the first limb length: ")
l2 = input("Enter the second limb length: ")

ex = float(input("Enter the x coordinate of the end effector: "))
ey = float(input("Enter the y coordinate of the end effector: "))

for i in range(iterations):
    x, y = (fk(l1,l2,d1,d2))
    print x,y
    errx,erry,mag = error(x,y,ex,ey)
    if( mag > tol):
        mg.append(mag);
        ii.append(i);
        delx = float(beta * errx)
        dely = float(beta * erry)
        df = np.matrix([[delx],[dely]])
        J = jacobian(l1,l2,d1,d2)
        deld = np.matmul(J,df)
        d1 = d1 + np.degrees(deld[0,0])
        d2 = d2 + np.degrees(deld[1,0])
    else:
        print " You have reached your destination"
        break
else:
    print "You have reached the end of your iterations"

plt.plot(ii,mg)
plt.xlabel('Iterations')
plt.ylabel('Magnitude of error')
plt.show()


