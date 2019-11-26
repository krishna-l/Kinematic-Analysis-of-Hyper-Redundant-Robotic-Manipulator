import numpy as np
from numpy.linalg import inv
import math
import matplotlib.pyplot as plt
beta = 0.1
tol = 0.1
lamb = 40
iterations = 100000
V = 20
thi1= []
thi2=[]
ti = []
T = 0;
Alpha = []
thdot = 0;
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
    JI =inv(np.matmul(J,Jt) + lamb*lamb*I)
    IJ = np.matmul(Jt,JI)
    return IJ;

def error(x,y,ex,ey):
    errorx = float(ex - x); errory = float(ey - y)
    magnitude = float(math.sqrt(math.pow(errorx, 2) + math.pow(errory, 2)))
    Vecx = V*(errorx/magnitude) ; Vecy = V*(errory/magnitude)
    return errorx , errory, magnitude,Vecx,Vecy;   

d1 = float(input("Enter the angle of the first joint: "))
d2 = float(input("Enter the angle of the second joint: "))
l1 = input("Enter the first limb length: ")
l2 = input("Enter the second limb length: ")

ex = float(input("Enter the destination x coordinate of the end effector: "))
ey = float(input("Enter the destination y coordinate of the end effector: "))

for i in range(iterations):
    x, y = (fk(l1,l2,d1,d2))
    print x,y
    errx,erry,mag,Vecx,Vecy = error(x,y,ex,ey)
    if( mag > tol):
        
        t = mag/V;
        delt = float(beta *t)
        Vec = np.matrix([[Vecx],[Vecy]])
        IJ = jacobian(l1,l2,d1,d2)
        deld = np.matmul(IJ,Vec)
        
        w1 = deld[1,0]
        w2 = deld[0,0]
        d1 = d1 + delt*(np.degrees(deld[0,0]))
        d2 = d2 + delt*(np.degrees(deld[1,0]))
        T = T+ delt
 
        
        thi1.append(w1)
        thi2.append(w2)
        ti.append(T)
       
        
    else:
        print " You have reached your destination"
        print "The Number of iterations: ",i 
        break
    
else:
    print "You have reached the end of your iterations"
plt.subplot(211)   
plt.plot(ti,thi1)
plt.xlabel('time')
plt.ylabel('w1')
#plt.plot(ti,Alpha)
plt.subplot(212)
plt.plot(ti,thi2)
plt.xlabel('time')
plt.ylabel('w2')
#plt.show()
