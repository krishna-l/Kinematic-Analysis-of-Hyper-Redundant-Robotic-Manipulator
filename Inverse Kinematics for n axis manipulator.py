import numpy as np
import math
from numpy.linalg import inv

beta = 0.1
tol = 0.1
lamb = 40
iterations =5000
V = 20

Tx=[]; R=[]; D=[];k=[]; Rx=[];Ts=[];k1=np.matrix(([0],[0],[1]));
ij=[];j=[];alpha=[];theta=[];d=[];a=[];Di=[];k2=np.matrix(([0],[0],[0]));w=[];ti=[]


def transformation(x):
    a1=x[0,0]
    alp1=x[0,1]
    d1=x[0,2]
    thet1=x[0,3]
    alp,thet=np.radians(alp1),np.radians(thet1)
    ct,st=float(np.cos(thet)),float(np.sin(thet))
    ca,sa=float(np.cos(alp)),float(np.sin(alp))
    sc=float(np.sin(thet)*np.cos(alp))
    ss=float(np.sin(thet)*np.sin(alp))
    cc=float(np.cos(thet)*np.cos(alp))
    cs=float(np.cos(thet)*np.sin(alp))
    T = np.matrix(([ct,-sc,ss,a1*ct],[st,cc,-cs,a1*st],[0,sa,ca,d1],[0,0,0,1]))
    return T;

def extractrot(t):
    r11=t[0,0];r12=t[0,1];r13=t[0,2];
    r21=t[1,0];r22=t[1,1];r23=t[1,2];
    r31=t[2,0];r32=t[2,1];r33=t[2,2];
    R=np.matrix(([r11,r12,r13],[r21,r22,r23],[r31,r32,r33]))
    return R;

def extracttra(t):
    dx11=t[0,3];dx12=t[1,3];dx13=t[2,3];
    DX=np.array([dx11,dx12,dx13])
    return DX;

def extractcor(c):
    x=c[0];y=c[1];z=c[2]
    return x,y,z;

def fk(table,num):
    for i in range(int(num)):
        if(i==0):
            Trx = transformation(table[i,])
        else:
            Trx = np.matmul(Trx,transformation(table[i,]))
    c=extracttra(Trx)
    x,y,z=extractcor(c)
    return x,y,z

def jacobian(x,num):
    Trx=[]; Tx=[]; Ts=[]; Rx=[]; R=[]; D=[]; Di=[]; jt=[]
    J=[]; Jt=[]; JI=[]; IJ=[]; I=[]; J=[]; j=[]; ij=[]
    table=x
    for i in range(int(num)):
        if(i==0):
            Trx = transformation(table[i,])
        else:
            Trx = np.matmul(Trx,transformation(table[i,]))
    for i in range(int(num)):
        Tx.append(transformation(table[i,]))
        if(i==0):
            Ts.append(Trx)
        else:
            Ts.append(np.matmul(inv(Tx[i-1]),Ts[i-1]))
        R.append(extractrot(Tx[i]))
        D.append(extracttra(Ts[i]))
        if(i==0):
            Rx.append(R[i])
        else:
            Rx.append(np.matmul(Rx[i-1],R[i]))
        k.append(np.matmul(Rx[i],k1))
        if(i==0):
            Di.append(D[i])
        else:
            Di.append(np.matmul(Rx[i-1],np.transpose(D[i])))
        if (i==0):
            ij.append(np.cross(np.transpose(k1),D[i]))
            j.append(np.hstack((ij[i],np.transpose(k2))))
        else:
            ij.append(np.cross(np.transpose(k[i]),Di[i]))
            j.append(np.hstack((ij[i],np.transpose(k2))))         
    jt=np.asarray(j)
    J=np.matrix(np.transpose(jt))
    Jt=[]
    Jt = J.transpose()
    I = np.identity(6)
    JI =inv(np.matmul(J,Jt)+(lamb*lamb*I))
    IJ = np.matmul(Jt,JI)
    return IJ

def error(x,y,z,dx,dy,dz):
    errorx = float(dx - x); errory = float(dy - y); errorz = float(dz - z)
    magnitude = float(math.sqrt(math.pow(errorx, 2) + math.pow(errory, 2) + math.pow(errorz, 2)))
    Vecx = V*(errorx/magnitude) ; Vecy = V*(errory/magnitude); Vecz = V*(errorz/magnitude)
    return errorx,errory,errorz,magnitude,Vecx,Vecy,Vecz;

num = input("Enter number of links: ")
for i in range(int(num)):
    print "Enter a for link",i+1
    aa = input()
    a.append(aa)
    print "Enter alpha for link",i+1
    al = input()
    alpha.append(al)
    print "Enter d for link",i+1
    dis = input()
    d.append(dis)
    print "Enter theta for link",i+1
    t = input()
    theta.append(float(t))
ta = np.matrix((a, alpha, d, theta))
table=ta.transpose()
print "The DH table is:"
print table
dx = input("Enter the X-axis destination end effector coordinates: ")
dy = input("Enter the Y-axis destination end effector coordinates: ")
dz = input("Enter the Z-axis destination end effector coordinates: ")
wx = 0#input("Enter the wx value")
wy = 0#input("Enter the wy value")
wz = 0#input("Enter the wz value")

for i in range(iterations):
    x,y,z=fk(table,num)
    print x,y,z
    errx,erry,errz,mag,Vecx,Vecy,Vecz = error(x,y,z,dx,dy,dz)
    Vec = np.matrix(([Vecx],[Vecy],[Vecz],[wx],[wy],[wz]))
    if( mag > tol):
        t = mag/V;
        delt = float(beta *t)
        T=t+delt
        ti.append(T)
        J1 = jacobian(table,num)
        deld = np.matmul(J1,Vec)
        for i in range(int(num)):
            table[i,3] = table[i,3] + delt*(np.degrees(deld[i,0]))
            #print table[i,3]
            
        
    else:
        print " You have reached your destination"
        print "No of iterations \n",i
        break
else:
    print "You have reached the end of your iterations"
