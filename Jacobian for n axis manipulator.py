import numpy as np
import math
from numpy.linalg import inv

Tx=[]; R=[]; D=[];k=[]; Rx=[];Ts=[];k1=np.matrix(([0],[0],[1]));
ij=[];j=[];Jt=[];Di=[]

def transformation(x):
    a1=x[0,0]
    alp1=x[0,1]
    d1=x[0,2]
    thet1=x[0,3]
    alp,thet=np.radians(alp1),np.radians(thet1)
    ct,st=round(float(np.cos(thet)),5),round(float(np.sin(thet)),5)
    ca,sa=round(float(np.cos(alp)),5),round(float(np.sin(alp)),5)
    sc=round(float(np.sin(thet)*np.cos(alp)),5)
    ss=round(float(np.sin(thet)*np.sin(alp)),5)
    cc=round(float(np.cos(thet)*np.cos(alp)),5)
    cs=round(float(np.cos(thet)*np.sin(alp)),5)
    T = np.matrix((([ct,-sc,ss,a1*ct],[st,cc,-cs,a1*st],[0,sa,ca,d1],[0,0,0,1])))
    return T;

def extractrot(t):
    r11=round(t[0,0],4);r12=round(t[0,1],4);r13=round(t[0,2],4);
    r21=round(t[1,0],4);r22=round(t[1,1],4);r23=round(t[1,2],4);
    r31=round(t[2,0],4);r32=round(t[2,1],4);r33=round(t[2,2],4);
    R=np.matrix(([r11,r12,r13],[r21,r22,r23],[r31,r32,r33]))
    return R;

def extracttra(t):
    dx11=round(t[0,3],5);dx12=round(t[1,3],5);dx13=round(t[2,3],5);
    DX=np.array([dx11,dx12,dx13])
    return DX;

alpha = []; theta = []; d = []; a = []
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
#print "The DH table is:"
#print table
for i in range(int(num)):
    if(i==0):
        Trx = transformation(table[i,])
    else:
        Trx = np.matmul(Trx,transformation(table[i,]))
       
for i in range(int(num)):
    print "Transformation matrix for link",i+1
    Tx.append(transformation(table[i,]))
    print Tx[i]
    if(i==0):
        Ts.append(Trx)
    else:
        Ts.append(np.matmul(inv(Tx[i-1]),Ts[i-1]))
    print "Combine transformations",i+1
    print Ts[i]
    R.append(extractrot(Tx[i]))
    print "Rotation matrix for link",i+1
    print R[i]
    D.append(extracttra(Ts[i]))
    print "Translation matrix for link",i+1
    print D[i]
    if(i==0):
        Rx.append(R[i])
    else:
        Rx.append(np.matmul(Rx[i-1],R[i]))
    print"RX Matrix",i+1
    print Rx[i]
print k1
for i in range(int(num)):
    k.append(np.matmul(Rx[i],k1))
    if(i==0):
        Di.append(D[i])
    else:
        Di.append(np.matmul(Rx[i-1],np.transpose(D[i])))
    print "ud"
    print Di[i]
    print "Updated translation matix"
    print Di[i]
    print "The corresponding k matrix is"
    print k[i]
    if (i==0):
        ij.append(np.cross(np.transpose(k1),D[i]))
    else:
        ij.append(np.cross(np.transpose(k[i]),Di[i]))
    print "Initial elements of column",i+1
    print ij[i]
    
for i in range(int(num)):
    if(i==0):
        j.append(np.hstack((ij[i],np.transpose(k1))))
    else:
        j.append(np.hstack((ij[i],np.transpose(k[i]))))
    #print j[i]
jt=np.asarray(j)
print "Jacobian Matrix is"
J=np.matrix(np.transpose(jt))
print J
Jt = J.transpose()
I = np.identity(6)
JI =inv(np.matmul(J,Jt)+(1600*I))
IJ = np.matmul(Jt,JI)
#print IJ

