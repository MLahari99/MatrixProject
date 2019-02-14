
import numpy as np
from sympy import *
import mpmath as mp
import matplotlib.pyplot as plt

def magnitude(A):
    return sqrt(((A.T)*A)[0])

def dir_vec(A,B):
    return B-A

def norm_vec(AB):
    return MatMul(omat,dir_vec(AB))

def plot_line(A,B,n):
    len = 10
    lam_1 = np.linspace(0,n,len)
    AB = np.zeros((2,10))

    for i in range(10):
        temp = A + lam_1[i]*(B-A)
        AB[:,i] = temp.T

    return AB

omat = Matrix([[0,1],[-1,0]],dtype='float')
dvec = Matrix([-1,1],dtype='float')

A = Matrix([2,-4],dtype='float')
C = ((-1/2)*A)
B = np.array([2,0],dtype='float')
r  = sqrt(4+magnitude(C)**2)                                                            #Finding the radius of the circle whose equation is given
r1 = 3                                                                                                                 #r1 - Radius of the given circle

X = Matrix([2,2],dtype='float')                                                          #Point of  contact of two circles
CX  = dir_vec(C,X)

k = symbols('k')
C1 = (k)*(CX)+C                                                                                         #C1 - Centre of the circle to be found

#Equations to solve for k
y = solve(((C1-C).T*(C1-C))[0] - (r1+r)**2,k)
z = solve(((C1-X).T*(C1-X))[0] - (r)**2,k)



p = symbols('p')
C1 = 2*(CX)+C                                                                                              #Computing C1 after finding the value of k
I = Matrix([p,0])                                                                                          #General point on x axis
w = solve(((C1-I).T*(C1-I))[0] - (r)**2,p)                                       #Solving for p given the point also lies on the circle

#Intercept points on X axis
I1 =np.array([w[0],0])
I2 = np.array([w[1],0])

#Intercept value
print('Intercept cut by the circle on the x axis :  ' ,round(sqrt(((I2-I1)*(I2-I1).T)[0]),3))

#To plot points and circles and lines

plt.xlim([-7,14])
plt.ylim([-5,10])

circle1=plt.Circle((C[0],C[1]),r,color='r',fill=False)
plt.gca().add_artist(circle1)
circle2=plt.Circle((C1[0],C1[1]),r1,color='b',fill=False)
plt.gca().add_artist(circle2)

plt.scatter(2,2,c='c',marker='o')
plt.text(2.1,2.0,'(2,2)')

plt.scatter(C[0],C[1],c='g',marker='o')
plt.text(C[0]*(1+0.1),C[1],"C({}, {})".format(round(C[0]),round(C[1])))

plt.scatter(C1[0],C1[1],c='y',marker='o')
plt.text(C1[0]*(1+0.1),C1[1],"C1({}, {})".format(round(C1[0]),round(C1[1])))

plt.scatter(I1[0],I1[1],c='m',marker='o')
plt.text(I1[0]*(1+0.05),-0.87,"X1({}, {})".format(round(I1[0],2),round(I1[1],2)))

plt.scatter(I2[0],I2[1],c='k',marker='o')
plt.text(I2[0]*(1+0.05),-0.6,"X2({}, {})".format(round(I2[0],2),round(I2[1],2)))

I = plot_line(I1,I2,1)
CC1 = plot_line(A,X,2)
plt.plot(CC1[0,:],CC1[1,:])
plt.plot(I[0,:],I[1,:])

plt.grid()
plt.show()













