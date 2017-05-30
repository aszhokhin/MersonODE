#!/usr/bin/env python
# -*- coding: utf-8 -*-

from subprocess import call
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

NDim=5
NPoints=8000
DatFl="D5EQ.DAT"
DTU="d5eq"
XORI=[.8,.9,.8,.9,.8]
TBE=0.0
TEN=20.0
TEMP=0.0

def d5eq(t,v,d):
 X1,X2,X3,X4,X5 = v
 F1=-2*X1+4*X2*X3+4*X4*X5
 F2=-9*X2+3*X1*X3
 F3=-5*X3-7*X1*X2+28.65
 F4=-5*X4-X1*X5
 F5=-X5-3*X1*X4
 d[0]=F1
 d[1]=F2
 d[2]=F3
 d[3]=F4
 d[4]=F5
 return d

def resultates(n,m,t1,t,u):
 t2=0.0
 fil=open(DatFl,'wb')
 for i in range(0,m):
  fil.write('%15.5f' % (t1+(t-t1)*i/(m-1)))
  for j in range(0,n):
   fil.write("     ")
   fil.write('%15.5f' % u[j][i])
  fil.write(" \n")
 fil.close()
 return 0

def soleq(n,m,x,u,z,v,dtu):
 ok=1
 a=1e-5
 hn=.01
 w=1e-5
 k=0
 l=1
 j=0
 t= 0.
 merson(t,z,x,n,a,hn,w,j,ok,dtu)
 for i in range(0,n):
  u[i][0]=x[i]
 s=(v-z)/(m-1)
 t=z
 for l in range(1,m):
  t1=t+s
  merson(t,t1,x,n,a,hn,w,j,ok,dtu)
  t = t1
  for i in range(0,n):
   u[i][l]=x[i]
 return 0

def merson(t,q,y,n,a,h,o, j,l,dtu):
 d0=np.zeros(n)
 d1=np.zeros(n)
 d2=np.zeros(n)
 d3=np.zeros(n)
 d4=np.zeros(n)
 r= 1e-13
 s=0.01
 l= 1
 z=0.

 for k in range(0,n):
  d4[k]=y[k]
 z= t
 s= h
 p=s
 iss=0
 ok=1
 while (ok==1):
  p=s
  c= q-z
  if abs(s)>abs(c):
   s=c
# end
   if abs(c/p)<r:
    break
   iss=1
# end
  for k in range(0,n):
   d0[k]=d4[k]
# end
  f=s/3.
  z=z+f

  for k in range(0,n):
   d1[k]=f*d3[k]
   d4[k]=d1[k]+d0[k]
# end
  dtu(z,d4,d3)
  for k in range(0,n):
   d1[k]=d1[k]*0.5
   d4[k]=f*0.5*d3[k]+d1[k]+d0[k]
# end
  dtu(z,d4,d3)
  z=z+f*0.5
  for k in range(0,n):
   d2[k]=f*4.5*d3[k]
   d4[k]=d2[k]*0.25 +0.75*d1[k]+d0[k]
# end
  dtu(z,d4,d3)
  z=z+s*0.5
  for k in range(0,n):
   d1[k]=f*2*d3[k]+d1[k]
   d4[k]=d1[k]*3-d2[k]+d0[k]
# end
  dtu(z,d4,d3)
  for k in range(0,n):
   d2[k]=-0.5*f*d3[k]-d2[k]+d1[k]*2.
   d4[k] = d4[k]-d2[k]
   d1[k]=abs(a*0.5*d4[k])
   d2[k]=abs(d2[k])
   if abs(d4[k])<=r:
    continue
#  end
   if d2[k]>d1[k]:
    c=s*0.5
    if (abs(c) >= o):
     for k in range(0,n):
      d4[k]=d0[k]
#   end
    z=z-s
    s=c
    iss=0
#  end
   if j==0:
    l=0
    break
#  end
   s=o
   if p<0.:
    s=-s
#  end
   if iss==1:
    break
#  end
   continue
# end
# end
  if iss==1:
   break
# end
  i=0
  for k in range(0,n):
   if (d2[k]>d1[k]/32.):
    i=1
#  end
# end
  if i==0:
   s=s*2
# end
  continue
# end
 h = p
 t=z
 for k in range(0,n):
  y[k]=d4[k]
# end
 TEMP=z
#end
 return 0

n=NDim
m=NPoints
x0=np.zeros(n)
x=np.zeros(n)
u=np.zeros((n,m))

x0=XORI
t1=TBE
t=TEN

for i in range(0,n):
 x[i]=x0[i]

soleq(n,m,x,u,t1,t,d5eq)
resultates(n,m,t1,t,u)

cmd = "gnuplot d5eq.gp"

ret = call(cmd, shell=True)
if ret !=0:
 print "mnd failed"
 sys.exit(1)





