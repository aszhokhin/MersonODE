package main

import (
        "fmt"
        "os"
        "math"
       )

func check(e error) {
    if e != nil {
        panic(e)
    }
}

  func d4eq(t float64, y []float64, f []float64) int{
   var c,g,x,e,dc,dg,dx,de float64
    c=y[0]; g=y[1]; x=y[2]; e=y[3]
    dg = e/(1.+e)*c/(1. +c+g)
    dx= g*x/(1.+g+x)
    dc = 1.-200*dg
    dg = 400.*dg-400.0*dx
    dx = 230.*dx-43.5*x
    de = 1100.*x*c/(1.+c)/(1.+g*20.)-20.*e
    f[0]=dc; f[1]=dg; f[2]=dx; f[3]=de
    return 0
}

func main() {
var n,m int
var t,t1 float64
var u [4][8000]float64
n=4; m=8000
t1=50.;t=80.
x:=[4]float64{.8, .9, .8, .9}

   soleq(n,m,x[0:],u[0:][0:],t1,t,d4eq)
   resultates(m,m,t1,t, u[0:][0:])
}

func resultates(n int, m int, t1 float64, t float64, u [][8000]float64) {
var i int
var t2,t3 float64
//var dtu fn
//fn=d4eq
f, err := os.Create("D4EQ.DAT")
check(err)
defer f.Close()

t2=t-t1
for i = 0; i < m; i++ { 
    t3=t1+float64(i)*t2/(float64(m)-1.)
    fmt.Fprintf(f,"%15.5f %15.5f %15.5f %15.5f  %15.5f\n", t3, u[0][i],u[1][i],u[2][i],u[3][i])}
f.Close()
}

   func soleq(n int,m int,x[] float64,u[][8000] float64,b float64,v float64,dtu func(float64, []float64, []float64) int) int {

    var a,hn float64
    var i,j,l int
    var s,t,w,t1 float64
    var ok int
    var z *float64
    z=&b
    ok=1
    a=.00001; hn=.01; w=.00001; l=1; j=0; t= 0.
    merson(&t,z,x[0:],n,&a,&hn,&w,&j,&ok,dtu)
    for i=0; i<n; i++ { u[i][0]=x[i]}
//    s=(v-*z)/(float64(m)-1.); t=*z
    s=(v-b)/(float64(m)-1.); t=b
    for l=1;l<m; l++ {t1=t+s
    merson(&t, &t1, x[0:], n, &a, &hn, &w, &j, &ok,dtu)
     t = t1
     for i=0; i< n; i++ { u[i][l]=x[i]}
    }
    return 0
}

 func merson(t *float64,q *float64,y []float64,n int,a *float64,h *float64,o *float64,j *int,l *int, dtu func(float64, []float64, []float64) int) int {
    var i,k,n2,n3,n4,kn,is,kn2,kn3,kn4 int
    var c,f,p,s,z  float64
    var d[500]float64
    r:= 1e-13
    s=.01
    *l= 1; n4= n*4; n3=n*3; n2= n*2
    for k=0;k<n;k++ {d[k+n4]=y[k];  z= *t;  s= *h; is= 0}
    ok := true

    for ok { p=s;  c= *q-z
    if math.Abs(s)>math.Abs(c) {s=c;if math.Abs(c/p)<r { break}; is=1}

    for k=0;k<n;k++ { d[k]=d[k+n4]} 
    f=s/3.
    dtu(z,d[n4:],d[n3:])
    z=z+f

    for k = 0; k<n; k++ {kn=k+n; kn3=k+n3; kn4=k+n4; d[kn]=f*d[kn3]; d[kn4]=d[kn]+d[k]}
    dtu(z, d[n4:],d[n3:])

    for k = 0; k < n; k++ {kn = k+n; kn3=k+n3; kn4=k+n4; d[kn] *= .5;d[kn4]=f*.5*d[kn3]+d[kn]+d[k] }
    dtu(z,d[n4:],d[n3:])

    z =z+f*.5
    for k= 0;k<n;k++ {kn = k+n; kn2 = k+n2; kn3 = k+n3; kn4 = k+n4; d[kn2]=f*4.5*d[kn3];d[kn4]=d[kn2]*.25 +d[kn]*.75+d[k]}
    dtu(z,d[n4:],d[n3:])

    z =z+s*.5
    for k=0;k<n;k++ {kn=k+n; kn2=k+n2; kn3=k+n3;kn4=k+n4; d[kn]=f*2*d[kn3]+d[kn];d[kn4]=d[kn]*3-d[kn2]+d[k]}
    dtu(z,d[n4:],d[n3:])

    for k = 0; k<n; k++ { kn=k+n; kn2 =k+n2; kn3 =k+n3; kn4=k+n4; d[kn2]=-.5*f*d[kn3]-d[kn2]+d[kn]*2.
       d[kn4] = d[k+n4]-d[kn2];  d[kn] =math.Abs(*a*.5*d[kn4]); d[kn2]=math.Abs(d[kn2])
       if math.Abs(d[kn4])<=r {continue}

       if d[kn2]>d[kn] {c=s*.5
                      if  math.Abs(c) >= *o {for k=0;k<n;k++ {d[k+n4]=d[k]}
                      z=z-s; s=c; is=0}
                      if *j==0 {*l=0; break}
                      s=*o
                      if p<0. { s=-s}
                      if is==1 {break}
                      continue}
    }

    if is==1 {break}
    i=0
    for k=0;k<n;k++ {if d[k+n2]>d[k+n]/32 {i=1}}
    if i==0 {s=s*2}
    continue}

    *h = p; *t=z
    for k=0;k<n;k++ {y[k]=d[k+n4]}
    return 0
}
