import std.stdio;
import std.math;
import std.array;
import std.algorithm;
import std.random;


 int dtu(double t, ref double[4]  v,ref  double[4] d)
{
  double c,g,x,e, dc,dg,dx,de;

   c=v[0]; g=v[1]; x=v[2]; e=v[3];

    dg = e/(1.+e)*c/(1. +c+g);
    dx= g*x/(1.+g+x);

    dc = 1.-200*dg;
    dg = 400*dg-400.0*dx;
    dx = 230*dx-43.5*x;
    de = 1100*x*c/(1+c)/(1.+g*20.) - 20.*e;

 d[0]=dc; d[1]=dg; d[2]=dx; d[3]=de;

 return 0;
}


 void resultates(int n, int m, double t1, double t,  ref double [8000][4] u)
{ int i,j;
  double t2;
  for(i=0; i< m; i++)
  { 
     write(t1+(t-t1)*i/(m-1)," ");
   for(j=0; j< n; j++) write(" ",u[j][i]);  write(" \n"); }
}

 int merson(double *t,double *q,ref double[4] y,int n,double *a,double *h,double *o, int *j,int *l)
{
    static double r = 1e-13;  int i1; double d1;
    double d[20];
    double d4[4];
    double d3[4];
    static double c,f,p,s,z; 
    static int i,k,n2,n3,n4,n31,n41,kn,iis,kn2,kn3,kn4;

    *l= 1; n4= n *4; n3=n*3; n2= n *2; n41=n4+1; n31=n3+1; i1= n;
    for(k=0;k<n;++k) d[k+n4]=y[k];  z= *t;  s= *h; iis= 0;


    while(1)
    { p= s;  c= *q-z;
    if (fabs(s)>fabs(c)){s=c;if(fabs(c/p)<r) break;iis=1;}

    for(k=0;k<n;++k) d[k]=d[k+n4]; f=s/3; 
    d4[0..n]=d[n4..n4+n];
    d3[0..n]=d[n3..n3+n];
    dtu(z,d4,d3);
    d[n3..n3+n]= d3[0..n];

    z=z+f;
    for (k = 0; k<n; ++k) {kn=k+n; kn3=k+n3; kn4=k+n4;
         d[kn]=f*d[kn3]; d[kn4]=d[kn]+d[k];}
    d4[0..n]=d[n4..n4+n];
    d3[0..n]=d[n3..n3+n];
    dtu(z,d4,d3);
    d[n3..n3+n]= d3[0..n];


    for (k = 0; k < n; ++k) {kn = k+n; kn3=k+n3; kn4=k+n4;
         d[kn] *= .5;d[kn4]=f*.5*d[kn3]+d[kn]+d[k]; }
    d4[0..n]=d[n4..n4+n];
    d3[0..n]=d[n3..n3+n];
    dtu(z,d4,d3);
    d[n3..n3+n]= d3[0..n];

    z =z+f*.5;  i1 = n;
    for(k= 0;k<n;++k) {kn = k+n; kn2 = k+n2; kn3 = k+n3; kn4 = k+n4;
        d[kn2]=f*4.5*d[kn3];d[kn4]=d[kn2]*.25 +d[kn]*.75+d[k];}
    d4[0..n]=d[n4..n4+n];
    d3[0..n]=d[n3..n3+n];
    dtu(z,d4,d3);
    d[n3..n3+n]= d3[0..n];

    z =z+s*.5;
    for(k=0;k<n;++k) {kn=k+n; kn2=k+n2; kn3=k+n3;kn4=k+n4; 
        d[kn]=f*2*d[kn3]+d[kn];d[kn4]=d[kn]*3-d[kn2]+d[k];}
    d4[0..n]=d[n4..n4+n];
    d3[0..n]=d[n3..n3+n];
    dtu(z,d4,d3);
    d[n3..n3+n]= d3[0..n];

    for (k = 0; k<n; ++k)
     { kn=k+n; kn2 =k+n2; kn3 =k+n3; kn4=k+n4;
         d[kn2]=-.5*f*d[kn3]-d[kn2]+d[kn]*2;
       d[kn4] = d[k+n4]-d[kn2];  d[kn] = (d1=*a*.5*d[kn4],fabs(d1)); d[kn2]=(d1=d[kn2],fabs(d1));

       if((d1 = d[kn4],fabs(d1))<=r) continue;

       if(d[kn2]>d[kn]) {c=s*.5;
                      if (fabs(c) >= *o) {for(k=0;k<n;++k)
                        d[k+n4]=d[k]; z=z-s; s=c; iis=0;}
                      if(*j==0){*l=0; break;}
                      s=*o;
                      if(p<0.) s=-s;
                      if(iis==1) break;
                      continue;}
    }

    if (iis==1) break;
    i=0;
    for (k=0;k<n;++k) if(d[k+n2]>d[k+n]/32) i=1;
    if(i==0) s=s*2;
    continue;
    }

    *h = p; *t=z; i1=n;
    for (k=0;k<n;++k) y[k]=d[k+n4];
    return 0;
}

  void soleq(int n, int m, ref double x[4], ref double u[4][8000], double b, double v)
  {
    static double a, hn;
    static int i, j, k, l;
    static double s, t, w, t1;
    static int ok;
    double* z;

    z=&b;
    ok=1;
    a=1e-5; hn=.01; w=1e-5; k=0; l=1; j=0; t= 0.;

    merson(&t,z, x,n,&a,&hn,&w,&j,&ok);
    for(i=0;i<n;++i) u[i][0]=x[i];
    s=(v-*z)/(m-1); t=*z;
    for(l=1;l<m;++l){t1=t+s;

    merson(&t, &t1, x, n, &a, &hn, &w, &j, &ok);

     t = t1;
     for(i=0;i< n;++i) u[i][l]=x[i];
    }
}


 void  main(){
  int n= 4;
  int m= 8000;
  int i,j;
  double t,t1;
  double[8000][4] u;
  double[4] x,x0;

  x0=[0.8,0.9,0.8,0.9];

  t1=50.;
  t=80.;
  for(i=0;i<n;++i)  x[i]=x0[i];

  for(i=0;i<n;++i)
  {  for(j=0;j<m;j++)  u[i][j]=0.; }

  soleq(n,m,x,u,t1, t);
  resultates(n,m, t1, t, u);
}


