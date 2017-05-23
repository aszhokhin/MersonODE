#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NPoints 8000
#define NDim 4
#define DatFl "D4EQ.DAT"
#define DTU d4eq
#define XORI {.8,.9,.8,.9}
#define TBE 50
#define TEN 80

 int main(){
  static int c4 = NDim;
  static int c400 = NPoints;
  static int i;
  static double t,u[NDim][NPoints],x[NDim],t1,x0[NDim]=XORI;
  extern  int resultates(int , int , double , double , double (*) [NPoints]);
  extern  void *DTU(double , double *, double *);
  extern  int soleq(int , int , double *, double (*) [NPoints] , double , double ,void *dtu(double , double *, double *));
  t1=TBE; t=TEN;
  for(i=0;i<c4;++i) x[i]=x0[i];
  soleq(c4, c400, x, u, t1, t,DTU);
  resultates(c4, c400, t1, t, u);
//  system("gnuplot d4eq.gp");
  return 0;
}

 int resultates(int n, int m, double t1, double t,  double u[][NPoints])
{ int i,j; FILE *str;
  double t2; str = fopen(DatFl, "wt");
  for(i=0; i< m; i++)
  {
   fprintf(str,"%15.5f ", t1+(t-t1)*i/(m-1));
   for(j=0; j< n; j++) fprintf(str,"%15.5f ", u[j][i]);
   fprintf(str," \n");
  }
     fclose(str);
     return 0;
}

    int soleq(int n, int m, double *x, double u[][NPoints], double b, double v, void *dtu(double , double *, double *))
{
    static double a, hn;
    static int i, j, k, l;
    static double s, t, w, t1;
    static long int ok;
    extern  int merson(double *,double *,double *,int ,double *,double *,double *,int *,long int *, void *dtu(double , double *, double *));
    double* z;
    z=&b;
    ok=1;
    a=1e-5; hn=.01; w=1e-5; k=0; l=1; j=0; t= 0.;
    merson(&t,z,x,n,&a,&hn,&w,&j,&ok,dtu);
    for(i=0;i<n;++i) u[i][0]=x[i];
    s=(v-*z)/(m-1); t=*z;
    for(l=1;l<m;++l){t1=t+s;
     merson(&t, &t1, x, n, &a, &hn, &w, &j, &ok,dtu);
     t = t1;
     for(i=0;i< n;++i) u[i][l]=x[i];
    }
    return 0;
}


 int merson(double *t,double *q,double *y,int n,double *a,double *h,double *o,
                        int *j,long int *l, void *dtu(double , double *, double *) )
{
    static double r = 1e-13;  int i1; double d1;
    static double c,d[500],f,p,s,z; 
    static int i,k,n2,n3,n4,n31,n41,kn,is,kn2,kn3,kn4;
    
    *l= 1; n4= n *4; n3=n*3; n2= n *2; n41=n4+1; n31=n3+1; i1= n;
    for(k=0;k<n;++k) d[k+n4]=y[k];  z= *t;  s= *h; is= 0;
    while(1)
    { p= s;  c= *q-z;
    if (fabs(s)>fabs(c)){s=c;if(fabs(c/p)<r) break;is=1;}

    for(k=0;k<n;++k) d[k]=d[k+n4]; f=s/3; dtu(z,d+n4,d+n3);

    z=z+f;
    for (k = 0; k<n; ++k) {kn=k+n; kn3=k+n3; kn4=k+n4;
         d[kn]=f*d[kn3]; d[kn4]=d[kn]+d[k];}
    dtu(z, d+n4,d+n3);

    for (k = 0; k < n; ++k) {kn = k+n; kn3=k+n3; kn4=k+n4;
         d[kn] *= .5;d[kn4]=f*.5*d[kn3]+d[kn]+d[k]; }
    dtu(z,d+n4,d+n3);

    z =z+f*.5;  i1 = n;
    for(k= 0;k<n;++k) {kn = k+n; kn2 = k+n2; kn3 = k+n3; kn4 = k+n4;
        d[kn2]=f*4.5*d[kn3];d[kn4]=d[kn2]*.25 +d[kn]*.75+d[k];}
    dtu(z,d+n4,d+n3);

    z =z+s*.5;
    for(k=0;k<n;++k) {kn=k+n; kn2=k+n2; kn3=k+n3;kn4=k+n4; 
        d[kn]=f*2*d[kn3]+d[kn];d[kn4]=d[kn]*3-d[kn2]+d[k];}
    dtu(z,d+n4,d+n3);

    for (k = 0; k<n; ++k)
     { kn=k+n; kn2 =k+n2; kn3 =k+n3; kn4=k+n4;
         d[kn2]=-.5*f*d[kn3]-d[kn2]+d[kn]*2;
       d[kn4] = d[k+n4]-d[kn2];  d[kn] = (d1=*a*.5*d[kn4],fabs(d1)); d[kn2]=(d1=d[kn2],fabs(d1));

       if((d1 = d[kn4],fabs(d1))<=r) continue;

       if(d[kn2]>d[kn]) {c=s*.5;
                      if (fabs(c) >= *o) {for(k=0;k<n;++k)
                        d[k+n4]=d[k]; z=z-s; s=c; is=0;}
                      if(*j==0){*l=0; break;}
                      s=*o;
                      if(p<0.) s=-s;
                      if(is==1) break;
                      continue;}
    }
    
    if (is==1) break;
    i=0;
    for (k=0;k<n;++k) if(d[k+n2]>d[k+n]/32) i=1;
    if(i==0) s=s*2;
    continue;
    }

    *h = p; *t=z; i1=n;
    for (k=0;k<n;++k) y[k]=d[k+n4];
    return 0;
}
#define c v[0]
#define g v[1]
#define x v[2]
#define e v[3]
    
#define dc d[0]
#define dg d[1]
#define dx d[2]
#define de d[3]
         
 void *d4eq(double t, double *v, double *d)
{    
    dg = e/(1.+e)*c/(1. +c+g);
    dx= g*x/(1.+g+x);

    dc = 1.-200*dg;
    dg = 400*dg-400.0*dx;
    dx = 230*dx-43.5*x;
    de = 1100*x*c/(1+c)/(1.+g*20.) - 20.*e;
        
 return 0;
}   
    
#define c c
#define g g
#define x x
#define e e

#define dc dc
#define dg dg
#define dx dx
#define de de

