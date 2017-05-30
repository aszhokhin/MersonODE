NDim=4;
NPoints=8000;
DatFil ="D4EQ.DAT";
TBE=20.;
TEN=50.;

XORI=[0.8,0.9,0.8,0.9];

TEMP=0

function f = dtu(t,r)
    c=r(1);g=r(2);x=r(3);e=r(4);
    dg_dt = e/(1.0+e)*c/(1.0 +c+g);
    dx_dt = g*x/(1.0+g+x);
    dc_dt = 1.0 - 200.0*dg_dt;
    dg_dt = 400.0*dg_dt-400.0*dx_dt;
    dx_dt = 230.0*dx_dt - 43.565*x;
    de_dt = 1100.0*x*c/(1+c)/(1.0+g*20.0) - 20.0*e;
    f(1)=dc_dt;
    f(2)=dg_dt;
    f(3)=dx_dt;
    f(4)=de_dt;
endfunction



function result(n,m,t1,t2,u)
fl = fopen ("D4EQ.DAT", "w");
 for i=1:m
  fprintf (fl,"%15.5f %15.5f %15.5f %15.5f %15.5f \n",t1+(t2-t1)*i/(m-1),u(1,i),u(2,i),u(3,i),u(4,i));
 end
fclose(fl)
endfunction

function u = soleq(n,m,x,u,b,v)
 t1=b; t=0.0; ok=true;
 i=j=l=0; a=0.00001; hn=0.01;
 w=0.00001; l=1;    j=0;
 t= 0.0; TEMP=t;
 x=merson(t,t1,x,n,a,hn,w,j,ok);
 TEMP=t=t1;

 for i=1:n #1
  u(i,1)=x(i);
 end #1
 s=(v-b)/(m-1.); t=b;

for l=2:m #2
 t1=t+s; TEMP=t;
 x= merson(t, t1, x, n, a, hn, w, j, ok);
 TEMP=t=t1;
 for i=1:n #3
  u(i,l)=x(i);
 end #3
end#2

endfunction # soleq

function y=merson(t,q,y,n,a,h,o, j,l)
 d0=zeros (1,n, "double"); d1=zeros (1,n, "double"); d2=zeros (1,n, "double");
 d3=zeros (1,n, "double"); d4=zeros (1,n, "double");
 r= 1e-13; s=0.01; l= 1; z=0.;

 for k=1:n
  d4(k)=y(k);
 end
  z= t;  s= h;  p=s;  is= 0; ok=true;

while ok==true
 p=s;c= q-z;
    if abs(s)>abs(c)
     s=c;
          if abs(c/p)<r
            break
          end
         is=1;
        end

    for k=1:n
     d0(k)=d4(k);
    end
    f=s/3.;
#    print(d0)
    z=z+f;

    for k=1:n
     d1(k)=f*d3(k);
     d4(k)=d1(k)+d0(k);
    end
    d3=dtu(z,d4,d3);

    for k=1:n
     d1(k)=d1(k)*0.5;
     d4(k)=f*0.5*d3(k)+d1(k)+d0(k);
    end
    d3=dtu(z,d4,d3);

    z=z+f*0.5;
    for k=1:n
     d2(k)=f*4.5*d3(k);
     d4(k)=d2(k)*0.25 +0.75*d1(k)+d0(k);
    end
    d3=dtu(z,d4,d3);

    z=z+s*0.5;
    for k=1:n
     d1(k)=f*2*d3(k)+d1(k);
     d4(k)=d1(k)*3-d2(k)+d0(k);
    end
    d3=dtu(z,d4,d3);

    for k=1:n
     d2(k)=-0.5*f*d3(k)-d2(k)+d1(k)*2.;
     d4(k) = d4(k)-d2(k);
     d1(k)=abs(a*0.5*d4(k));
     d2(k)=abs(d2(k));
     if abs(d4(k))<=r
      continue
     end

     if d2(k)>d1(k)
      c=s*0.5;
      if (abs(c) >= o)
       for k=1:n
         d4(k)=d0(k);
       end
       z=z-s;
       s=c;
       is=0;
      end

      if j==0
       l=0;
       break
      end
      s=o;
      if p<0.
       s=-s;
      end
      if is==1
       break
      end
       continue
      end
     end

    if is==1
     break
    end
    i=0;

    for k=1:n
     if (d2(k)>d1(k)/32.)
      i=1;
     end
    end

    if i==0
      s=s*2;
    end
    continue
   end

    h = p;  t=z;
    for k=1:n
     y(k)=d4(k);
    end
 TEMP=z;

endfunction

n=NDim;
m=NPoints;
t1=TBE;
t=TEN;
x=XORI;
u=zeros(n,m)*0.1;
u=soleq(n,m,x,u,t1,t);
result(n,m,t1,t,u);




