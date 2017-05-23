TEMP=Float64(0.0)

function d3eq(t::Float64,r::Array{Float64,1},f::Array{Float64,1})
    (x, y, z) = r
    dx_dt=-0.2*x+(z-1)*y
    dy_dt=-1.08*y+(z+1)*x
    dz_dt=-0.1*z-2*x*y+5
     f1 = [dx_dt, dy_dt, dz_dt]
     f[1:3]=f1[1:3]
end

function result(n::Int64, m::Int64, t1::Float64, t2::Float64, u::Array{Float64})

f = open("D3EQ.DAT","w")
t=(t2-t1)

 for i=1:m
  t=t1+(i-1)*(t2-t1)/(m-1.)
  println(f,@sprintf("%15.5f %15.5f %15.5f %15.5f",t,u[1,i],u[2,i],u[3,i]))
 end
close(f)
end

function soleq(n::Int,m::Int,x::Array{Float64,1},u::Array{Float64,2},b::Float64,v::Float64,dtu)
    t1=b
    t=0.0
    ok=true
    i=j=l=0
    a=0.00001
    hn=0.01
    w=0.00001
    l=1
    j=0
    t= 0.0
    TEMP=t
    merson(t,t1,x,n,a,hn,w,j,ok,dtu)
    TEMP=t=t1
   for i=1:n #1
     u[i,1]=x[i]
   end #1

    s=(v-b)/(m-1.)
    t=b

  for l=2:m #2
     t1=t+s
     TEMP=t
     merson(t, t1, x, n, a, hn, w, j, ok,dtu)
     TEMP=t=t1

    for i=1:n #3
      u[i,l]=x[i]
    end #3
  end#2
  u
end # soleq


function merson(t,q,y,n,a,h,o, j,l,dtu)

 d0=zeros(Float64,n)
 d1=zeros(Float64,n)
 d2=zeros(Float64,n)
 d3=zeros(Float64,n)
 d4=zeros(Float64,n)

 r= 1e-13
 s=0.01
 l= 1
 z=0.

 for k=1:n
  d4[k]=y[k]
 end
  z= t
  s= h
  p=s
  is= 0

    ok=true

    while ok==true
     p=s
     c= q-z
    if abs(s)>abs(c)
     s=c
          if abs(c/p)<r
            break
          end
         is=1
        end

    for k=1:n
     d0[k]=d4[k]
    end
    f=s/3.
#    print(d0)
    z=z+f


    for k=1:n
     d1[k]=f*d3[k]
     d4[k]=d1[k]+d0[k]
    end
    dtu(z,d4,d3)

    for k=1:n
     d1[k]=d1[k]*0.5
     d4[k]=f*0.5*d3[k]+d1[k]+d0[k]
    end
    dtu(z,d4,d3)

    z=z+f*0.5
    for k=1:n
     d2[k]=f*4.5*d3[k]
     d4[k]=d2[k]*0.25 +0.75*d1[k]+d0[k]
    end
    dtu(z,d4,d3)

    z=z+s*0.5
    for k=1:n
     d1[k]=f*2*d3[k]+d1[k]
     d4[k]=d1[k]*3-d2[k]+d0[k]
    end
    dtu(z,d4,d3)

    for k=1:n
     d2[k]=-0.5*f*d3[k]-d2[k]+d1[k]*2.
     d4[k] = d4[k]-d2[k]
     d1[k]=abs(a*0.5*d4[k])
     d2[k]=abs(d2[k])
     if abs(d4[k])<=r
      continue
     end

     if d2[k]>d1[k]
      c=s*0.5
      if (abs(c) >= o)
       for k=1:n
         d4[k]=d0[k]
       end
       z=z-s
       s=c
       is=0
      end

      if j==0
       l=0
       break
      end
      s=o
      if p<0.
       s=-s
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
    i=0

    for k=1:n
     if (d2[k]>d1[k]/32.)
      i=1
     end
    end

    if i==0
      s=s*2
    end
    continue
   end

    h = p
    t=z
    for k=1:n
     y[k]=d4[k]
    end
 TEMP=z
end


#function rlt()
 t1=10.
 t2=50.
 n=3
 m=8000
 x=[0.8,0.9,0.8]
 u=zeros(n,m)
 soleq(n,m,x,u,t1,t2,d3eq)
 result(n,m,t1,t2,u)
#end
