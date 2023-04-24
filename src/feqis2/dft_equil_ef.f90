! fast sine transform from Zhongde Wang, Signal Processing 19 (1990) 91-102 Elsevier


	subroutine fst1(f,g,d,l,m,n)
	
	dimension f(n-1),g(n-1),d(n-1),l(n-1)
	m1=m-1
	n1=n/2
	n2=n1-1
	do 10 i=1,n2
	g(i)=f(i)+f(n-i)
10	f(n-i)=f(i)-f(n-i)
	g(n1)=f(n1)	
	call fst3(g,d,m1,n1)
	n0=0
20	n0=n1+n0
n1=n1/2
m1=m1-1
n2=n1-1
n3=n0+n1
n4=n3+n3
do 30 i=n0+1,n3-1
g(i)=f(i)+f(n4-i)
30	f(n4-i)=f(n4-i)-f(i)
g(n3)=f(n3)
	call fst3(g(n0+1),d,m1,n1)
	if (n1.ge.2) goto 20
	g(n-1)=f(n-1)
	do 40 i=1,n-1
40	f(l(i))=g(i)
return
end


subroutine fst2(f,d,m,n)
dimension f(n),d(1)
n1=n/2
n2=2
j3=1
do 10 i=1,n1
j3=-j3
i2=i+i
i1=i2-1
t=f(i1)
f(i1)=(t+f(i2))*d(n1+i-j3)
10	f(i2)=t-f(i2)
do 50 i=1,m-1
n1=n1/2
n7=n2
n2=n2+n2
n6=0
if (n1.eq.1) j3=0
do 40 j=1,n1
j3=-j3
j2=n1+j
b=d(j2)+d(j2)
n3=n6+1
n6=n6+n2
n8=n6-n7
n5=n8-1
n4=n8+n8
t=f(n8)
f(n8)=(t+f(n6))*d(j2-j3)
f(n6)=t-f(n6)
do 20 k=n3,n5
k1=n7+k
t=f(k)
f(k)=t+f(k1)
20 f(k1)=b*(t-f(k1))
do 30 k=n3,n5
30 f(n4-k)=f(n4-k)+f(k)
40 continue
50 continue
f(n)=d(1)*f(n)
return
end



subroutine fst3(f,d,m,n)
dimension f(n),d(1)
j1=1
n1=n
j3=-1
f(n)=d(1)*f(n)
do 40 i=1,m-1
n2=n1
n1=n1/2
n6=0
do 30 j=1,j1
j2=j1+j
j3=-j3
a=d(j2-j3)
b=d(j2)+d(j2)
n3=n6+1
n6=n2+n6
n5=n6-n1
n4=n5+n5
do 10 k=n3,n5-1
k1=n4-k
f(k)=f(k)+f(k1)
10 f(k1)=f(k1)*b
f(n5)=f(n5)*a
do 20 k=n3,n5
t=f(k+n1)
f(k+n1)=f(k)-t
20 f(k)=f(k)+t
30 continue
40 j1=j1+j1
50 do 60 i=1,j1
i1=i+i
i2=i1-1
j3=-j3
t=d(j1+i-j3)*f(i2)
f(i2)=t+f(i1)
60 f(i1)=t-f(i1)
return
end

subroutine fst4(f,d,m,n)
dimension f(n),d(1)
j1=1
n1=n
do 40 i =1,m
j2=j1
n2=n1
n1=n1/2
n7=0
do 30 j=1,j1
a=d(j2+j)
a=a+a
n6=n7
n7=n7+n2
n3=n6+1
n4=n3+n7
n5=n6+n1
do 10 k=n3,n5
k1=n4-k
f(k)=f(k)+f(k1)
10 f(k1)=a*f(k1)
do 20 k=n3,n5
t=f(k+n1)
f(k+n1)=f(k)-t
20 f(k)=f(k)+t
30 continue
40 j1=j1+j1
do 50 i=1,j1
i1=i+i
i2=i1-1
f(i1)=d(j1+j2)*f(i1)
50 f(i2)=D(j1+i1)*f(i2)
return
end


subroutine coefs(d,n1)
dimension d(n1),c(15)
double precision c,t
c(1)=sqrt(0.5)
d(1)=c(1)
m=2
m0=1
10 do 20 i=m0,m-1
i1=i+i
i2=i1+1
c(i1)=sqrt(0.5*(1.+c(i)))
c(i2)=sqrt(0.5*(1.-c(i)))
d(i1)=c(i1)
20 d(i2)=c(i2)
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
if (m.ge.16) goto 30
goto 10
30 do 40 i=m0,m-1
i1=i+i
i2=i1+1
t=sqrt(0.5*(1.+c(i)))
d(i+i)=t
t=sqrt(0.5*(1.-c(i)))
40 d(i+i+1)=t
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
50 do 60 i=m0,m-1
d(i+i)=sqrt(.5*(1.+d(i)))
60 d(i+i+1)=sqrt(.5*(1.-d(i)))
m0=m0+m0
m=m0+m0
if (m.ge.n1) goto 70
goto 50
70 do 80 i=n1,2,-1
80 d(i)=d(i-1)
return
end


subroutine hdmod(f,g,n)
dimension f(n),g(1)
if (n.le.2) return
n1=n/2
n0=n1
n2=1
10 i1=n
i2=-n2
do 20 i=1,n1
i1=i1-n2
i2=i2+n2
do 20 j=1,n2
20 g(i2+j)=f(i1+j)
n3=0
n4=n2+n2
i1=n
do 30 i=1,n1
n3=n3+n2
i1=i1-n4
i2=n0-n3
do 30 j=1,n2
30 f(i1+j)=f(i2+j)
i1=-n2
i2=-n2
do 40 i=1,n1
i1=i1+n4
i2=i2+n2
do 40 j=1,n2
40 f(i1+j)=g(i2+j)
n2=n2+n2
n1=n1/2
if (n1.ge.2) goto 10
return
end


subroutine ivhdm(f,g,n)
dimension f(n),g(1)
if (n.le.2) return
n1=n/4
i1=1
10 n2=n1+n1
n3=n2+n2
i2=0
do 40 i=1,i1
i3=i2
i2=i2+n3
i4=i3+n2
do 20 j=1,n2
20 g(j)=f(i3+j+j)
do 30 j=2,n2
30 f(i3+j)=f(i3+j+j-1)
do 40 j=1,n1
j1=j+j
j2=j1-1
f(i4+j)=G(j2)
40 f(i4+j2)=g(j1)
n1=n1/2
i1=i1+i1
if (n1.ge.1) goto 10
return
end

subroutine reord(l,n)
dimension l(n)
l(1)=1
l(2)=3
n2=2
10 n1=n2
n2=n1+n1
do 20 i=n1,1,-1
i1=i+i
l(i1)=n2-l(i)
20 l(i1-1)=l(i)
if (n1.lt.n/2) goto 10
n0=0
n2=n/2
n3=n2
30 n2=n2/2
do 40 i=1,n2
i1=n0+i+i-1
40 l(n3+i)=l(i1)+l(i1)
if (n2.eq.1) return
n0=n3
n3=n3+n2
goto 30
return
end










