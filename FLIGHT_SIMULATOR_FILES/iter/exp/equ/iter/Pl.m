%working
%
%
echo off

%fid=fopen('tab_q.dat','r')
%[ntab,count]=fscanf(fid,'%d',1)
%   for i=1:ntab,
%[qq,count]=fscanf(fid,'%g',2);
%pst(i)=qq(1);
%qt(i)=qq(2);
%   end
%fid=fclose(fid)

%fid=fopen('fpol.wr','r')
%[fpt,count]=fscanf(fid,'%g',ntab);
%fid=fclose(fid)
%
%  for i=1:ntab,
%fpt(i)=10.*fpt(i)/(4.*pi);
%  end


fid=fopen('tabppf.wr','r')
[ntab,count]=fscanf(fid,'%d',1)
   for i=1:ntab,
[arr,count]=fscanf(fid,'%g',3);
psitab(i)=arr(1);
pptab(i)=arr(2);
fptab(i)=arr(3);
   end
fid=fclose(fid)


%
%fid=fopen('ddp_corsica.wr','r')
fid=fopen('ddp.wr','r')
%
%
[ipl,count]=fscanf(fid,'%d',1)
%
[q,count]=fscanf(fid,'%g',ipl);

[f,count]=fscanf(fid,'%g',ipl);

[df,count]=fscanf(fid,'%g',ipl);

[psf,count]=fscanf(fid,'%g',ipl);

[sq,count]=fscanf(fid,'%g',ipl);

[dp,count]=fscanf(fid,'%g',ipl);
[cf_av,count]=fscanf(fid,'%g',ipl);
[b2_av,count]=fscanf(fid,'%g',ipl);
%
fid=fclose(fid)
%
figure(1)
%clf
   for i=1:ipl-1,

      ac(i)=(i+0.5-1)/(ipl-1);
      qc(i)=0.5*q(i)/pi;
      fc(i)=f(i);
      sqm(i)=sq(i);
      dfm(i)=df(i);
      dpm(i)=dp(i);
      am(i)=i-1;
      psc(i)=1.-0.5*(psf(i)+psf(i+1));
      psm(i)=1.-psf(i);
   end
   psm(ipl)=1.;
   psc(ipl)=1.;
   dfm(ipl)=df(ipl);
   dpm(ipl)=dp(ipl);
   fc(ipl)=f(ipl);
   qc(ipl)=0.5*q(ipl)/pi;
   ac(ipl)=1.;
%
%
clear a
   for i=1:ipl,
  a(i)=(i-1)/(ipl-1);
   end
%






%title('1D equilibrium plots')

   subplot(2,2,1)
%plot(pst,qt,'k.')

   hold on

plot(psc,qc,'m')
ylabel('q')
xlabel('\psi')

hold on

subplot(2,2,2)
hold on

plot(psc,f,'m')

ylabel('f')
xlabel('\psi')

subplot(2,2,3)
hold on

plot(psm,dfm,'m')
%plot(psitab,fptab,'g.')

ylabel('F*dF/d\psi')
xlabel('\psi')

subplot(2,2,4)

plot(psm,dpm,'m')
%plot(psitab,pptab,'g.')

hold on

ylabel('dp/d\psi')
xlabel('\psi')

hold on


 figure(2)


 subplot(2,2,1)
plot(ac,qc,'m')
title('safety factor')
ylabel('q')
xlabel('sqrt(\phi)')

hold on

subplot(2,2,3)
hold on

title('averaged toroidal current dens.')
plot(a,cf_av,'m')

ylabel('(J*B)_a_v')
xlabel('sqrt(\phi)')

subplot(2,2,2)
hold on

title('averaged magnetic field')
plot(ac,sqrt(b2_av),'m')

ylabel('B')
xlabel('sqrt(\phi)')

subplot(2,2,4)
hold on

title('poloidal current')
plot(ac,f,'m')

hold on

ylabel('F')
xlabel('sqrt(\phi)')






