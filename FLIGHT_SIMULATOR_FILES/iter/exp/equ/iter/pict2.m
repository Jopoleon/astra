% working
clf
clc
clear all
echo off
%
fid=fopen('out.wr','r+');
%
[ni,count]=fscanf(fid,'%d',1);
[nj,count]=fscanf(fid,'%d',1);
[ni1,count]=fscanf(fid,'%d',1);
[nj1,count]=fscanf(fid,'%d',1);
[ni2,count]=fscanf(fid,'%d',1);
[nj2,count]=fscanf(fid,'%d',1);
[nxb,count]=fscanf(fid,'%d',1);
%
[r,count]=fscanf(fid,'%g',ni);
[z,count]=fscanf(fid,'%g',nj);
[psi,count]=fscanf(fid,'%g',[ni,nj]);
psi=psi';
[cur,count]=fscanf(fid,'%g',[ni,nj]);
cur=cur';
[ipr,count]=fscanf(fid,'%d',[ni,nj]);
ipr=ipr';
%
[rm,count]=fscanf(fid,'%g',1);
[zm,count]=fscanf(fid,'%g',1);
[psim,count]=fscanf(fid,'%g',1);
[rx,count]=fscanf(fid,'%g',1);
[zx,count]=fscanf(fid,'%g',1);
[psix,count]=fscanf(fid,'%g',1);
[psib,count]=fscanf(fid,'%g',1);
%
[rxb,count]=fscanf(fid,'%g',nxb);
[zxb,count]=fscanf(fid,'%g',nxb);
%
[psii,count]=fscanf(fid,'%g',[ni,nj]);
psii=psii';
[psie,count]=fscanf(fid,'%g',[ni,nj]);
psie=psie';
[psin,count]=fscanf(fid,'%g',[ni,nj]);
%
fid=fclose(fid);
%
fid=fopen('outp.wr','r+')
%
[nr_pl,count]   = fscanf(fid,'%d',1);
[nt_pl,count]   = fscanf(fid,'%d',1);
[nr1_pl,count]  = fscanf(fid,'%d',1);
[nt1_pl,count]  = fscanf(fid,'%d',1);
[nr2_pl,count]  = fscanf(fid,'%d',1);
[nt2_pl,count]  = fscanf(fid,'%d',1);
[ipl_pl,count]  = fscanf(fid,'%d',1);
[r_pl,count]    = fscanf(fid,'%g',[nr_pl,nt_pl]);
[z_pl,count]    = fscanf(fid,'%g',[nr_pl,nt_pl]);
[g_pl,count]    = fscanf(fid,'%g',[nr_pl,nt_pl]);
[cur_pl,count]  = fscanf(fid,'%g',[nr_pl,nt_pl]);
[psi_pl,count]  = fscanf(fid,'%g',[nr_pl,nt_pl]);
[psii_pl,count]  = fscanf(fid,'%g',[nr_pl,nt_pl]);
[psie_pl,count]  = fscanf(fid,'%g',[nr_pl,nt_pl]);
%[time,count]=fscanf(fid,'%g',1);
%
fid=fclose(fid)
ss=size(psi_pl)
%
fid=fopen('pascon.wr','r+')
%
[nc,count]    = fscanf(fid,'%d',1);
[ncpfc,count] = fscanf(fid,'%d',1);
[nfw,count]   = fscanf(fid,'%d',1);
[nbp,count]   = fscanf(fid,'%d',1);
[nvv,count]   = fscanf(fid,'%d',1);
[rc,count]    = fscanf(fid,'%g',nc);
[zc,count]    = fscanf(fid,'%g',nc);
%
fid=fclose(fid)
%
ss=size(psi)

fid=fopen('limpnt.dat','r+')

[nl,count]   = fscanf(fid,'%d',1);
 for il=1:nl,
[arrd,count]   = fscanf(fid,'%g',2);
rl(il)=arrd(1); 
zl(il)=arrd(2);
end
size(arrd)
fid=fclose(fid);

%
fid=fopen('propoi.wr','r+')
%
[npro,count] = fscanf(fid,'%d',1);
[rpro,count] = fscanf(fid,'%g',npro);
[zpro,count] = fscanf(fid,'%g',npro);
[fpro,count] = fscanf(fid,'%g',npro);
%
fid=fclose(fid);
%
%
%
%
size(psie)
%pause  %mesh, press any key to continue
%
for j=1:nt1_pl,
  rb(j)  = r_pl(ipl_pl,j);
  zb(j)  = z_pl(ipl_pl,j);
  rbb(j)  = r_pl(nr_pl,j);
  zbb(j)  = z_pl(nr_pl,j);
  rv(j)  = r_pl(2,j);
  zv(j)  = z_pl(2,j);
  rv1(j) = r_pl(3,j);
  zv1(j) = z_pl(3,j);
end
%
%for j=1:nc-ncpfc,
%  rves(j)  = rc(ncpfc+j);
%  zves(j)  = zc(ncpfc+j);
%end
%
for i=1:ipl_pl,
for j=1:nt1_pl,
  rpl(i,j)   = r_pl(i,j);
  zpl(i,j)   = z_pl(i,j);
  pspl(i,j)  = 0.;
  curpl(i,j) = cur_pl(i,j);
end
end
%
%
figure(1)
%
clf
%axis('image')
hold on
colormap(jet)
mesh(rpl,zpl,pspl)
plot(rx,zx,'k+')
view(0,90)

h=plot(rb,zb,'m')
set(h,'linewidth',2)
plot(rc,zc,'r+')
%plot(rves,zves,'r*')
plot(rl,zl,'g.')
%plot(rloo,zloo,'kx')
%plot(rpro,zpro,'bx')
axis('image')

%title(['Plasma domain and adapted to magnetic surfaces grid.',' time=',num2str(time)])
 
figure(2)
clf
hold on
colormap(jet)
mesh(rpl,zpl,curpl)
view(0,90)
%title(['Density of toroidal  current',' time=',num2str(time)])

h=plot(rb,zb,'m')
set(h,'linewidth',2)
plot(rc,zc,'r+')
%plot(rves,zves,'r*')
plot(rl,zl,'g.')
 axis('image')

%plot(rloo,zloo,'kx')
%plot(rpro,zpro,'bx')
%
%
figure(3)
clf
hold on
colormap(jet)
mesh(r,z,psi-100)
plot(rx,zx,'k+')
view(0,90)
%title(['Whole grid and plasma domain.',' time=',num2str(time)])
%
h=plot(rb,zb,'m')
set(h,'linewidth',2)
plot(rc,zc,'r+')
%plot(rves,zves,'r*')
plot(rl,zl,'g.')
%plot(rloo,zloo,'kx')
%plot(rpro,zpro,'bx')
axis('image')

%
figure(4)
clf
hold on
colormap(jet)
mesh(r,z,cur)
view(0,90)
%title(['Density of toroidal current.',' time=',num2str(time)])
%
h=plot(rb,zb,'m')
set(h,'linewidth',2)
plot(rc,zc,'r+')
%plot(rves,zves,'r*')
plot(rl,zl,'g.')
%plot(rloo,zloo,'kx')
%plot(rpro,zpro,'bx')
axis('image')

%
figure(5)
clf
hold on
colormap(jet)
mesh(r,z,cur)
view(0,90)
%title(['Density of toroidal current.',' time=',num2str(time)])
axis('image')

figure(6)
clf
hold on
colormap(jet)
contour(r,z,psie,80)
plot(rb,zb,'k')
plot(rbb,zbb,'-k.')
%plot(rves,zves,'r*')
%view(0,90)
title('psi_e_x_t')
axis('equal')

figure(7)
clf
axis('equal')
hold on
colormap(jet)
%contour(r,z,psi,80)
try
contour(r_pl,z_pl,psi_pl,20)
catch
end
plot(rb,zb,'k')
plot(rbb,zbb,'-k.')
%plot(rves,zves,'r*')
plot(rx,zx,'k+')
plot(rl,zl,'k+')
%title(['Poloidal flux contours.',' time=',num2str(time)])
axis equal
%
figure(71)
clf
axis('equal')
hold on
colormap(jet)
contour(r,z,psi,180)
%contour(r_pl,z_pl,psi_pl,20)
plot(rb,zb,'k')
plot(rbb,zbb,'-k.')
%plot(rves,zves,'r*')
plot(rx,zx,'k+')
plot(rl,zl,'k+')
%title(['Poloidal flux contours.',' time=',num2str(time)])
axis equal
%
%------------------------------------------------
%  the first and the second nearest to mag. axes
%  mag. surfaces:
%
%figure(6)
%clf
%plot(rv,zv,'m')
%hold on
%plot(rv1,zv1)
%axis('image')
%hold on
%pause  %mesh, press any key to continue
%------------------------------------------------



fid=fopen('../../../fort.4444','r'); %r z j
a1=fscanf(fid,'%f');r=a1(1:65);
z=a1(66:66+65-1);
a1=a1(66+65:end);
jg2=reshape(a1,[65 65]);
fclose(fid);
fid=fopen('../../../fort.4445','r'); %psiext
a2=fscanf(fid,'%f');
psieg2=reshape(a2,[65 65]);
fclose(fid);
fid=fopen('../../../fort.4447','r'); %psi
a3=fscanf(fid,'%f');
psig2=reshape(a3,[65 65]);
fclose(fid);









