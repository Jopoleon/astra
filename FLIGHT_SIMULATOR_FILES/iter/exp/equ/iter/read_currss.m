fid=fopen('ppind_mat.wr');
aa=fscanf(fid,'%f');
ncoils=aa(1);
i=1+1;
for k=1:ncoils
for j=1:ncoils
try
mcoils(j,k)=aa(i)*2*pi;   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);

fid=fopen('res_mat.wr');
aa=fscanf(fid,'%f');
i=1;
for k=1:ncoils
for j=1:ncoils
try
rcoils(j,k)=aa(i)*2*pi;   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);

fid=fopen('induc_matrix.dat');
aa=fscanf(fid,'%f');
ncoils=aa(1);
i=1+1;
for k=1:ncoils
for j=1:ncoils
try
dmcoils(j,k)=aa(i);   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);

fid=fopen('resistance_matrix.dat');
aa=fscanf(fid,'%f');
i=2;
for k=1:ncoils
for j=1:ncoils
try
drcoils(j,k)=aa(i);   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);

fdsfsdsd


figure
plot(diag(mcoils)./diag(rcoils))

axis([0 116 0 1])

fshdkfsdfsd

clear all

fid=fopen('fort.78');
aa=fscanf(fid,'%f');
i=1
for k=1:10000
try
iplx(k)=aa(i:i).';
i=i+1;
ipl(k)=aa(i:i).';
i=i+1;
fpp(:,k)=aa(i:i+41-1).';
i=i+41;
catch
end
end
fclose(fid);

ncoils=13;
fid=fopen('fort.75');
aa=fscanf(fid,'%f');
i=1
for k=1:10000
try
ccoils(:,k)=aa(i:i+ncoils-1).';
i=i+ncoils;
catch
end
end
ncoils=131-13;

fid=fopen('fort.76');
aa=fscanf(fid,'%f');
i=1
for k=1:10000
try
wcoils(:,k)=aa(i:i+ncoils-1).';
i=i+ncoils;
psifb(k)=aa(i).';
i=i+1;
ipin(k)=aa(i).';
i=i+1;
ipout(k)=aa(i).';
i=i+1;
catch
end
end

fclose(fid);

fid=fopen('ppind_mat.wr');
aa=fscanf(fid,'%f');
ncoils=aa(1);
i=1+1;
for k=1:ncoils
for j=1:ncoils
try
mcoils(j,k)=aa(i)*2*pi;   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);

fid=fopen('res_mat.wr');
aa=fscanf(fid,'%f');
i=1;
for k=1:ncoils
for j=1:ncoils
try
rcoils(j,k)=aa(i)*2*pi;   % inductances from spider have a 2 pi ?
i=i+1;
catch
end
end
end
fclose(fid);






fid=fopen('psiplcoils.wr');
aa=fscanf(fid,'%f');
ncoils=aa(1);
ncoils=aa(2);
i=1+1;
i=1+1;
for j=1:ncoils
try
pmcoils(j,1)=aa(i)*2*pi;   % plasma-to-coil. to be divided by plasma current
i=i+1;
catch
end
end

fclose(fid);
i=1
fid=fopen('psicoilspl.wr');
aa=fscanf(fid,'%f');
for j=1:50
try
pcoilsm(j,1)=aa(i)*2*pi;   % coil-to-plasma. to be divided by coil current
i=i+1;
catch
end
end

fclose(fid);

fid=fopen('currents.wr');
aa=fscanf(fid,'%f');
ncoils=aa(1);
ncoilz=aa(2);
i=1+1;
i=3;
for j=1:ncoilz
try
currenz(j,1)=aa(i);   % plasma-to-coil. to be divided by plasma current
i=i+1;
catch
end
end

fclose(fid);

coilress=[  22400.  ...
               5000  ... 
                  5000   ...
                 8680   ...
              8680   ...
                 23450   ...
                  23450   ...
                 17180   ...
                550    ...
                550   ...
                  13.6    ...
                 13.6  280/118*[13:130]];

tges=[1.85 1.68 1.68 3.95 3.95 3.66 3.66 1.51 0.53 0.53 0 0 0*[13:130]];

coilw=[530 81 81 35+31+27 35+31+27 39+47 39+47 28 28 5 5 1 1 1*[13:130]];


