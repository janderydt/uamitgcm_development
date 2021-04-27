function compare_outputfiles

froot1 = '/Volumes/mainJDeRydt/UaMITgcm_v2/cases/PTDC_805/mitgcm_run/input/';
nx1 = 180;  
ny1 = 360;
nz1 = 66;
nt1 = 1896;
froot2 = '/Volumes/mainJDeRydt/UaMITgcm_v2/cases/UaMITgcm_development/AS_PROPHET_999_Jan/mitgcm_run/input/';
nx2 = 240;
ny2 = 480;
nz2 = 90;
nt2 = 790;

ieee = 'b';
prec = 'real*8';

xyfiles = {};%{'bathymetry.shice','shelfice_topo.bin','pload.mdjwf'};
for ii=1:length(xyfiles)
    fid1 = fopen([froot1,xyfiles{ii}],'r');
    data1 = fread(fid1,prec,ieee);
    data1 = reshape(data1,nx1,ny1);
    fid2 = fopen([froot2,xyfiles{ii}],'r');
    data2 = fread(fid2,prec,ieee);
    data2 = reshape(data2,nx2,ny2);
    figure; hold on;
    subplot(1,2,1); pcolor(data1'); shading flat; axis equal;
    subplot(1,2,2); pcolor(data2'); shading flat; axis equal;
    sgtitle(xyfiles{ii},'interpreter','none');
end

xyzfiles = {'T_ini.bin','S_ini.bin'};
for ii=1:length(xyzfiles)
    fid1 = fopen([froot1,xyzfiles{ii}],'r');
    data1 = fread(fid1,prec,ieee);
    data1 = reshape(data1,nx1,ny1,nz1);
    fid2 = fopen([froot2,xyzfiles{ii}],'r');
    data2 = fread(fid2,prec,ieee);
    data2 = reshape(data2,nx2,ny2,nz2);
    % plot xy slice at -500m
    figure; hold on;
    subplot(1,2,1); pcolor(squeeze(data1(:,:,1))'); shading flat; axis equal;
    data2 = squeeze(data2(:,:,1));
    subplot(1,2,2); pcolor(data2'); shading flat; axis equal;
    sgtitle(xyzfiles{ii},'interpreter','none');
end

return

prec='real*4';

yztfiles = {'OBWt.bin','OBWs.bin','OBWu.bin','OBWv.bin'};
for ii=1:length(yztfiles)
    fid1 = fopen([froot1,yztfiles{ii}],'r');
    data1 = fread(fid1,prec,ieee);
    data1 = reshape(data1,ny1,nz1,nt1);
    fid2 = fopen([froot2,yztfiles{ii}],'r');
    data2 = fread(fid2,prec,ieee);
    data2 = reshape(data2,ny2,nz2,nt2);
    % plot xz slice at t=0
    figure; hold on;
    subplot(1,2,1); pcolor(squeeze(data1(:,:,1))'); shading flat; axis equal;
    caxis([-1 1]);
    data2 = squeeze(data2(:,:,1));
    subplot(1,2,2); pcolor(data2'); shading flat; axis equal;
    caxis([-1 1]);
    sgtitle(yztfiles{ii},'interpreter','none');
end

xztfiles = {};%{'OBSt.bin','OBSs.bin','OBSu.bin','OBSv.bin'};
for ii=1:length(xztfiles)
    fid1 = fopen([froot1,xztfiles{ii}],'r');
    data1 = fread(fid1,prec,ieee);
    data1 = reshape(data1,nx1,nz1,nt1);
    fid2 = fopen([froot2,xztfiles{ii}],'r');
    data2 = fread(fid2,prec,ieee);
    data2 = reshape(data2,nx2,nz2,nt2);
    % plot xz slice at t=0
    figure; hold on;
    subplot(1,2,1); pcolor(squeeze(data1(:,:,1))'); shading flat; axis equal;
    caxis([-1 1]);
    data2 = squeeze(data2(:,:,3));
    subplot(1,2,2); pcolor(data2'); shading flat; axis equal;
    caxis([-1 1]);
    sgtitle(xztfiles{ii},'interpreter','none');
end








