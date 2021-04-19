function gen_mesh

x = ncread('MISOMIP2_OceanA2_Geometries.nc','x');
y = ncread('MISOMIP2_OceanA2_Geometries.nc','y');

xright = max(x);
xleft = min(x);
ybot = min(y);
ytop = max(y);

PAS = 851

% choose grid size, as well as domain decomposition
ny = 432;
nx = 216;
npx = 8;
npy = 18;

% this assumes that x,y from netcdf file are corners of cells -- we want coord centres
x_mesh = linspace(min(x),max(x),nx+1);
y_mesh = linspace(min(y),max(y),ny+1);
x_mesh_mid = .5 * (x_mesh(1:end-1)+x_mesh(2:end));
y_mesh_mid = .5 * (y_mesh(1:end-1)+y_mesh(2:end));
diffx = diff(x_mesh);
diffy = diff(y_mesh);

[x_mesh_mid y_mesh_mid] = meshgrid(x_mesh_mid,y_mesh_mid);

lendiff_x = ceil(nx/npx)*npx - nx
lendiff_y = ceil(ny/npy)*npy - ny

save griddata.mat x_mesh_mid y_mesh_mid 

disp(['mesh size: ' num2str(nx + lendiff_x) ' by ' num2str(ny + lendiff_y)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/exports/csce/datastore/geos/users/dgoldber/ice_data/ThwaitesDataProphet/CODE');
density_ice = 917;
density_oce = 1027;
air_depth = 0;
X = x_mesh_mid;
Y = y_mesh_mid;

bed = interpBedmachineAntarctica(X,Y,'bed','2020-07-15');
thick = interpBedmachineAntarctica(X,Y,'thickness','2020-07-15');
surf = interpBedmachineAntarctica(X,Y,'surface','2020-07-15');
mask = interpBedmachineAntarctica(X,Y,'mask','2020-07-15');

% following creates a smoothed version of surface elevation (but conserving
% ice free cells)
surf2 = surf;
surf2(thick<5 & surf2>0)=.5; 
surf2(surf2==0 | mask==1 | mask==0)=nan;
gaussFilter = fct_GaussianFilter([2 2], 1, 0);
[surf2,im_conv,count,NaNcount] = fct_convNaN(surf2, gaussFilter, 'same', .5);
surf2(mask==0 | mask==1)=nan;
surf = surf2;
surf(isnan(surf))=0;

% just to make simpler to define nz
bed = max(bed,-1800);

% rather than use bedmachine thickness we ensure our thickness is
% hydrostatically balanced with the choices of density
thick_floatation = (density_ice * air_depth - density_oce*surf) / (density_ice - density_oce);
thick_floatation(surf==0)=0;
base_floatation = surf - thick_floatation;
base = max(bed,base_floatation);
thick = surf-base;
thick(surf==0)=0;

% add overflow regions to binary files to fit on processors
gx = lendiff_x
gy = lendiff_y
base = [[base zeros(ny,gx)];zeros(gy,nx+gx)];
binwrite('icetopo_init.bin',base');
thick = [[thick zeros(ny,gx)];zeros(gy,nx+gx)];
binwrite('icethick.bin',thick');
bed = [[bed zeros(ny,gx)];zeros(gy,nx+gx)];
binwrite('topog.bin',bed');


disp('low bed')
min(min(bed))

%readDataLarge_pahol(PAS,false,1996,2015,...
%                    true,true,true);
%rdmds_init(nx,ny,gx,gy,PAS,0);
gendata(nx,ny,gx,gy,PAS,0);

disp(['tot_x: ' num2str(gx+nx)]);
disp(['tot_y: ' num2str(gy+ny)]);

binwrite('delX.bin',[diffx ones(1,gx)]);
binwrite('delY.bin',[diffy ones(1,gy)]);
