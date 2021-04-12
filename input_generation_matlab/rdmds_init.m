%This program creates open boundary prescription files for the PIG
%experiment based on the 10 yrs spin-up run using OBCS with 
%U,V = 0, Tref, Sref

%following variables 
%-----3D fields-----
% T Temperature (C)
% S Salinity (psu)
% U u-velocity (m/s)
% PH ocean pressure (or atm geopotential)

function rdmds_init(nx_in,ny_in,gx,gy,PAS,smoothnum)

close all

%Set the grid size;
nx = nx_in;    delx = 1;   X = nx*delx;
ny = ny_in;    dely = 1;   Y = ny*dely;
nz = 72;     delz = 25;  Z = nz*delz;

rhoConst = 1030;

load griddata.mat

%x axis 
x = x_mesh_mid(1,:);
%y axis
y = y_mesh_mid(:,1);
%z axis [m]
z = 0:delz:Z;
zmid = .5 * (z(1:end-1)+z(2:end));

dzmod = 25*ones(1,nz);
zface = [0 cumsum(dzmod)];
zmid = .5 * (zface(1:end-1)+zface(2:end));

%Note: 
%Ensure that the volume flux is zero at the open boundary

%% Print init files for T, S


z2 = 0:delz:Z;
z2 = zface;

load initpahol.mat;

T_init = zeros(ny,nx,nz);
S_init = zeros(ny,nx,nz);

for i=1:ny; 
    for j=1:nx; 
        T_init(i,j,:) = interp1(z,squeeze(Tinit(i,j,:)),zmid);
        S_init(i,j,:) = interp1(z,squeeze(Sinit(i,j,:)),zmid);
    end
end

T_init_prof = zeros(1,1,nz);
S_init_prof = zeros(1,1,nz);
for k=1:nz;
    Tlayer = T_init(:,:,k);
    Slayer = S_init(:,:,k);
    T_init_prof(k) = mean(Tlayer(~isnan(Tlayer)));
    S_init_prof(k) = mean(Slayer(~isnan(Slayer)));
end

size(T_init_prof)
squeeze(T_init_prof)
size(isnan(T_init_prof))

T_init_prof(isnan(T_init_prof))=T_init_prof(1,1,max(find(~isnan(squeeze(T_init_prof)))));
S_init_prof(isnan(S_init_prof))=S_init_prof(1,1,max(find(~isnan(squeeze(S_init_prof)))));

T_init  = repmat(T_init_prof,[nx+gx ny+gy 1]);
S_init  = repmat(S_init_prof,[nx+gx ny+gy 1]);


fid_T=fopen(['theta.init.' num2str(PAS)],'w','b');fwrite(fid_T,permute(T_init,[1 2 3]),'real*8');fclose(fid_T);
fid_T=fopen(['salt.init.' num2str(PAS)],'w','b');fwrite(fid_T,permute(S_init,[1 2 3]),'real*8');fclose(fid_T);


load bdryDatapahol.mat

n_months = size(Sbotbdry,3);

z=double(-z);
z2=double(-z2);

S_south = Sbotbdry; S_south = permute(S_south,[2 1 3]);
T_south = Tbotbdry; T_south = permute(T_south,[2 1 3]);
V_south = Vbotbdry; V_south = permute(V_south,[2 1 3]);
U_south = Ubotbdry; U_south = permute(U_south,[2 1 3]);

S_south(120:end,:,:)=nan;
T_south(120:end,:,:)=nan;
V_south(120:end,:,:)=nan;
U_south(120:end,:,:)=nan;

%should use nearest neighbour not S_col and T_col... will cause spurious
%t/s and convection

%Do this instead;
S_Int = zeros(nx+gx,nz,n_months);
T_Int = zeros(nx+gx,nz,n_months);
U_Int = zeros(nx+gx,nz,n_months);
V_Int = zeros(nx+gx,nz,n_months);
[X,Y] = meshgrid(linspace(1,nx,nx),z); Y=Y';X=X';
zmid = .5 * (z2(1:end-1)+z2(2:end));
[X2,Y2] = meshgrid(linspace(1,nx,nx),zmid); Y2=Y2';X2=X2';
for t=1:n_months
    S_time = smooth_obc(S_south(:,:,t),smoothnum);
    T_time = smooth_obc(T_south(:,:,t),smoothnum);
    U_time = smooth_obc(U_south(:,:,t),smoothnum);
    V_time = smooth_obc(V_south(:,:,t),smoothnum);
    ibad = isnan(S_time);
    S_interp = griddata(X(~ibad), Y(~ibad), S_time(~ibad), X2, Y2, 'nearest'); S_interp(isnan(S_interp))=-9999;
    S_Int(1:nx,:,t) = S_interp;
    T_interp = griddata(X(~ibad), Y(~ibad), T_time(~ibad), X2, Y2, 'nearest'); T_interp(isnan(T_interp))=-9999;
    T_Int(1:nx,:,t) = T_interp;
    U_interp = griddata(X(~ibad), Y(~ibad), U_time(~ibad), X2, Y2, 'nearest'); U_interp(isnan(U_interp))=-9999;
    U_Int(1:nx,:,t) = U_interp;
    
    V_interp = griddata(X(~ibad), Y(~ibad), V_time(~ibad), X2, Y2, 'nearest'); V_interp(isnan(V_interp))=-9999;
    V_Int(1:nx,:,t) = V_interp;
    
    t
end


binwrite(['vvel_' num2str(PAS)  '.obs'],V_Int);
binwrite(['uvel_' num2str(PAS)  '.obs'],U_Int);
binwrite(['salt_' num2str(PAS)  '.obs'],S_Int);
binwrite(['temp_' num2str(PAS)  '.obs'],T_Int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_west = Sleftbdry; S_west = permute(S_west,[2 1 3]);
T_west = Tleftbdry; T_west = permute(T_west,[2 1 3]);
U_west = Uleftbdry; U_west = permute(U_west,[2 1 3]);
V_west = Vleftbdry; V_west = permute(V_west,[2 1 3]);

[X,Y] = meshgrid(linspace(1,ny,ny),z); Y=Y';X=X';
[X2,Y2] = meshgrid(linspace(1,ny,ny),zmid); Y2=Y2';X2=X2';
S_Int = zeros(ny+gy,nz,n_months);
T_Int = zeros(ny+gy,nz,n_months);
U_Int = zeros(ny+gy,nz,n_months);
V_Int = zeros(ny+gy,nz,n_months);
for t=1:n_months
    S_time = smooth_obc(S_west(:,:,t),smoothnum);
    T_time = smooth_obc(T_west(:,:,t),smoothnum);
    U_time = smooth_obc(U_west(:,:,t),smoothnum);
    V_time = smooth_obc(V_west(:,:,t),smoothnum);
    ibad = isnan(S_time);
    S_interp = griddata(X(~ibad), Y(~ibad), S_time(~ibad), X2, Y2, 'nearest'); S_interp(isnan(S_interp))=-9999;
    S_Int(1:ny,:,t) = S_interp;
    T_interp = griddata(X(~ibad), Y(~ibad), T_time(~ibad), X2, Y2, 'nearest'); T_interp(isnan(T_interp))=-9999;
    T_Int(1:ny,:,t) = T_interp;
    U_interp = griddata(X(~ibad), Y(~ibad), U_time(~ibad), X2, Y2, 'nearest'); U_interp(isnan(U_interp))=-9999;
    U_Int(1:ny,:,t) = U_interp; 
    V_interp = griddata(X(~ibad), Y(~ibad), V_time(~ibad), X2, Y2, 'nearest'); V_interp(isnan(V_interp))=-9999;
    V_Int(1:ny,:,t) = V_interp;
    t
end


binwrite(['vvel_' num2str(PAS)  '.obw'],V_Int);
binwrite(['uvel_' num2str(PAS)  '.obw'],U_Int);
binwrite(['salt_' num2str(PAS)  '.obw'],S_Int);
binwrite(['temp_' num2str(PAS)  '.obw'],T_Int);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load surfaceDatapahol.mat

Hflx = zeros(nx+gx,ny+gy,n_months);
Sflx = zeros(nx+gx,ny+gy,n_months);
taux = zeros(nx+gx,ny+gy,n_months);
tauy = zeros(nx+gx,ny+gy,n_months);
sst = zeros(nx+gx,ny+gy,n_months);
sss = zeros(nx+gx,ny+gy,n_months);


% need to ensure that all fields are zero where there is ice mass
thick = binread('icethick.bin',8,nx+gx,ny+gy);

for i=1:n_months;
    Hflx(1:nx,1:ny,i) = Hflux(:,:,i).*double(thick(1:nx,1:ny)<=0);

% oceFWflux is kg/m^2/s DOWN. EXF input is m/s UP.

    Sflx(1:nx,1:ny,i) = (-Sflux(:,:,i)/rhoConst).*double(thick(1:nx,1:ny)<=0);
    taux(1:nx,1:ny,i) = Utau(:,:,i).*double(thick(1:nx,1:ny)<=0);
    tauy(1:nx,1:ny,i) = Vtau(:,:,i).*double(thick(1:nx,1:ny)<=0);
    sst(1:nx,1:ny,i) = SST(:,:,i).*double(thick(1:nx,1:ny)<=0);
    sss(1:nx,1:ny,i) = SSS(:,:,i).*double(thick(1:nx,1:ny)<=0);
    
end

binwrite(['Hflux_' num2str(PAS) '.bin'],Hflx);
binwrite(['Sflux_' num2str(PAS) '.bin'],Sflx);
binwrite(['ocetaux_' num2str(PAS) '.bin'],taux);
binwrite(['ocetauy_' num2str(PAS) '.bin'],tauy);
binwrite(['sss_' num2str(PAS) '.bin'],sss);
binwrite(['sst_' num2str(PAS) '.bin'],sst);


%%%%%%%%%%%%%%%%%%%

if (false)
% S_north = Stopbdry; S_north = permute(S_north,[2 1 3]);
% T_north = Ttopbdry; T_north = permute(T_north,[2 1 3]);
% V_north = Vtopbdry; V_north = permute(V_north,[2 1 3]);
% U_north = Utopbdry; U_north = permute(U_north,[2 1 3]);
% 
% %should use nearest neighbour not S_col and T_col... will cause spurious
% %t/s and convection
% 
% %Do this instead;
% S_Int = zeros(nx+gx,nz,n_months);
% T_Int = zeros(nx+gx,nz,n_months);
% U_Int = zeros(nx+gx,nz,n_months);
% V_Int = zeros(nx+gx,nz,n_months);
% [X,Y] = meshgrid(linspace(1,nx,nx),z); Y=Y';X=X';
% zmid = .5 * (z2(1:end-1)+z2(2:end));
% [X2,Y2] = meshgrid(linspace(1,nx,nx),zmid); Y2=Y2';X2=X2';
% for t=1:n_months
%     S_time = S_north(:,:,t);
%     T_time = T_north(:,:,t);
%     U_time = U_north(:,:,t);
%     V_time = V_north(:,:,t);
%     ibad = isnan(S_time);
%     S_interp = griddata(X(~ibad), Y(~ibad), S_time(~ibad), X2, Y2, 'nearest');
%     S_Int(1:nx,:,t) = S_interp;
%     T_interp = griddata(X(~ibad), Y(~ibad), T_time(~ibad), X2, Y2, 'nearest');
%     T_Int(1:nx,:,t) = T_interp;
%     U_interp = griddata(X(~ibad), Y(~ibad), U_time(~ibad), X2, Y2, 'nearest'); U_interp(isnan(U_interp))=0;
%     U_Int(1:nx,:,t) = U_interp;
%     V_interp = griddata(X(~ibad), Y(~ibad), V_time(~ibad), X2, Y2, 'nearest'); V_interp(isnan(V_interp))=0;
%     V_Int(1:nx,:,t) = V_interp;
%     t
% end
end

