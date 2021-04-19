%Verion of gendata.m modified by Vero
%This is a matlab script that generates the input data


% the configuation approximately the ISOMIP experiment no. 1
% require matlab functions for equation of state

function gendata(nx_in,ny_in,gx,gy,PAS,mincol)
% Dimensions of grid
nx=nx_in; 
ny=ny_in;
nz=72;
delz = 25;
prec = 8;
thresh = 0.2

rho_ice = 917;

hfacMin = 0.04;
mwct = hfacMin * delz * 2;

% dlat = 0.125/4; dy=dlat;
% dlon = 0.125; dx=dlon;

eos = 'jmd95z';

acc = 'real*8';

bathy = binread('topog.bin',8,nx+gx,ny+gy);

dz = delz*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);

zmid = zc;

T0 = binread(['theta.init.' num2str(PAS)],prec,nx+gx,ny+gy,nz);
S0 = binread(['salt.init.' num2str(PAS)],prec,nx+gx,ny+gy,nz);
tref = squeeze(T0(1,1,:));
sref = squeeze(S0(1,1,:));



rhoConst = 1030;
%talpha = 2e-4;
%sbeta  = 7.4e-4;

t    = tref;
s    = sref;
gravity = 9.81;
k=1;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
p = abs(zc)*gravity*rhoConst*1e-4;
dp = p;
kp = 0;

while rms(dp) > 1e-13
  phiHydF(k) = 0;
  p0 = p;
  kp = kp+1
  for k = 1:nz
    switch eos
     case 'linear'
      drho = rhoConst*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k)))-rhoConst;
     case 'jmd95z'
      drho = densjmd95(s(k),t(k),p(k))-rhoConst;
      rho0(k) = drho;
     case 'mdjwf'
      drho = densmdjwf(s(k),t(k),p(k))-rhoConst;
     otherwise
      error(sprintf('unknown EOS: %s',eos))
    end
    phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
    phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
  end
  switch eos
   case 'mdjwf'
    p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
  end
  dp = p-p0;
end

p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity;

fid = fopen('icethick.bin','r','b');
thick = fread(fid,inf,'real*8');
fclose(fid);
thick=(reshape(thick, [nx+gx ny+gy]));
shelficemass=thick*rho_ice;


topo = zeros(nx+gx,ny+gy);

bathy(bathy>-delz)=0;
bathy2 = bathy;

% ensure that there will no adjustments to bathy due to hfacmin
for ix=1:nx
  for iy=1:ny
          cellFace = mod(abs(bathy2(ix,iy)),delz);
          if (cellFace<hfacMin*delz & cellFace>0);
              bathy2(ix,iy) = -floor(abs(bathy2(ix,iy))/delz)*delz;
          end
  end
end
bathy = bathy2;

% find equilibrium ice shelf topo corresponding to mass

for ix=1:(nx+gx)
  for iy=1:(ny+gy)

     mass = shelficemass (ix,iy);
     massFuncC = rhoConst * (phiHydC/gravity + zc);
    
     massFuncF = rhoConst * (phiHydF/gravity + zgp1);

     k = max (find ( massFuncF < mass ));
     if (k>nz)
         topo(ix,iy) = bathy(ix,iy)+mwct;
     else
        if (isempty(k))
         k=0;
        end
         
     
        if (k>0)
         if (mass < massFuncC(k))
          ztopo = -zg(k) - (mass-massFuncF(k)) * delz/2 / (massFuncC(k)-massFuncF(k));
          topo(ix,iy) = max(ztopo,bathy(ix,iy)+mwct);
         else
          ztopo = -zc(k) - (mass-massFuncC(k)) * delz/2 / (massFuncF(k+1)-massFuncC(k));
          topo(ix,iy) = max(ztopo,bathy(ix,iy)+mwct);
         end
        end
     end
     
  end
end

etainit = zeros(size(topo));

% new topography: icetopo rounded to the nearest k * deltaZ
%                 eta_init set to make difference


icetopo2 = topo;

for ix=1:nx
  for iy=1:ny
    k=max(find(abs(zg)<abs(icetopo2(ix,iy))));
    if isempty(k)
      k=0;
    else
      
      dr = 1-(-zg(k) - icetopo2(ix,iy))/delz;

      if (dr > thresh | ((-zg(k+1)-bathy(ix,iy))<hfacMin))
          if (k==1);
          % top level -- cannot make R_shelfice=0. Must leave in between
          % cell faces until a merge down occurs
           icetopo2(ix,iy) = icetopo2(ix,iy);
           etainit(ix,iy) = 0;
          else;
          % bring Ro_surf *up* to closest grid face & make etainit negative
          % to compensate
           icetopo2(ix,iy) = -zg(k);
           etainit(ix,iy) = (dr-1)*delz;
          end
      else
          % bring Ro_surf *down* to closest grid face & make etainit pos
          % to compensate
          icetopo2(ix,iy) = -zg(k+1);
          etainit(ix,iy) = (dr)*delz;
      end

       
    end
  end
end

load griddata

ind = min(find(x_mesh_mid(1,:)>-1.6e6)); 
bathy(ind:end,1) = 0;
disp(['south ind: ' num2str(ind-1)]);

ind = min(find(y_mesh_mid(:,1)>-3.435e5)); 
bathy(1,ind:end) = 0;
disp(['west ind: ' num2str(ind-1)]);

bathy(:,end)=0;
bathy(end,:)=0;

if (mincol>0);
    
    col = icetopo2+etainit-bathy;
    bathy(bathy<0 & col<mincol) = 0;
    
    con_mask = bathy<0 ;
    CC = bwconncomp(con_mask,4);
    pp = CC.PixelIdxList;
    if (length(pp)>1);
     for i=2:length(pp);
      disp(['length: ' num2str(length(pp{i}))]);
      for k=1:length(pp{i});
       bathy(pp{i}(k)) = 0;
      end;
     end;
    end;

end

% need to be careful not to "trap" inflow at bdries




fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,icetopo2,'real*8'); fclose(fid);
fid = fopen('etainit.bin','w','b'); fwrite(fid,etainit,'real*8'); fclose(fid);
fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,shelficemass,'real*8'); fclose(fid);
fid = fopen('bathy_mod.bin','w','b'); fwrite(fid,bathy,'real*8'); fclose(fid);


