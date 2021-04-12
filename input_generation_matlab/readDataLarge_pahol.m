function readDataLarge_pahol (PAS,loop,start_year,end_year,...
    calc_init,calc_bounds, calc_surf)

if (ispc)
    if (PAS>0)
     nc = ['Z:/dgoldber/pahol_output/PAS_' num2str(PAS) '/run/'];
    else
     nc = ['R:/ice_data/toshi_data/ASE'];
    end
else
    if (PAS>0)
     nc = ['/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/pahol_output/PAS_' num2str(PAS) '/run/'];
    else
     nc = ['/home/dgoldber/ice_data/toshi_data/ASE/']
    end
end

lon = double(ncread([nc 'stateUvel.nc'],'LONGITUDE'));
lat = double(ncread([nc 'stateUvel.nc'],'LATITUDE'));
z = double(ncread([nc 'stateUvel.nc'],'DEPTH'));

load   griddata

delx = x_mesh_mid(2)-x_mesh_mid(1);
dely = y_mesh_mid(2)-y_mesh_mid(1);
x_mesh_mid = x_mesh_mid(1,:)';
y_mesh_mid = y_mesh_mid(:,1);

x_botbdry = x_mesh_mid;
y_botbdry = (y_mesh_mid(1)) * ones(length(x_mesh_mid),1);
x_lbdry = (x_mesh_mid(1)) * ones(length(y_mesh_mid),1);
y_lbdry = y_mesh_mid;
x_topbdry = x_mesh_mid;
y_topbdry = (y_mesh_mid(end)) * ones(length(x_mesh_mid),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (PAS>0)
 if (~loop);
    start_month = (start_year - 1955) * 12 
    end_month = (end_year - 1955 + 1) * 12 + 1
 else
    start_month = (start_year - 1955) * 12 + 1
    end_month = (end_year - 1955 + 1) * 12 
 end
else
 if (~loop); 
    start_month = (start_year - 1991) * 12
    end_month = (end_year - 1991 + 1) * 12 + 1
 else
    start_month = (start_year - 1991) * 12 + 1
    end_month = (end_year - 1991 + 1) * 12 
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xgridp ygridp] = meshgrid(x_mesh_mid,y_mesh_mid);
if (calc_init)

% initialise T/S by averaging over 1st year
% will be averaged horizontally as well in rdmds_init.m
    


[phi,lambda]=polarstereo_inv(xgridp,ygridp,[],[],-71,0);
lambda = lambda + 360;

nSz = length(x_mesh_mid);
mSz = length(y_mesh_mid);

ntot = 12;

Tinit = zeros(mSz,nSz,length(z));
Sinit = zeros(mSz,nSz,length(z));
Uinit = zeros(mSz,nSz,length(z));
Vinit = zeros(mSz,nSz,length(z));

[latInterp lonInterp] = meshgrid(lat,lon);

 for n=1:ntot;
    
  start_month-1+n
  Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n],[600 384 50 1]));
  Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n],[600 384 50 1]));
  Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n],[600 384 50 1]));
  Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n],[600 384 50 1]));
  disp('got here')

  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
  
  for k=1:length(z)
  Tinit(:,:,k) = Tinit(:,:,k) + 1/ntot * interp2(lat,lon,Ttemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  Sinit(:,:,k) = Sinit(:,:,k) + 1/ntot * interp2(lat,lon,Stemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  tlayer = Tinit(:,:,k);
  I = ~isnan(tlayer);
  disp(['month ' num2str(n) '; level ' num2str(k) ' ' num2str(mean(tlayer(I)))]);
  end
  
  
  
  
 end

save initpahol.mat Sinit Tinit z xgridp ygridp


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[latps lonps] = meshgrid(lat,lon);


if (calc_surf);
    
   [phi,lambda]=polarstereo_inv(xgridp,ygridp,[],[],-71,0);
   lambda = lambda + 360;
    
   ntot = end_month-start_month+1
   Hflux = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   Sflux = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   Utau = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   Vtau = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   SSS = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   SST = zeros(length(x_mesh_mid),length(y_mesh_mid),end_month-start_month+1);
   
   for n=1:ntot;
    if (PAS==0 & n>ntot-2)
     n2 = ntot-2;
    else
     n2 = n;
    end
    n
    HfluxTmp = double(ncread([nc 'state2D.nc'],'oceQnet',[1 1 start_month-1+n2],[600 384 1]));
    Sfluxtmp = double(ncread([nc 'state2D.nc'],'oceFWflx',[1 1 start_month-1+n2],[600 384 1]));
    Utautmp = double(ncread([nc 'state2D.nc'],'oceTAUX',[1 1 start_month-1+n2],[600 384 1]));
    Vtautmp = double(ncread([nc 'state2D.nc'],'oceTAUY',[1 1 start_month-1+n2],[600 384 1]));
    ssstmp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 1 1]));
    ssttmp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 1 1]));
    
    
    tmp = interp2(lat,lon,HfluxTmp,phi,lambda)'; tmp(isnan(tmp))=0;
    Hflux(:,:,n) = tmp;
    
    tmp = interp2(lat,lon,Sfluxtmp,phi,lambda)'; tmp(isnan(tmp))=0;
    Sflux(:,:,n) = tmp;
    
    tmp = interp2(lat,lon,ssstmp,phi,lambda)'; tmp(isnan(tmp) | tmp==0)=34;
    SSS(:,:,n) = tmp;
    
    tmp = interp2(lat,lon,ssttmp,phi,lambda)'; tmp(isnan(tmp) | tmp==0)=0;
    SST(:,:,n) = tmp;
    
    [xdummy ydummy taux tauy] = vec_ll2ps(lonps,latps,Utautmp,Vtautmp,[],[]);
    tmp = interp2(lat,lon,taux,phi,lambda)'; tmp(isnan(tmp))=0;
    Utau(:,:,n) = tmp;
    tmp = interp2(lat,lon,tauy,phi,lambda)'; tmp(isnan(tmp))=0;
    Vtau(:,:,n) = tmp;
   end
    
   
   save surfaceDatapahol.mat Hflux Sflux Utau Vtau SSS SST
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (calc_bounds)

bdryX = x_botbdry; bdryY = y_botbdry;
Tbotbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Sbotbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Ubotbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Vbotbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
[phiLeft,lambdaLeft]=polarstereo_inv(bdryX,bdryY,[],[],-71,0);
lambdaLeft = lambdaLeft + 360;

ntot = end_month-start_month+1

for n=1:ntot;
  if (PAS==0 & n>ntot-2)
   n2 = ntot-2;
  else
   n2 = n;
  end
  Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
      
  for k=1:length(z);
      
      Tbotbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
      Sbotbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);

      [xdummy ydummy Ups Vps] = vec_ll2ps(lonps,latps,Utemp(:,:,k),Vtemp(:,:,k),[],[]);
      Ubotbdry(k,:,n) = interp2(lat,lon,Ups,phiLeft,lambdaLeft);
      Vbotbdry(k,:,n) = interp2(lat,lon,Vps,phiLeft,lambdaLeft);
      
      
      disp(['bot ' num2str([k n Tbotbdry(k,1,n)])])
  end
  
end


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bdryX = x_lbdry; bdryY = y_lbdry;
Tleftbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Sleftbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Uleftbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
Vleftbdry = zeros(length(z),length(bdryX),end_month-start_month+1);
[phiLeft,lambdaLeft]=polarstereo_inv(bdryX,bdryY,[],[],-71,0);
lambdaLeft = lambdaLeft + 360;


ntot = end_month-start_month+1;

for n=1:ntot;
  if (PAS==0 & n>ntot-2)
   n2 = ntot-2;
  else
   n2 = n;
  end
  Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
      
  for k=1:length(z);
      
      Tleftbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
      Sleftbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);
      [xdummy ydummy Ups Vps] = vec_ll2ps(lonps,latps,Utemp(:,:,k),Vtemp(:,:,k),[],[]);
      Uleftbdry(k,:,n) = interp2(lat,lon,Ups,phiLeft,lambdaLeft);
      Vleftbdry(k,:,n) = interp2(lat,lon,Vps,phiLeft,lambdaLeft);
      
      disp(['left ' num2str([k n Tleftbdry(k,1,n)])])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% bdryX = x_topbdry; bdryY = y_topbdry;
% Ttopbdry = zeros(length(z),length(bdryX),end_month-start_month+3);
% Stopbdry = zeros(length(z),length(bdryX),end_month-start_month+3);
% Utopbdry = zeros(length(z),length(bdryX),end_month-start_month+3);
% Vtopbdry = zeros(length(z),length(bdryX),end_month-start_month+3);
% [phiLeft,lambdaLeft]=polarstereo_inv(bdryX,bdryY,[],[],-71,0);
% lambdaLeft = lambdaLeft + 360;
% 
% gradLat = zeros(length(phiLeft),2);
% for i=1:length(phiLeft);
%     [latTempA lonTemp] = polarstereo_inv(bdryX(i)+100,bdryY(i),[],[],-71,0);
%     [latTempB lonTemp] = polarstereo_inv(bdryX(i)-100,bdryY(i),[],[],-71,0);
%     [latTempC lonTemp] = polarstereo_inv(bdryX(i),bdryY(i)+100,[],[],-71,0);
%     [latTempD lonTemp] = polarstereo_inv(bdryX(i),bdryY(i)-100,[],[],-71,0);
%     gradLat(i,1) = (latTempA-latTempB)/200;
%     gradLat(i,2) = (latTempC-latTempD)/200;
%     gradLat(i,:) = gradLat(i,:) / sqrt(sum(gradLat(i,:).^2));
% end
% gradLon = [gradLat(:,2) -gradLat(:,1)];
% 
% ntot = end_month-start_month+1;
% 
% for n=1:0;
%   Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n],[600 384 50 1]));
%   Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n],[600 384 50 1]));
%   Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n],[600 384 50 1]));
%   Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n],[600 384 50 1]));
%   
%   Ttemp(Stemp==0)=nan;
%   Vtemp(Stemp==0)=nan;
%   Utemp(Stemp==0)=nan;
%   Stemp(Stemp==0)=nan;
%       
%   for k=1:length(z);
% %     for i=1:length(bdryX);
%       
%       Ttopbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
%       Stopbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);
%       u = interp2(lat,lon,Utemp(:,:,k),phiLeft,lambdaLeft);
%       v = interp2(lat,lon,Vtemp(:,:,k),phiLeft,lambdaLeft);
%       
%       uPolar = u .* gradLon(:,1) + v .* gradLat(:,1);
%       vPolar = u .* gradLon(:,2) + v .* gradLat(:,2);
%       Utopbdry(k,:,n) = uPolar;
%       Vtopbdry(k,:,n) = vPolar;
%       
%       disp(['left ' num2str([k n Ttopbdry(k,1,n)])])
% %     end
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save bdryDatapahol.mat Tleftbdry Sleftbdry Uleftbdry Vleftbdry ... 
     Tbotbdry Sbotbdry Ubotbdry Vbotbdry ...
     z

end
