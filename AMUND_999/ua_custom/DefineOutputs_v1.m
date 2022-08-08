function UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)


v2struct(F);
time=CtrlVar.time; 
    
% check if folder 'UaOutputsDirectory' exists, if not create
if exist(fullfile(cd,UserVar.UaMITgcm.UaOutputDirectory),'dir')~=7
    mkdir(UserVar.UaMITgcm.UaOutputDirectory) ;
end

% write output in matlab or netcdf format
if strcmp(UserVar.UaMITgcm.UaOutputFormat,'matlab')

    DaysSinceStart = num2str(sprintf('%04d',round(time*365.25)));
    FileName=sprintf([UserVar.UaMITgcm.UaOutputDirectory,'/',CtrlVar.Experiment,'_',...
        UserVar.UaMITgcm.StartYear,UserVar.UaMITgcm.StartMonth,'_',DaysSinceStart]);
    fprintf(' Saving data in %s \n',FileName);
    save(FileName,'UserVar','CtrlVar','MUA','F','GF');
    
elseif strcmp(UserVar.UaMITgcm.UaOutputFormat,'netcdf')
    
    %% add here for netcdf output - use MISOMIP code
    
else
    
    error('Unknown case for writing Ua output');
    return
    
end

%% initialize coordinates
lonCMIT = UserVar.UaMITgcm.MITgcmCGridlon; % 2d arrays
latCMIT = UserVar.UaMITgcm.MITgcmCGridlat;
dlat = latCMIT(1,2)-latCMIT(1,1); dlon =  lonCMIT(2,1)-lonCMIT(1,1);
[nx,ny] = size(lonCMIT);
[UaBlat,UaBlon]=psxy2ll(MUA.Boundary.x,MUA.Boundary.y,-71,0);

%% Generate b and B fields for MITgcm
% Obtain bed via linear interpolation from Ua B field
FB = scatteredInterpolant(x,y,F.B,'linear');
B_forMITgcm = FB(UserVar.UaMITgcm.MITgcmCGridX(:),UserVar.UaMITgcm.MITgcmCGridY(:));
% For MITgcm tracer points outside Ua domain replace via DefineGeometry
I_outsideUa = find(~inpoly([UserVar.UaMITgcm.MITgcmCGridX(:),UserVar.UaMITgcm.MITgcmCGridY(:)],[MUA.Boundary.x(:),MUA.Boundary.y(:)]));
MITmesh.coordinates = [UserVar.UaMITgcm.MITgcmCGridX(I_outsideUa) UserVar.UaMITgcm.MITgcmCGridY(I_outsideUa)];
[~,~,~,~,B_forMITgcm(I_outsideUa),~]=DefineGeometry(UserVar,CtrlVar,MITmesh,[],'B');

% We use a linear interpolation to map the Ua draft onto the MITgcm grid. Note that more sophisticated
% methods can be implemented, such as 'data binning'. If the MITgcm tracer points are a subset of the Ua nodes then
% interpolation is not required. 
Fb = scatteredInterpolant(x,y,F.b,'linear');
b_forMITgcm = Fb(UserVar.UaMITgcm.MITgcmCGridX(:),UserVar.UaMITgcm.MITgcmCGridY(:));
b_forMITgcm(I_outsideUa) = 0;

%% Generate mask for MITgcm: 
%% =====THIS IS SPECIFIC TO THE AMUND SETUP=====%%

% generate triangulation of MITgcm points
MITgcm_connectivity = delaunay(lonCMIT(:),latCMIT(:));

%% First identify open ocean nodes
% temporary mask: 0 for nodes outside Ua domain, 1 inside Ua boundary
I_BL = inpoly([lonCMIT(:)-dlon/2 latCMIT(:)-dlat/2],[UaBlon(:) UaBlat(:)]); 
I_TL  = inpoly([lonCMIT(:)-dlon/2 latCMIT(:)+dlat/2],[UaBlon(:) UaBlat(:)]); 
I_TR  = inpoly([lonCMIT(:)+dlon/2 latCMIT(:)+dlat/2],[UaBlon(:) UaBlat(:)]); 
I_BR = inpoly([lonCMIT(:)+dlon/2 latCMIT(:)-dlat/2],[UaBlon(:) UaBlat(:)]); 
I_outsideUa = find(I_BL+I_TL+I_TR+I_BR==0);
Mask_tmp = 0*lonCMIT+1;
Mask_tmp(I_outsideUa) = 0;
% for graph, only retain MITgcm elements that completely fall outside Ua domain boundary
Mask_ele = sum(Mask_tmp(MITgcm_connectivity),2);
TRI = MITgcm_connectivity(find(Mask_ele==0),:);
TRI_OO = TRI;
% elements crossing the ice front:
Iele_crossingicefront = find(Mask_ele>0 & Mask_ele<3);

% create undirected graph
G = graph(TRI,TRI(:,[2 3 1]));
% calculate the connected components of the graph
bins=conncomp(G) ;

% define indices of OceanBoundaryNodes at northern boundary for flooding
OceanBoundaryNodes = find(latCMIT == max(latCMIT(:)));

% initialise arrays
Nnum = zeros(nx*ny,1);
OpenOceanNodes = zeros(nx*ny,1);
Nnum(OceanBoundaryNodes) = 1;

% loop through ocean boundary nodes until each one has been checked for
% connected floating nodes, once this is done for all boundary nodes the
% only floating nodes left should be lakes
while sum(Nnum)>0
    
    NodeSeed = find(Nnum,1,'first');
    ID=bins(NodeSeed) ;
    % list of all connected nodes to this ocean boundary node
    nodes=find(bins==ID);
    % add these to the OpenOcean list
    OpenOceanNodes(nodes) = 1;
    % also remove these from the list of boundary nodes to save time where
    % one ice shelf has multiple nodes on the ocean boundary
    Nnum(nodes) = 0;
    
end

%% Now identify ice shelf nodes
% temporary mask: 0 for nodes outside Ua domain or with b>B+1, 1 inside remainder of the boundary
Mask_tmp = 0*lonCMIT+1;
Mask_tmp([I_outsideUa;find(b_forMITgcm>B_forMITgcm+1)]) = 0;

% for graph, only retain elements that are not fully grounded
Mask_ele = sum(Mask_tmp(MITgcm_connectivity),2);
% conservative approach near ice front:
Iele_fullyfloatingnearIF = intersect(find(Mask_ele==0),Iele_crossingicefront);
% include all not-fully grounded elements elsewhere
Iele_notfullygrounded = setdiff(find(Mask_ele<3),Iele_crossingicefront);
TRI = MITgcm_connectivity([Iele_fullyfloatingnearIF(:);Iele_notfullygrounded(:)],:);

% create undirected graph
G=graph(TRI,TRI(:,[2 3 1]));
% calculate the connected components of the graph
bins=conncomp(G) ;

% initialise arrays
Nnum = zeros(nx*ny,1);
OpenOceanIceShelfNodes = zeros(nx*ny,1);
Nnum(OceanBoundaryNodes) = 1;

% loop through ocean boundary nodes until each one has been checked for
% connected floating nodes, once this is done for all boundary nodes the
% only floating nodes left should be lakes
while sum(Nnum)>0
    
    NodeSeed = find(Nnum,1,'first');
    ID=bins(NodeSeed) ;
    % list of all connected nodes to this ocean boundary node
    nodes=find(bins==ID);
    % add these to the OpenOcean list
    OpenOceanIceShelfNodes(nodes) = 1;
    % also remove these from the list of boundary nodes to save time where
    % one ice shelf has multiple nodes on the ocean boundary
    Nnum(nodes) = 0;
    
end

Mask_tmp = 0*lonCMIT + 2;
Mask_tmp(OpenOceanNodes==1) = 0;
IceShelfNodes = OpenOceanIceShelfNodes; IceShelfNodes(OpenOceanNodes==1)=0;
Mask_tmp(IceShelfNodes==1) = 1;

%% Now do some cleaning up:
%1. check that for ice shelf elements, at least 1 corner is afloat
I_IS = find(Mask_tmp==1);
lonLMIT = lonCMIT(I_IS)-dlon/2;
lonRMIT = lonCMIT(I_IS)+dlon/2;
latLMIT = latCMIT(I_IS)-dlat/2;
latUMIT = latCMIT(I_IS)+dlat/2;
tic
I_err = [];
for ii=1:numel(I_IS)
    [xC,yC] = ll2psxy([latLMIT(ii) latUMIT(ii) latUMIT(ii) latLMIT(ii)],...
        [lonLMIT(ii) lonLMIT(ii) lonRMIT(ii) lonRMIT(ii)],-71,0);
    if ~any(FB(xC(:),yC(:))-Fb(xC(:),yC(:))<0)
        I_err = [I_err; I_IS(ii)];
    end
end
toc
Mask_tmp(I_err) = 2; %set to grounded
B_forMITgcm(I_err) = b_forMITgcm(I_err);

Mask = Mask_tmp;

%% In UaMITgcm mask convention is different:
% 2 - open ocean
% 1 - ice shelf
% 0 - grounded
Mask(Mask==2)=3;
Mask(Mask==0)=2;
Mask(Mask==3)=0;
mask_forMITgcm = Mask(:);

%% finally, we perform a few consistency checks:
% -> b_forMITgcm = 0 for open ocean
% -> B_forMITgcm = b_forMITgcm whereever the mask indicates that the
% ice is grounded
% -> B_forMITgcm < b_forMITgcm whereever the mask indicates that the
% ice is afloat. If this is not the case, then the bed is carved
% away
% -> b_forMITgcm<0 wherever the mask indicates that the ice is afloat.
%    If this is not the case, the mask is edited to say the ice is
%    grounded, and the draft is set to the bathymetry.

disp('Running consistency checks');

disp(['There are ',num2str(nnz((mask_forMITgcm==2).*(b_forMITgcm~=0))),' open-ocean points with nonzero ice shelf draft']);
b_forMITgcm(mask_forMITgcm==2) = 0;

disp(['There are ',num2str(nnz((mask_forMITgcm==0).*(B_forMITgcm~=b_forMITgcm))),' grounded points where ice draft does not equal bedrock depth']);
B_forMITgcm(mask_forMITgcm==0) = b_forMITgcm(mask_forMITgcm==0);

disp(['There are ',num2str(nnz((mask_forMITgcm==1).*(B_forMITgcm>=b_forMITgcm))),' ice shelf points with negative water column thickness']);
Ierr = find((mask_forMITgcm==1).*(B_forMITgcm>=b_forMITgcm));
B_forMITgcm(Ierr) = b_forMITgcm(Ierr)-1;

disp(['There are ',num2str(nnz((mask_forMITgcm==1).*(b_forMITgcm>0))),' ice shelf points with positive draft']);
Ipos = find((mask_forMITgcm==1).*(b_forMITgcm>0));
b_forMITgcm(Ipos) = B_forMITgcm(Ipos);
mask_forMITgcm(Ipos) = 0;

%% make sure that output fields are 2D
mask_forMITgcm = reshape(mask_forMITgcm,nx,ny);
B_forMITgcm = reshape(B_forMITgcm,nx,ny);
b_forMITgcm = reshape(b_forMITgcm,nx,ny);

figure; hold on;

pcolor(lonCMIT-dlon/2,latCMIT-dlat/2,mask_forMITgcm);
triplot(TRI_OO,lonCMIT,latCMIT,'-g')
[UaBlat,UaBlon]=psxy2ll(MUA.Boundary.x,MUA.Boundary.y,-71,0);
plot(UaBlon,UaBlat,'-xk');
CtrlVar.PlotGLs = 0;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF);
[latGL,lonGL] = psxy2ll(xGL,yGL,-71,0);
plot(lonGL,latGL,'-xk');

%% now save B, b and mask
save([UserVar.UaMITgcm.UaOutputDirectory,'/',UserVar.UaMITgcm.UaDraftFileName],'B_forMITgcm','b_forMITgcm','mask_forMITgcm','-v6');
end