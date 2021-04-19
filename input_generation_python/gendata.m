function gendata

%% generates required input data for Amundsen Sea run with Bedmachine bathymetry,
%% OBCS restoring at domain boundaries and RBCS restoring at surface

%% 1. Define some constants
rhoConst = 1024;
hFacMin = 0.05;
hFacMinDr = 0;
eosType = 'MDJWF';

%% 2. Generate MITgcm grid
% Grid dimensions. These are all things set in the input/data namelist.
nx = 240;    % first part of delX
dx = 1000;   % second part of delX
ny = 480;    % first part of delY
dy = 1000;   % second part of delY
nz = [80, 10];     % first part of delZ
dz = [20, 40];   % second part of delZ

% Build vertical grid
delz=[];
for ii=1:length(nz)
    delz = [delz dz(ii)*ones(1,nz(ii))];
end
zgp1 = [0,cumsum(delz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);

% Build horizontal grid
XC = linspace(-1.7e6+dx/2,-1.7e6+(nx-1/2)*dx,dx);
YC = linspace(-7e5+dy/2,-7e5+(ny-1/2)*dy,dy);
[Xm,Ym]=ndgrid(XC,YC);
       

%% 3. Use Bedmachine interpolants to generate geometry files
BMdir = '/Volumes/mainJDeRydt/Antarctic_datasets/BedMachine_Antarctica';
if exist('BedMachineGriddedInterpolants_2020-07-15.mat')~=2
    BMfile = [BMdir,'/BedMachineAntarctica-2020-07-15.nc'];
    [Fs,Fb,FB,Frho,~]=BedMachineToUaGriddedInterpolants(BMfile,1);
    save([BMdir,'/BedMachineGriddedInterpolants_2020-07-15'],'Fs','Fb','FB','Frho','-v7.3');
else
    if isempty(Fs)
        load([BMdir,'/BedMachineGriddedInterpolants_2020-07-15.mat']);
    end
end
bathy = FB(Xm,Ym);


% 3. Generate initial T&S fields


% 4. Generate OBCS restoring fields (northern and western boundaries)


% 5. Generate RBCS restoring fields (surface)

% # Some additional stuff about the forcing
% obcs_forcing_data = 'Kimura' # either 'Kimura' or 'Holland'
% constant_forcing = False# False # if set to True, the forcing from options.startDate+options.spinup_time will be taken
% 
% # read information about startDates, spinup time and simulation time from the options
% options = Options()
% ini_year = int(options.startDate[:4]) # this is the start year and should take the spinup into account
% ini_month = int(options.startDate[4:6]) # this is the start month and should take the spinup into account
% spinup = int(options.spinup_time) # in months
% totaltime = int(options.total_time) # in months
% 
% # generate array of months/years for OBCS forcing
% class OBCSForcingArray:
% 
%     def __init__ (self):
%         # first, initialize variables
%         self.nt = totaltime
%         self.years,self.months = np.zeros(totaltime), np.zeros(totaltime)
% 
%         # assign years and months for forcing
%         if constant_forcing:
%             print 'Constant OBCS forcing turned ON'
%             out = raw_input('You have chosen constant OBCS forcing. Enter the date code (eg 199201):').strip()
%             # make sure input is a valid date
%             valid_date = len(out)==6
%             try:
%                 int(valid_date)
%             except(ValueError):
%                 valid_date = False
%             if not valid_date:
%                 print 'Error: invalid date code ' + out
%                 sys.exit()
%             # assign input to array
%             self.years = self.years + int(out[:4])
%             self.months = self.months + int(out[4:6])
%         else:
%             print 'Time-varying OBCS forcing turned ON'
%             self.years = self.years + ini_year + np.floor(np.arange(totaltime)/12)
%             self.months = self.months + np.mod(np.arange(totaltime),12) + 1
% 
