
# Create the initial topography, initial conditions,
# and boundary conditions for MITgcm.
######################################################

import numpy as np
from scipy.io import loadmat
from scipy import interpolate
import sys
import xarray as xr
import matplotlib.pyplot as plt
from pyproj import Proj
import pandas as pd

# Get mitgcm_python in the path and load some functions
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/tools/')
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/coupling/')
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/MITgcm_67k/utils/python/MITgcmutils')
from mitgcm_python.file_io import write_binary
from mitgcm_python.utils import z_to_xyz, calc_hfac
from mitgcm_python.make_domain import do_digging, do_zapping
from mitgcm_python.ics_obcs import calc_load_anomaly
from set_parameters import Options

# Global parameters
# These are all things set in the input/data namelist.
nx = 240    # first part of delX
dx = 1000   # second part of delX
ny = 480    # first part of delY
dy = 1000   # second part of delY
nz = [80, 10]   # first part of delZ
dz = [20, 40]   # second part of delZ
eosType = 'MDJWF'
rhoConst = 1024.
hFacMin = 0.05
hFacMinDr = 0.

# Some additional stuff about the forcing data
forcing_data = 'Holland' # either 'Kimura' or 'Holland'
constant_forcing = False # False # if set to True, the forcing from options.startDate+options.spinup_time will be taken

# read information about startDates, spinup time and simulation time from the options
options = Options()
ini_year = int(options.startDate[:4]) # this is the start year and should take the spinup into account
ini_month = int(options.startDate[4:6]) # this is the start month and should take the spinup into account
totaltime = int(options.total_time) # in months

# ============================================================================================
# gather info on forcing data
class ForcingInfo:

    def __init__ (self):
        # first, initialize variables
        self.nt = totaltime
        self.years,self.months = np.zeros(totaltime), np.zeros(totaltime)

        # assign years and months for forcing
        if constant_forcing:
            print 'Constant OBCS forcing turned ON - write some code here to work out the averaging'
            sys.exit()
        else:
            print 'Time-varying OBCS forcing turned ON'
            self.years = self.years + ini_year + np.floor(np.arange(totaltime)/12)
            self.months = self.months + np.mod(np.arange(totaltime),12) + 1
        # assign forcing data
        if forcing_data == 'Kimura':
            print 'Using Kimura data for restoring and initial conditions. I do not have the surface data ' \
                  'so quit here'
            sys.exit()
            self.Tfile = ''
            self.Sfile = ''
            self.BC = loadmat('/Volumes/mainJDeRydt/UaMITgcm_v2/MIT_InputData/Kimura_OceanBC.mat')
        elif forcing_data == 'Holland':
            print 'Using Holland data for restoring and initial conditions. OBCS requires a mat file ' \
                  ' that contains T, S, U and V at the boundaries. Matlab processing scripts for OBCS restoring need to ' \
                  'be converted to Python at some point...'
            self.Tfile = '/Volumes/mainJDeRydt/Antarctic_datasets/PIG_Thwaites_DotsonCrosson/MITgcm_AS/PHolland_2019/PAS_851/run/stateTheta.nc'
            self.Sfile = '/Volumes/mainJDeRydt/Antarctic_datasets/PIG_Thwaites_DotsonCrosson/MITgcm_AS/PHolland_2019/PAS_851/run/stateSalt.nc'
            self.BC = loadmat('/Volumes/mainJDeRydt/UaMITgcm_v2/MIT_InputData/Holland_OceanBC_PAS_851.mat')
        else: 
            print 'Error: input data for restoring and initial conditions not found'

        BCyears = np.where(self.BC['year'][-1,:] == ini_year)
        BCmonths = np.where(self.BC['month'][-1,:] == ini_month)
        startIndex = np.int(np.intersect1d(BCyears,BCmonths))
        self.startIndex = startIndex
        self.BC['Theta'] = self.BC['Theta'][:,:,startIndex:startIndex+totaltime]
        self.BC['Salt'] = self.BC['Salt'][:,:,startIndex:startIndex+totaltime]
        self.BC['Ups'] = self.BC['Ups'][:,:,startIndex:startIndex+totaltime]
        self.BC['Vps'] = self.BC['Vps'][:,:,startIndex:startIndex+totaltime]
        self.BC['year'] = self.BC['year'][:,startIndex:startIndex+totaltime]
        self.BC['month'] = self.BC['month'][:,startIndex:startIndex+totaltime]

# ============================================================================================
# BasicGrid object to hold some information about the grid - just the variables we need to create all the initial conditions, with the same conventions as the mitgcm_python Grid object where needed. This way we can call calc_load_anomaly without needing a full Grid object.
class BasicGrid:

    def __init__ (self):
        # Build vertical grid
        self.z_edges = [0]
        for x in xrange(0,len(dz)):
            self.newedges = -np.arange(dz[x],(nz[x]+1)*dz[x],dz[x])
            self.z_edges = np.concatenate((self.z_edges,self.z_edges[-1]+np.array(self.newedges)),axis=None)
        self.z = 0.5*(self.z_edges[:-1] + self.z_edges[1:])
        self.dz = -self.z_edges[1:] + self.z_edges[:-1]
        # Build horizontal grid
        self.x = np.arange(-1.7e6+dx/2,-1.7e6+(nx-1/2)*dx,dx)
        self.y = np.arange(-7e5+dy/2,-7e5+(ny-1/2)*dy,dy)
        self.x2d, self.y2d = np.meshgrid(self.x, self.y)
        # Save grid dimensions
        self.nx = nx
        self.ny = ny
        self.nz = np.sum(nz)
        # convert regional grid to lat and lon for interpolation
        p = Proj('+init=EPSG:3031')
        MITx2d, MITy2d = np.meshgrid(self.x, self.y)
        self.lon2d, self.lat2d = p(MITx2d, MITy2d, inverse=True)
        self.lon2d = self.lon2d + 360

    # Calculate hFacC given the bathymetry and ice shelf draft.
    # Save to the object.
    def save_hfac (self, bathy, draft):
        self.hfac = calc_hfac(bathy, draft, self.z_edges, hFacMin=hFacMin, hFacMinDr=hFacMinDr)

    # Compatibility function with Grid.
    def get_hfac (self, gtype='t'):
        if gtype != 't':
            print 'Error (BasicGrid.get_hfac): hfac only exists on tracer grid'
            sys.exit()
        return self.hfac
    
# end BasicGrid object

# ============================================================================================
# Calculate the topography and write to binary files.
def make_topo (grid, ua_topo_file, bathy_file, draft_file, prec=64, dig_option='none'):

    # Read Bedmachine data
    BM = xr.open_dataset(
        '/Volumes/mainJDeRydt/Antarctic_datasets/BedMachine_Antarctica/BedMachineAntarctica-2020-07-15.nc')
    X2d, Y2d = np.meshgrid(BM.x, BM.y)

    # Reduce the size of dataset
    xmin, xmax, ymin, ymax = np.min(grid.x), np.max(grid.x), np.min(grid.y), np.max(grid.y)
    xmin = xmin - 1e4  # take a bit more for interpolation
    xmax = xmax + 1e4
    ymin = ymin - 1e4
    ymax = ymax + 1e4

    cond2d = ((X2d >= xmin) & (X2d <= xmax) & (Y2d >= ymin) & (Y2d <= ymax))
    X1d = X2d[cond2d]
    Y1d = Y2d[cond2d]

    draft = BM.surface.values - BM.thickness.values
    bathy = BM.bed.values
    mask = BM.mask.values # open ocean 0, rock outcrops 1, grounded ice 2, ice shelf 3
    draft1d = draft[cond2d]
    bathy1d = bathy[cond2d]
    mask1d = mask[cond2d]

    # Mask out grounded ice and outcrops in bathy
    bathy1d[(mask1d == 1) | (mask1d == 2)] = np.nan
    # Mask out open ocean, grounded ice and outcrops in draft
    draft1d[(mask1d == 0) | (mask1d == 1) | (mask1d == 2)] = np.nan

    # Interpolate onto MITgcm grid and restore NaN values to zero
    MITx1d, MITy1d = np.reshape(grid.x2d, grid.nx * grid.ny), np.reshape(grid.y2d, grid.nx * grid.ny)
    draft = interpolate.griddata((X1d, Y1d), draft1d, (MITx1d, MITy1d), method='linear',
                                 fill_value=np.nan)
    draft[np.isnan(draft)] = 0
    draft = np.reshape(draft, (grid.ny, grid.nx))
    bathy =  interpolate.griddata((X1d, Y1d), bathy1d, (MITx1d, MITy1d), method='linear',
                                 fill_value=np.nan)
    bathy[np.isnan(bathy)] = 0
    bathy = np.reshape(bathy, (grid.ny, grid.nx))

    #N=1
    #plt.pcolor(MITx2d[0:-1:N,0:-1:N],MITy2d[0:-1:N,0:-1:N],draft2d[0:-1:N,0:-1:N])
    #plt.colorbar()
    #plt.show()

    if dig_option == 'none':
        print 'Not doing digging as per user request'
    elif dig_option == 'bathy':
        print 'Digging bathymetry which is too shallow'
        bathy = do_digging(bathy, draft, grid.dz, grid.z_edges, hFacMin=hFacMin, hFacMinDr=hFacMinDr, dig_option='bathy')
    elif dig_option == 'draft':
        print 'Digging ice shelf drafts which are too deep'
        draft = do_digging(bathy, draft, grid.dz, grid.z_edges, hFacMin=hFacMin, hFacMinDr=hFacMinDr, dig_option='draft')

    print 'Zapping ice shelf drafts which are too thin'
    draft = do_zapping(draft, draft!=0, grid.dz, grid.z_edges, hFacMinDr=hFacMinDr)[0]

    # Calculate hFacC and save to the grid for later
    grid.save_hfac(bathy, draft)

    # Write to file
    write_binary(bathy, bathy_file, prec=prec)
    write_binary(draft, draft_file, prec=prec)

# ============================================================================================
# Returns temperature and salinity profiles, varying with depth, to be used for initial and boundary conditions.
def ts_profile(x,y,z,obcs):

    sizetz = (obcs.nt,np.sum(nz))
    t_profile, s_profile, u_profile, v_profile = np.zeros(sizetz), np.zeros(sizetz), np.zeros(sizetz), np.zeros(sizetz)

    L = np.sqrt((x-obcs.BC['x'][:,0])**2+(y-obcs.BC['y'][:,0])**2)
    IL = np.nanargmin(L)
    
    for i in range(0,obcs.nt):    
        findtime = np.in1d(obcs.BC['year'],obcs.years[i]) & np.in1d(obcs.BC['month'],obcs.months[i])
        Itime = np.where(findtime)
        Itime = Itime[0][0]
        t_profile[i,:] = np.interp(-z,-obcs.BC['depth'][:,0],obcs.BC['Theta'][IL,:,Itime])
        s_profile[i,:] = np.interp(-z,-obcs.BC['depth'][:,0],obcs.BC['Salt'][IL,:,Itime])
        u_profile[i,:] = np.interp(-z,-obcs.BC['depth'][:,0],obcs.BC['Ups'][IL,:,Itime])
        v_profile[i,:] = np.interp(-z,-obcs.BC['depth'][:,0],obcs.BC['Vps'][IL,:,Itime])

    return t_profile, s_profile, u_profile, v_profile

# ============================================================================================
# Creates OBCS for the southern/western boundary
def make_obcs (grid, forcing, obcs_temp_file_S, obcs_salt_file_S, obcs_uvel_file_S, obcs_vvel_file_S, obcs_temp_file_W, obcs_salt_file_W, obcs_uvel_file_W, obcs_vvel_file_W, prec):
    
    sizetzx = (forcing.nt,np.sum(nz),nx)
    sizetzy = (forcing.nt,np.sum(nz),ny)

    ## Southern boundary
    OBS_t, OBS_s, OBS_u, OBS_v = np.zeros(sizetzx), np.zeros(sizetzx), np.zeros(sizetzx), np.zeros(sizetzx)
    
    for i in xrange(0,nx):
        x = grid.x[i]
        y = grid.y[0]-dy/2
        t_profile, s_profile, u_profile, v_profile = ts_profile(x,y,grid.z,forcing)
        OBS_t[:,:,i] = t_profile
        OBS_s[:,:,i] = s_profile
        OBS_u[:,:,i] = u_profile
        OBS_v[:,:,i] = v_profile

    # Write the files
    # No need to mask out the land because MITgcm will do that for us
    write_binary(OBS_t, obcs_temp_file_S, prec)
    write_binary(OBS_s, obcs_salt_file_S, prec)
    write_binary(OBS_u, obcs_uvel_file_S, prec)
    write_binary(OBS_v, obcs_vvel_file_S, prec)

    # Remove variables from workspace
    OBS_t, OBS_s, OBS_u, OBS_v = None, None, None, None

    ## Western boundary
    OBW_t, OBW_s, OBW_u, OBW_v = np.zeros(sizetzy), np.zeros(sizetzy), np.zeros(sizetzy), np.zeros(sizetzy)
    
    for i in xrange(0,ny):
        x = grid.x[0]-dx/2
        y = grid.y[i]
        t_profile, s_profile, u_profile, v_profile = ts_profile(x,y,grid.z,forcing)
        OBW_t[:,:,i] = t_profile
        OBW_s[:,:,i] = s_profile
        OBW_u[:,:,i] = u_profile
        OBW_v[:,:,i] = v_profile

    # Write the files
    write_binary(OBW_t, obcs_temp_file_W, prec)
    write_binary(OBW_s, obcs_salt_file_W, prec)
    write_binary(OBW_u, obcs_uvel_file_W, prec)
    write_binary(OBW_v, obcs_vvel_file_W, prec)

    # Remove variable from workspace
    OBW_t, OBW_s, OBW_u, OBW_v = None, None, None, None

# ============================================================================================
# Creates RBCS files for surface restoring
def make_rbcs(grid, forcing, rbcs_temp_file, rbcs_salt_file, rbcs_tempmask_file, rbcs_saltmask_file, prec):

    sizetyx = (totaltime, grid.ny, grid.nx)
    sizezyx = (grid.nz, grid.ny, grid.nx)

    # define forcing data
    T = xr.open_dataset(forcing.Tfile)
    S = xr.open_dataset(forcing.Sfile)
    lon, lat = T.LONGITUDE, T.LATITUDE
    lon2d, lat2d = np.meshgrid(lon, lat)
    lon1d = lon2d.ravel()
    lat1d = lat2d.ravel()

    # regional grid in lat and lon for interpolation
    MITlon2d, MITlat2d = grid.lon2d, grid.lat2d

    # Now extract T & S arrays at the surface from the forcing data.
    # New T & S arrays have dimensions (nt, ny, nx),
    # The time slices correspond to restoring conditions at the end of each calendar month, averaged over
    # the preceding month (will be interpolated onto regular interval of 30 days later on - this is required by RBCS)
    # The horizontal slices are regridded to the local MITgcm grid
    k=0 # a dummy index
    RBsurf_t_full = np.zeros(sizetyx) # define T & S arrays with correct dimensions
    RBsurf_s_full = np.zeros(sizetyx)
    RBsurf_time = np.empty(totaltime-1,dtype='datetime64[ns]') # define time array
    for t in range(forcing.startIndex, forcing.startIndex+totaltime-1):
        # slice full T and S arrays to surface layer and time t
        Tsurf_slice = T.isel(DEPTH=slice(1),TIME=t)
        Ssurf_slice = S.isel(DEPTH=slice(1),TIME=t)
        #Tsurf = Tsurf_slice.THETA.values.ravel()
        #Ssurf = Ssurf_slice.SALT.values.ravel()
        # eliminate zeros as these will lead to unrealistic T & S values during interpolation
        #Tsurf = Tsurf[Ssurf != 0]
        #x = lon1d[Ssurf != 0]
        #y = lat1d[Ssurf != 0]
        #Ssurf = Ssurf[Ssurf != 0]
        # interpolate using linear barycentric interpolation where triangles exist,
        # and allow extrapolation with nearest values where triangles do not exist. The latter only occurs in very small
        # areas close to the ice front
        #t_xx = interpolate.griddata((x, y), Tsurf, (MITlon2d, MITlat2d), method='linear',fill_value=np.nan)
        #t_ss = interpolate.griddata((x, y), Tsurf, (MITlon2d, MITlat2d), method='nearest')
        #t_xx[(np.isnan(t_xx))] = t_ss[(np.isnan(t_xx))]
        #s_xx = interpolate.griddata((x, y), Ssurf, (MITlon2d, MITlat2d), method='linear',fill_value=np.nan)
        #s_ss = interpolate.griddata((x, y), Ssurf, (MITlon2d, MITlat2d), method='nearest')
        #s_xx[(np.isnan(s_xx))] = s_ss[(np.isnan(s_xx))]
        # assign regridded timeslice to full array
        #RBsurf_t_full[k, :, :] = t_xx
        #RBsurf_s_full[k, :, :] = s_xx
        RBsurf_time[k] = Tsurf_slice.TIME.values
        # print some info and step dummy index
        print 'Done ',k+1,' out of ',totaltime
        k=k+1
        # N = 1
        # plt.pcolor(MITx2d[0:-1:N, 0:-1:N], MITy2d[0:-1:N, 0:-1:N], RBsurf_s[k, 0:-1:N, 0:-1:N])
        # plt.colorbar()
        # plt.show()
        # sys.exit()

    # We use RBCS in /Non-cyclic data, multiple files/ mode
    # RBCS has a constant forcing period so we interpolate RBsurf_t and RBsurf_s
    # linearly to fit equally spaced time intervals
    # original timeseries and corresponding indices
    torig = pd.to_datetime(RBsurf_time)
    df = pd.DataFrame(np.arange(0,totaltime-1),index=torig,columns=list('I'))
    df = df.resample('D').mean().interpolate('linear')
    teq = df.resample('30D').nearest()

    Ind = (np.floor(teq['I']), np.ceil(teq['I']))
    Ind = np.asarray(Ind)
    Fracs = (Ind[1, :] - teq['I'], teq['I'] - Ind[0, :])
    Fracs = np.asarray(Fracs)
    Fracs[(Fracs == 0) & (Ind[0, :] == Ind[1, :])] = 0.5
    Ind = Ind.astype('int')

    # Interpolate to regular intervals
    for t in range(0,np.size(Ind,1)):
        RBsurf_t = np.zeros(sizezyx)
        RBsurf_s = np.zeros(sizezyx)
        RBsurf_t_reg = Fracs[0, t] * RBsurf_t_full[Ind[0, t], :, :] + Fracs[1, t] * RBsurf_t_full[Ind[1, t], :, :]
        RBsurf_s_reg = Fracs[0, t] * RBsurf_s_full[Ind[0, t], :, :] + Fracs[1, t] * RBsurf_s_full[Ind[1, t], :, :]
        # generate 3D arrays for T, S at each time interval
        RBsurf_t[0, :, :] = RBsurf_t_reg
        RBsurf_s[0, :, :] = RBsurf_s_reg
        # Write binary files for T, S and Mask at each time interval.
        # No need to take care of mask as MITgcm will do this.
        write_binary(RBsurf_t, rbcs_temp_file+str(t).zfill(10)+'.data', prec)
        write_binary(RBsurf_s, rbcs_salt_file+str(t).zfill(10)+'.data', prec)

    # Write binary Mask files
    Mask = np.zeros(sizezyx)
    Mask[0,:,:]=1
    write_binary(Mask, rbcs_tempmask_file, prec)
    write_binary(Mask, rbcs_saltmask_file, prec)

    # Remove variables from workspace
    RBsurf_t_full, RBsurf_s_full, RBsurf_t, RBsurf_s = None, None, None, None

# ============================================================================================
# Creates initial T & S fields and calculates pressure loading
def make_ics(grid, forcing, ini_temp_file, ini_salt_file, pload_file, prec):

    # T & S arrays for ICs have dimensions (nz, ny, nx)
    sizezyx = (grid.nz, grid.ny, grid.nx)
    ics_T = np.zeros(sizezyx)  # define T & S arrays with correct dimensions
    ics_S = np.zeros(sizezyx)
    ics_T_vertinterp = np.zeros((grid.nz, nx * ny))
    ics_S_vertinterp = np.zeros((grid.nz, nx * ny))

    # define forcing data
    T = xr.open_dataset(forcing.Tfile)
    S = xr.open_dataset(forcing.Sfile)
    lon, lat, z = T.LONGITUDE, T.LATITUDE, T.DEPTH
    nx, ny, nz = lon.shape[0], lat.shape[0], z.shape[0]
    lon2d, lat2d = np.meshgrid(lon, lat)
    lon1d = lon2d.ravel()
    lat1d = lat2d.ravel()

    # Initial T&S conditions are from forcing dataset at time forcing.startIndex
    # slice full T and S arrays to correct time
    ics_Tslice = T.isel(TIME=forcing.startIndex).THETA.values
    ics_Sslice = S.isel(TIME=forcing.startIndex).SALT.values

    # linearly interpolate vertical levels first
    FI = interpolate.interp1d(z,np.arange(0,nz),kind='linear')
    I = FI(-grid.z)
    Ind = (np.floor(I), np.ceil(I))
    Ind = np.asarray(Ind)
    Fracs = (Ind[1, :] - I, I - Ind[0, :])
    Fracs = np.asarray(Fracs)
    Fracs[(Fracs == 0) & (Ind[0, :] == Ind[1, :])] = 0.5
    Ind = Ind.astype('int')
    ics_Tslice = np.reshape(ics_Tslice, (nz, nx * ny))
    ics_Sslice = np.reshape(ics_Sslice, (nz, nx * ny))
    for t in range(0, np.size(Ind, 1)):
        ics_T_vertinterp[t, :] = Fracs[0, t] * ics_Tslice[Ind[0, t], :] + Fracs[1, t] * ics_Tslice[Ind[1, t], :]
        ics_S_vertinterp[t, :] = Fracs[0, t] * ics_Sslice[Ind[0, t], :] + Fracs[1, t] * ics_Sslice[Ind[1, t], :]

    # Now use linear barycentric interpolation for horizontal slices, allowing for
    # extrapolation with nearest values where triangles do not exist
    # First eliminate zeros as these will lead to unrealistic T & S values during interpolation
    for k in range(0,3):#grid.nz):
        Tslice = ics_T_vertinterp[k,:]
        Sslice = ics_S_vertinterp[k,:]
        x, y = lon1d, lat1d
        Tslice = Tslice[Sslice != 0]
        x = x[Sslice != 0]
        y = y[Sslice != 0]
        Sslice = Sslice[Sslice != 0]
        t_xx = interpolate.griddata((x, y), Tslice, (grid.lon2d, grid.lat2d), method='linear',fill_value=np.nan)
        t_ss = interpolate.griddata((x, y), Tslice, (grid.lon2d, grid.lat2d), method='nearest')
        t_xx[(np.isnan(t_xx))] = t_ss[(np.isnan(t_xx))]
        s_xx = interpolate.griddata((x, y), Sslice, (grid.lon2d, grid.lat2d), method='linear',fill_value=np.nan)
        s_ss = interpolate.griddata((x, y), Sslice, (grid.lon2d, grid.lat2d), method='nearest')
        s_xx[(np.isnan(s_xx))] = s_ss[(np.isnan(s_xx))]
        # assign regridded horizontal slice to full array
        ics_T[k,:,:] = t_xx
        ics_S[k,:,:] = s_xx
        # print some info and step dummy index
        print 'Done ', k + 1, ' out of ', grid.nz
        #N = 1
        #plt.pcolor(grid.x2d[0:-1:N, 0:-1:N], grid.y2d[0:-1:N, 0:-1:N], ics_T[k, 0:-1:N, 0:-1:N])
        #plt.colorbar()
        #plt.show()
        #sys.exit()

    # Write the files. Don't worry about masking ice and bedrock cells, as MITgcm takes care of this
    write_binary(ics_T, ini_temp_file, prec=prec)
    write_binary(ics_S, ini_salt_file, prec=prec)

    # Calculate the pressure load anomaly
    calc_load_anomaly(grid, pload_file, option='precomputed', ini_temp_file=ini_temp_file, ini_salt_file=ini_salt_file,
                      eosType=eosType, rhoConst=rhoConst, prec=prec, check_grid=False)


############## USER INPUT HERE #########################
# Path to MITgcm input/ directory
input_dir = options.mit_case_dir+'input/'

print 'Building grid'
grid = BasicGrid()

print 'Reading obcs data'
forcing = ForcingInfo()

#print 'Creating topography'
make_topo(grid, './ua_custom/DataForMIT.mat', input_dir+'bathymetry.shice', input_dir+'shelfice_topo.bin', prec=64, dig_option='bathy')

#print 'Creating initial and boundary conditions'
#make_obcs(grid, forcing, input_dir+'OBSt.bin', input_dir+'OBSs.bin', input_dir+'OBSu.bin', input_dir+'OBSv.bin', input_dir+'OBWt.bin', input_dir+'OBWs.bin', input_dir+'OBWu.bin', input_dir+'OBWv.bin', prec=32)
#make_rbcs(grid, forcing, input_dir+'RBsurf_t', input_dir+'RBsurf_s', input_dir+'RBsurf_tmask', input_dir+'RBsurf_smask', prec=32)
make_ics(grid, forcing, input_dir+'T_ini.bin', input_dir+'S_ini.bin', input_dir+'pload.mdjwf', prec=32)
    
