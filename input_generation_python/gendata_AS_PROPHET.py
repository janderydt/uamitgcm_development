
# Create the initial topography, initial conditions,
# and boundary conditions for MITgcm.
######################################################

import numpy as np
from scipy.io import loadmat
from scipy import interpolate
import sys
import xarray as xr
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pyproj import Proj
import pandas as pd
from dateutil.relativedelta import relativedelta

# Get mitgcm_python in the path and load some local functions
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/tools/')
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/coupling/')
sys.path.append('/Volumes/mainJDeRydt/UaMITgcm_v2/UaMITgcm_source/MITgcm_67k/utils/python/MITgcmutils')
from mitgcm_python.file_io import write_binary
from mitgcm_python.utils import z_to_xyz, calc_hfac
from mitgcm_python.make_domain import do_digging, do_zapping
from mitgcm_python.ics_obcs import calc_load_anomaly
from set_parameters import Options
from input_generation_python import interp_functions
from convert_ll2psxy import uv_2psxy

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
time_interpolation = 'endofmonth' # time interpolation method for restoring files. Options are 'regularintervals', 'endofmonth'

# read information about startDates, spinup time and simulation time from the options
options = Options()
ini_year = str(options.startDate[:4]) # this is the start year and should take the spinup into account
ini_month = str(options.startDate[4:6]) # this is the start month and should take the spinup into account
starttime = pd.to_datetime(ini_year+'-'+ini_month+'-01', format = '%Y-%m-%d')
simulationtime_months = options.total_time
endtime = starttime + relativedelta(months=simulationtime_months)
simulationtime_seconds = (endtime-starttime).total_seconds()

# ============================================================================================
# gather info on forcing data
class ForcingInfo:

    def __init__ (self):

        if constant_forcing:
            print 'Constant OBCS forcing turned ON - write some code here to work out the averaging'
            sys.exit()
        else:
            print 'Time-varying OBCS forcing turned ON'

        # assign forcing data
        if forcing_data == 'Kimura':
            print 'Using Kimura data for restoring and initial conditions.'
            sys.exit('I do not have the surface data so quit here')
            self.Tfile = ''
            self.Sfile = ''
            self.Ufile = ''
            self.Vfile = ''
            #self.BC = loadmat('/Volumes/mainJDeRydt/UaMITgcm_v2/MIT_InputData/Kimura_OceanBC.mat')
        elif forcing_data == 'Holland':
            print 'Using Holland data for restoring and initial conditions.'
            forcingdir = '/Volumes/mainJDeRydt/Antarctic_datasets/PIG_Thwaites_DotsonCrosson/MITgcm_AS/PHolland_2019/PAS_851/run/'
            self.Tfile = forcingdir + 'stateTheta.nc'
            self.Sfile = forcingdir + 'stateSalt.nc'
            self.Ufile = forcingdir + 'stateUVEL.nc'
            self.Vfile = forcingdir + 'stateVVEL.nc'
            #self.BC = loadmat('/Volumes/mainJDeRydt/UaMITgcm_v2/MIT_InputData/Holland_OceanBC_PAS_851.mat')
        else: 
            print 'Error: input data for restoring and initial conditions not found'

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
        self.dx = dx
        self.y = np.arange(-7e5+dy/2,-7e5+(ny-1/2)*dy,dy)
        self.dy = dy
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
    draft = interpolate.griddata((Y1d, X1d), draft1d, (MITy1d, MITx1d), method='linear',
                                 fill_value=np.nan)
    draft[np.isnan(draft)] = 0
    draft = np.reshape(draft, (grid.ny, grid.nx))

    bathy =  interpolate.griddata((Y1d, X1d), bathy1d, (MITy1d, MITx1d), method='linear',
                                 fill_value=np.nan)
    bathy[np.isnan(bathy)] = 0
    bathy = np.reshape(bathy, (grid.ny, grid.nx))

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

    #N = 1
    #plt.pcolor(grid.x2d[0:-1:N, 0:-1:N], grid.y2d[0:-1:N, 0:-1:N], draft[0:-1:N, 0:-1:N]-bathy[0:-1:N, 0:-1:N])
    #plt.colorbar()
    #plt.show()

    # Calculate hFacC and save to the grid for later
    grid.save_hfac(bathy, draft)

    # Write to file
    write_binary(bathy, bathy_file, prec=prec)
    write_binary(draft, draft_file, prec=prec)

# ============================================================================================
# Creates OBCS for the southern/western boundary
def make_obcs (grid, forcing, obcs_temp_file_S, obcs_salt_file_S, obcs_uvel_file_S, obcs_vvel_file_S, obcs_temp_file_W, obcs_salt_file_W, obcs_uvel_file_W, obcs_vvel_file_W, spinup, prec):

    # OBCS package requires xz and yz slices at regular time intervals, but the forcing data only exists at
    # irregular intervals, i.e. at the end of each month. There are three possible solutions:
    # 1. We set OBCS.ForcingPeriod = 2592000 (30 days) and interpolate the forcing data to regular time intervals.
    # 2. We take the total runtime of the simulation and divide by the number of months
    # 3. We assume a Gregorian Calendar with a regular cycle of leap years, for which OBCS.ForcingPeriod can be
    # computed as 0.25*(366 + 3*365)*86400/12 = 2629800
    # The 2nd or 3rd approach are often identical, and the easiest/cheapest

    if (time_interpolation != 'endofmonth') & (time_interpolation != 'regularintervals'):
        sys.exit('Uknown time interpolation method')

    # define forcing data
    T = xr.open_dataset(forcing.Tfile)
    S = xr.open_dataset(forcing.Sfile)
    U = xr.open_dataset(forcing.Ufile)
    V = xr.open_dataset(forcing.Vfile)
    lon, lat, z, time = T.LONGITUDE, T.LATITUDE, T.DEPTH, T.TIME
    lon2d, lat2d = np.meshgrid(lon, lat)

    # define boundary coordinates of our regional MITgcm configuration
    x1dW, y1dW = np.meshgrid(grid.x[0]-grid.dx/2,grid.y) # 'western' boundary, i.e. constant x
    x1dS, y1dS = np.meshgrid(grid.x,grid.y[0]-grid.dy/2) # 'southern' boundary, i.e. constant y

    # convert boundary coordinates to lat and lon for interpolation
    p = Proj('+init=EPSG:3031')
    lon1dW, lat1dW = p(x1dW, y1dW, inverse=True)
    lon1dW = lon1dW + 360
    lon1dS, lat1dS = p(x1dS, y1dS, inverse=True)
    lon1dS = lon1dS + 360

    # reduce the size of the forcing dataset to save memory
    lonmin, lonmax = np.min(np.concatenate([lon1dW, lon1dS.T])), np.max(np.concatenate([lon1dW, lon1dS.T]))
    latmin, latmax = np.min(np.concatenate([lat1dW, lat1dS.T])), np.max(np.concatenate([lat1dW, lat1dS.T]))
    lonmin, latmin = lonmin, latmin - 0.1  # take a bit more for interpolation
    lonmax, latmax = lonmax, latmax + 0.1

    cond1d = ((lon >= lonmin) & (lon <= lonmax))
    Ilon = np.where(cond1d)[0]
    cond1d = ((lat >= latmin) & (lat <= latmax))
    Ilat = np.where(cond1d)[0]
    Iz = np.where((-z > grid.z[-1]-300))[0]
    Ilonmin, Ilatmin, Izmin = np.min(Ilon).astype('int'), np.min(Ilat).astype('int'), np.min(Iz).astype('int')
    Ilonmax, Ilatmax, Izmax = np.max(Ilon).astype('int'), np.max(Ilat).astype('int'), np.max(Iz).astype('int')

    # Now loop through time dimension of forcing data and extract T,S,U and V at boundary
    if spinup: # if spinup then only save the first time slice
        nt = 1
    else:
        nt = np.shape(time)[0]

    # First initialize some variables
    OBW_t = np.zeros((nt, grid.nz, np.shape(lon1dW)[0]))
    OBW_s = np.zeros((nt, grid.nz, np.shape(lon1dW)[0]))
    OBW_u = np.zeros((nt, grid.nz, np.shape(lon1dW)[0]))
    OBW_v = np.zeros((nt, grid.nz, np.shape(lon1dW)[0]))
    OBS_t = np.zeros((nt, grid.nz, np.shape(lon1dS)[1]))
    OBS_s = np.zeros((nt, grid.nz, np.shape(lon1dS)[1]))
    OBS_u = np.zeros((nt, grid.nz, np.shape(lon1dS)[1]))
    OBS_v = np.zeros((nt, grid.nz, np.shape(lon1dS)[1]))

    for t in range(0,nt):
        # Extract T, S, U & V arrays from the forcing data for reduced spatial domain, horizontal level and time
        T_slice = T.isel(LONGITUDE=range(Ilonmin, Ilonmax+1), LATITUDE=range(Ilatmin, Ilatmax+1), DEPTH=range(Izmin, Izmax+1), TIME=t)
        S_slice = S.isel(LONGITUDE=range(Ilonmin, Ilonmax+1), LATITUDE=range(Ilatmin, Ilatmax+1), DEPTH=range(Izmin, Izmax+1), TIME=t)
        U_slice = U.isel(LONGITUDE=range(Ilonmin, Ilonmax+1), LATITUDE=range(Ilatmin, Ilatmax+1), DEPTH=range(Izmin, Izmax+1), TIME=t)
        V_slice = V.isel(LONGITUDE=range(Ilonmin, Ilonmax+1), LATITUDE=range(Ilatmin, Ilatmax+1), DEPTH=range(Izmin, Izmax+1), TIME=t)

        lon1d = lon2d[Ilat[:,None], Ilon[None,:]].ravel()
        lat1d = lat2d[Ilat[:, None], Ilon[None, :]].ravel()
        T_slice = T_slice.THETA.values
        S_slice = S_slice.SALT.values
        U_slice = U_slice.UVEL.values
        V_slice = V_slice.VVEL.values

        OBW_t_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dW)[0]))
        OBW_s_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dW)[0]))
        OBW_u_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dW)[0]))
        OBW_v_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dW)[0]))
        OBS_t_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dS)[1]))
        OBS_s_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dS)[1]))
        OBS_u_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dS)[1]))
        OBS_v_temp = np.zeros((np.shape(Iz)[0], np.shape(lon1dS)[1]))

        # interpolate horizontal slices
        for iz in range(0,np.shape(T_slice)[0]):
            t_y = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dW, lat1dW, T_slice[iz,:,:].ravel())
            s_y = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dW, lat1dW, S_slice[iz,:,:].ravel())
            u_y = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dW, lat1dW, U_slice[iz,:,:].ravel())
            v_y = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dW, lat1dW, V_slice[iz,:,:].ravel())
            OBW_t_temp[iz, :] = t_y.ravel()
            OBW_s_temp[iz, :] = s_y.ravel()
            OBW_u_temp[iz, :] = u_y.ravel()
            OBW_v_temp[iz, :] = v_y.ravel()

            t_x = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dS, lat1dS, T_slice[iz, :, :].ravel())
            s_x = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dS, lat1dS, S_slice[iz, :, :].ravel())
            u_x = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dS, lat1dS, U_slice[iz, :, :].ravel())
            v_x = interp_functions.horizontal_interp_nonan(lon1d, lat1d, lon1dS, lat1dS, V_slice[iz, :, :].ravel())
            OBS_t_temp[iz, :] = t_x.ravel()
            OBS_s_temp[iz, :] = s_x.ravel()
            OBS_u_temp[iz, :] = u_x.ravel()
            OBS_v_temp[iz, :] = v_x.ravel()

        # interpolate vertical levels
        (Int,Fracs) = interp_functions.vertical_interp(z, -grid.z)

        fracs0W_2d = np.repeat(Fracs[0, :, np.newaxis], np.shape(OBW_t_temp)[1], axis=1)
        fracs1W_2d = np.repeat(Fracs[1, :, np.newaxis], np.shape(OBW_t_temp)[1], axis=1)
        fracs0S_2d = np.repeat(Fracs[0, :, np.newaxis], np.shape(OBS_t_temp)[1], axis=1)
        fracs1S_2d = np.repeat(Fracs[1, :, np.newaxis], np.shape(OBS_t_temp)[1], axis=1)

        OBW_t[t, :, :] = fracs0W_2d * OBW_t_temp[Int[0, :], :] + fracs1W_2d * OBW_t_temp[Int[1, :], :]
        OBW_s[t, :, :] = fracs0W_2d * OBW_s_temp[Int[0, :], :] + fracs1W_2d * OBW_s_temp[Int[1, :], :]
        OBW_u[t, :, :] = fracs0W_2d * OBW_u_temp[Int[0, :], :] + fracs1W_2d * OBW_u_temp[Int[1, :], :]
        OBW_v[t, :, :] = fracs0W_2d * OBW_v_temp[Int[0, :], :] + fracs1W_2d * OBW_v_temp[Int[1, :], :]
        OBS_t[t, :, :] = fracs0S_2d * OBS_t_temp[Int[0, :], :] + fracs1S_2d * OBS_t_temp[Int[1, :], :]
        OBS_s[t, :, :] = fracs0S_2d * OBS_s_temp[Int[0, :], :] + fracs1S_2d * OBS_s_temp[Int[1, :], :]
        OBS_u[t, :, :] = fracs0S_2d * OBS_u_temp[Int[0, :], :] + fracs1S_2d * OBS_u_temp[Int[1, :], :]
        OBS_v[t, :, :] = fracs0S_2d * OBS_v_temp[Int[0, :], :] + fracs1S_2d * OBS_v_temp[Int[1, :], :]

        # print some info
        print 'Done ',t+1,' out of ',nt,' time slices'

        #L, Z = np.meshgrid(lon1dW, grid.z)
        #plt.pcolor(Z, L, np.squeeze(OBW_t[0, :, :]))
        #plt.colorbar()
        #plt.show()
        #sys.exit()

    # Remove variables from workspace
    OBS_t_temp, OBS_s_temp, OBS_u_temp, OBS_v_temp = None, None, None, None
    fracs0W_2d, fracs0S_2d, fracs1S_2d, fracs1W_2d = None, None, None, None
    T_slice, S_slice, U_slice, V_slice = None, None, None, None

    ## WHAT TO DO WITH TIME INTERPOLATION???
    if ((nt != 1) & (time_interpolation == 'regularintervals')):
        # interpolate to regular time intervals of 30 days
        (Int, Fracs) = interp_functions.time_interp_reg(time.TIME.values, 30)

        # Western Boundary
        fracs0W_3d = np.repeat(Fracs[0, :, np.newaxis], np.shape(OBW_t)[1], axis=1)
        fracs0W_3d = np.repeat(fracs0W_3d[:, : ,np.newaxis], np.shape(OBW_t)[2], axis=2)
        fracs1W_3d = np.repeat(Fracs[1, :, np.newaxis], np.shape(OBW_t)[1], axis=1)
        fracs1W_3d = np.repeat(fracs1W_3d[:, :, np.newaxis], np.shape(OBW_t)[2], axis=2)

        OBW_t = fracs0W_3d * OBW_t[Int[0, :], :, :] + fracs1W_3d * OBW_t[Int[1, :], :, :]
        OBW_s = fracs0W_3d * OBW_s[Int[0, :], :, :] + fracs1W_3d * OBW_s[Int[1, :], :, :]
        OBW_u = fracs0W_3d * OBW_u[Int[0, :], :, :] + fracs1W_3d * OBW_u[Int[1, :], :, :]
        OBW_v = fracs0W_3d * OBW_v[Int[0, :], :, :] + fracs1W_3d * OBW_v[Int[1, :], :, :]

        # Southern Boundary
        fracs0S_3d = np.repeat(Fracs[0, :, np.newaxis], np.shape(OBS_t)[1], axis=1)
        fracs0S_3d = np.repeat(fracs0S_3d[:, :, np.newaxis], np.shape(OBS_t)[2], axis=2)
        fracs1S_3d = np.repeat(Fracs[1, :, np.newaxis], np.shape(OBS_t)[1], axis=1)
        fracs1S_3d = np.repeat(fracs1S_3d[:, :, np.newaxis], np.shape(OBS_t)[2], axis=2)

        OBS_t = fracs0S_3d * OBS_t[Int[0, :], :, :] + fracs1S_3d * OBS_t[Int[1, :], :, :]
        OBS_s = fracs0S_3d * OBS_s[Int[0, :], :, :] + fracs1S_3d * OBS_s[Int[1, :], :, :]
        OBS_u = fracs0S_3d * OBS_u[Int[0, :], :, :] + fracs1S_3d * OBS_u[Int[1, :], :, :]
        OBS_v = fracs0S_3d * OBS_v[Int[0, :], :, :] + fracs1S_3d * OBS_v[Int[1, :], :, :]

    elif ((nt != 1) & (time_interpolation == 'endofmonth')):
        # total simulation time in secs:
        print 'You have chosen to provide restoring files at the end of each month. Since OBCS needs a fixed' \
              'restoring period, you will need to set OBCS.ForcingPeriod to one of the values below (whichever' \
              'is most appropriate: \n'\
              '1. Assuming a regular leap year cycle: OBCS.ForcingPeriod = 2629800 \n'\
              '2. Divide total simulation time into equal intervals assuming 12 months a year: ', \
                simulationtime_seconds/simulationtime_months
        Int = np.zeros((1,1))

    else:
        Int = np.zeros((1, 1))

    # Convert u and v from lat lon to ps
    Ull = np.reshape(OBW_u,(np.shape(Int)[1]*grid.nz, np.shape(lon1dW)[0]))
    Vll = np.reshape(OBW_v,(np.shape(Int)[1]*grid.nz, np.shape(lon1dW)[0]))
    Ups, Vps = uv_2psxy(Ull,Vll,np.squeeze(lon1dW),np.squeeze(lat1dW))
    OBW_u = np.reshape(Ups, (np.shape(Int)[1], grid.nz, np.shape(lon1dW)[0]))
    OBW_v = np.reshape(Vps, (np.shape(Int)[1], grid.nz, np.shape(lon1dW)[0]))

    Ull = np.reshape(OBS_u, (np.shape(Int)[1] * grid.nz, np.shape(lon1dS)[1]))
    Vll = np.reshape(OBS_v, (np.shape(Int)[1] * grid.nz, np.shape(lon1dS)[1]))
    Ups, Vps = uv_2psxy(Ull,Vll,np.squeeze(lon1dS),np.squeeze(lat1dS))
    OBS_u = np.reshape(Ups, (np.shape(Int)[1], grid.nz, np.shape(lon1dS)[1]))
    OBS_v = np.reshape(Vps, (np.shape(Int)[1], grid.nz, np.shape(lon1dS)[1]))

    # Write the files
    # No need to mask out the land because MITgcm will do that for us
    write_binary(OBW_t, obcs_temp_file_W, prec)
    write_binary(OBW_s, obcs_salt_file_W, prec)
    write_binary(OBW_u, obcs_uvel_file_W, prec)
    write_binary(OBW_v, obcs_vvel_file_W, prec)
    write_binary(OBS_t, obcs_temp_file_S, prec)
    write_binary(OBS_s, obcs_salt_file_S, prec)
    write_binary(OBS_u, obcs_uvel_file_S, prec)
    write_binary(OBS_v, obcs_vvel_file_S, prec)

    # Remove variables from workspace
    OBW_t, OBW_s, OBW_u, OBW_v = None, None, None, None
    fracs0W_3d, fracs1W_3d = None, None
    OBS_t, OBS_s, OBS_u, OBS_v = None, None, None, None
    fracs0S_3d, fracs1S_3d = None, None

# ============================================================================================
# Creates RBCS files for surface restoring
def make_rbcs(grid, forcing, rbcs_temp_file, rbcs_salt_file, rbcs_tempmask_file, rbcs_saltmask_file, spinup, prec):

    # define forcing data
    T = xr.open_dataset(forcing.Tfile)
    S = xr.open_dataset(forcing.Sfile)
    lon, lat, time = T.LONGITUDE, T.LATITUDE, T.TIME
    lon2d, lat2d = np.meshgrid(lon, lat)
    lon1d = lon2d.ravel()
    lat1d = lat2d.ravel()

    # RBCS are from forcing dataset at time 'starttime'
    # slice full T and S arrays to correct time
    starttime_p1 = starttime + relativedelta(months=+1)  # shift starttime to first day of next month to correspond to time
    # conventions in forcing data
    forcingtime = pd.to_datetime(time.TIME.values, format='%Y-%m-%d')
    startIndex = np.where(starttime_p1 == forcingtime)[0]  # find index of starttime within time array of forcing data

    # Now loop through time dimension of forcing data and extract T,S,U and V at boundary
    if spinup:  # if spinup then only save the first time slice
        nt = 1
    else:
        nt = np.shape(time)[0]

    sizetyx = (nt, grid.ny, grid.nx)
    sizezyx = (grid.nz, grid.ny, grid.nx)

    # New T & S arrays have dimensions (nt, ny, nx),
    # The time slices in the forcing data correspond to T & S at the end of each calendar month, averaged over
    # the preceding month. RBCS requires restoring conditions at regular intervals (30 days in our case) so
    # we will need to do some interpolation from calendar dates to regular 30-day intervals.
    # Horizontal slices are regridded to the local MITgcm grid
    k=0 # a dummy index
    RBsurf_t_full = np.zeros(sizetyx) # define T & S arrays with correct dimensions
    RBsurf_s_full = np.zeros(sizetyx)
    RBsurf_time = np.empty(nt,dtype='datetime64[ns]') # define time array

    for t in range(startIndex, startIndex+nt):
        # slice full T and S arrays to surface layer and time t
        Tsurf_slice = T.isel(DEPTH=slice(1),TIME=t)
        Ssurf_slice = S.isel(DEPTH=slice(1),TIME=t)
        Tsurf = Tsurf_slice.THETA.values.ravel()
        Ssurf = Ssurf_slice.SALT.values.ravel()
        t_xx = interp_functions.horizontal_interp_nonan(lon1d,lat1d,grid.lon2d,grid.lat2d,Tsurf)
        s_xx = interp_functions.horizontal_interp_nonan(lon1d,lat1d,grid.lon2d,grid.lat2d,Ssurf)
        RBsurf_t_full[k, :, :] = t_xx
        RBsurf_s_full[k, :, :] = s_xx
        RBsurf_time[k] = Tsurf_slice.TIME.values
        # print some info and step dummy index
        print 'Done ',k+1,' out of ',nt
        k=k+1
        # N = 1
        #         # plt.pcolor(MITx2d[0:-1:N, 0:-1:N], MITy2d[0:-1:N, 0:-1:N], RBsurf_s[k, 0:-1:N, 0:-1:N])
        #         # plt.colorbar()
        #         # plt.show()
        #         # sys.exit()

    ## WHAT TO DO WITH TIME INTERPOLATION???
    if ((nt != 1 ) & (time_interpolation == 'regularintervals')):

        # We use RBCS in /Non-cyclic data, multiple files/ mode
        # RBCS has a constant forcing period so we interpolate RBsurf_t and RBsurf_s
        # linearly to fit equally spaced time intervals
        # original timeseries and corresponding indices
        # interpolate to regular time intervals of 30 days
        (Ind,Fracs) = interp_functions.time_interp_reg(RBsurf_time,30)

        # Interpolate to regular time intervals
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
            write_binary(RBsurf_t, rbcs_temp_file+'.'+str(t).zfill(10)+'.data', prec)
            write_binary(RBsurf_s, rbcs_salt_file+'.'+str(t).zfill(10)+'.data', prec)

    elif ((nt != 1) & (time_interpolation == 'endofmonth')):
        # total simulation time in secs:
        print 'You have chosen to provide restoring files at the end of each month. Since OBCS needs a fixed' \
              'restoring period, you will need to set OBCS.ForcingPeriod to one of the values below (whichever' \
              'is most appropriate: \n' \
              '1. Assuming a regular leap year cycle: OBCS.ForcingPeriod = 2629800 \n' \
              '2. Divide total simulation time into equal intervals assuming 12 months a year: ' + \
              simulationtime_seconds / simulationtime_months
        for t in range(0, nt):
            RBsurf_t = np.zeros(sizezyx)
            RBsurf_s = np.zeros(sizezyx)
            # generate 3D arrays for T, S at each time interval
            RBsurf_t[0, :, :] = RBsurf_t_full[t, :, :]
            RBsurf_s[0, :, :] = RBsurf_s_full[t, :, :]
            # Write binary files for T, S and Mask at each time interval.
            # No need to take care of mask as MITgcm will do this.
            write_binary(RBsurf_t, rbcs_temp_file + '.' + str(t).zfill(10) + '.data', prec)
            write_binary(RBsurf_s, rbcs_salt_file + '.' + str(t).zfill(10) + '.data', prec)

    else:
        RBsurf_t = np.zeros(sizezyx)
        RBsurf_s = np.zeros(sizezyx)
        # generate 3D arrays for T, S at each time interval
        RBsurf_t[0, :, :] = RBsurf_t_full[0, :, :]
        RBsurf_s[0, :, :] = RBsurf_s_full[0, :, :]
        # Write binary files for T, S and Mask at each time interval.
        # No need to take care of mask as MITgcm will do this.
        write_binary(RBsurf_t, rbcs_temp_file, prec)
        write_binary(RBsurf_s, rbcs_salt_file, prec)

    # Write binary Mask files
    Mask = np.zeros(sizezyx)
    Mask[0,:,:]=1

    write_binary(Mask, rbcs_tempmask_file, prec)
    write_binary(Mask, rbcs_saltmask_file, prec)

    # Remove variables from workspace
    RBsurf_t_full, RBsurf_s_full, RBsurf_t, RBsurf_s = None, None, None, None

# ============================================================================================
# Creates initial T & S fields and calculates pressure loading
def make_ics(grid, forcing, ini_temp_file, ini_salt_file, pload_file, spinup, prec):

    if spinup:
        # set constant T & S
        sizezyx = (grid.nz, grid.ny, grid.nx)
        ics_T = np.zeros(sizezyx) - 1 # define T & S arrays with correct dimensions
        ics_S = np.zeros(sizezyx) + 34.2

    else:
        # define forcing data
        T = xr.open_dataset(forcing.Tfile)
        S = xr.open_dataset(forcing.Sfile)
        lon, lat, z, time = T.LONGITUDE, T.LATITUDE, T.DEPTH, T.TIME
        nx, ny, nz = lon.shape[0], lat.shape[0], z.shape[0]
        lon2d, lat2d = np.meshgrid(lon, lat)
        lon1d = lon2d.ravel()
        lat1d = lat2d.ravel()

        # T & S arrays for ICs have dimensions (nz, ny, nx)
        sizezyx = (grid.nz, grid.ny, grid.nx)
        ics_T = np.zeros(sizezyx)  # define T & S arrays with correct dimensions
        ics_S = np.zeros(sizezyx)
        ics_T_vertinterp = np.zeros((grid.nz, nx * ny))
        ics_S_vertinterp = np.zeros((grid.nz, nx * ny))

        # Initial T&S conditions are from forcing dataset at time 'starttime'+12months
        # slice full T and S arrays to correct time
        starttime_p1 = starttime + relativedelta(months=+1) #shift starttime to first day of next month to correspond to time
                                                        #conventions in forcing data
        forcingtime = pd.to_datetime(time.TIME.values, format = '%Y-%m-%d')
        startIndex = np.where(starttime_p1 == forcingtime)[0] # find index of starttime within time array of forcing data

        ics_Tslice = T.isel(TIME=np.arange(startIndex,startIndex+13)).THETA.values
        ics_Sslice = S.isel(TIME=np.arange(startIndex,startIndex+13)).SALT.values
        ics_Tslice = np.mean(ics_Tslice, axis=0)
        ics_Sslice = np.mean(ics_Sslice, axis=0)

        ics_Tslice = np.reshape(ics_Tslice, (nz, nx * ny))
        ics_Sslice = np.reshape(ics_Sslice, (nz, nx * ny))

        # Interpolation of vertical levels
        (Ind,Fracs) = interp_functions.vertical_interp(z, -grid.z)

        for t in range(0, np.size(Ind, 1)):
            ics_T_vertinterp[t, :] = Fracs[0, t] * ics_Tslice[Ind[0, t], :] + Fracs[1, t] * ics_Tslice[Ind[1, t], :]
            ics_S_vertinterp[t, :] = Fracs[0, t] * ics_Sslice[Ind[0, t], :] + Fracs[1, t] * ics_Sslice[Ind[1, t], :]

        for k in range(0,grid.nz):
            Tslice = ics_T_vertinterp[k,:]
            Sslice = ics_S_vertinterp[k,:]
            # Use linear barycentric interpolation for horizontal slices, allowing for
            # extrapolation with nearest values where triangles do not exist
            t_xx = interp_functions.horizontal_interp_nonan(lon1d, lat1d, grid.lon2d, grid.lat2d, Tslice)
            s_xx = interp_functions.horizontal_interp_nonan(lon1d, lat1d, grid.lon2d, grid.lat2d, Sslice)
            # assign regridded horizontal slice to full array
            ics_T[k, :, :] = t_xx
            ics_S[k, :, :] = s_xx
            # print some info and step dummy index
            print 'Done ', k + 1, ' out of ', grid.nz

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

print 'Creating topography'
#make_topo(grid, './ua_custom/DataForMIT.mat', input_dir+'bathymetry.shice', input_dir+'shelfice_topo.bin', prec=64, dig_option='bathy')

print 'Reading info on forcing data'
forcing = ForcingInfo()

print 'Creating initial conditions'
#make_ics(grid, forcing, input_dir+'T_ini.bin', input_dir+'S_ini.bin', input_dir+'pload.mdjwf', spinup=1, prec=64)

print 'Creating restoring conditions at open boundaries'
# note that forcing files have prec=32 by default when used in combination with EXF package (exf_iprec=32)
#make_obcs(grid, forcing, input_dir+'OBSt.bin', input_dir+'OBSs.bin', input_dir+'OBSu.bin', input_dir+'OBSv.bin', input_dir+'OBWt.bin', input_dir+'OBWs.bin', input_dir+'OBWu.bin', input_dir+'OBWv.bin', spinup=1, prec=32)

print 'Creating restoring conditions at surface'
# note that forcing files have prec=32 by default when used in combination with EXF package (exf_iprec=32)
make_rbcs(grid, forcing, input_dir+'rbcs_surf_T.bin', input_dir+'rbcs_surf_S.bin', input_dir+'rbcs_mask_T.bin', input_dir+'rbcs_mask_S.bin', spinup=1, prec=32)
