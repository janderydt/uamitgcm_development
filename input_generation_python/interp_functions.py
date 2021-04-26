# ============================================================================================
import numpy as np
from scipy import interpolate
import pandas as pd

# ============================================================================================
def time_interp_reg(original_timeseries, regular_interval):
    """ Input1: Original timeseries, assumed to be numpy array with dtype='datetime64[ns]'
        Input2: Regular output interval, provided as integer in units of days,
            e.g. regular_interval=30 will interpolate onto regular 30 day intervals

        Function will calculate regular output times that fit the original time interval, starting at
        original_timeseries[0]

        Output consists of index array 'ind' and fractions 'fracs', where
            t_reg[i] = fracs[0,i]*t_orig[ind[0,i]] + fracs[1,i]*t_orig[ind[1,i]]
    """
    torig = pd.to_datetime(original_timeseries)
    nt = np.shape(torig)[0]
    df = pd.DataFrame(np.arange(0, nt), index=torig, columns=list('I'))
    df = df.resample('D').mean().interpolate('linear')
    teq = df.resample(np.str(regular_interval)+'D').nearest()

    ind = (np.floor(teq['I']), np.ceil(teq['I']))
    ind = np.asarray(ind)
    fracs = (ind[1, :] - teq['I'], teq['I'] - ind[0, :])
    fracs = np.asarray(fracs)
    fracs[(fracs == 0) & (ind[0, :] == ind[1, :])] = 0.5
    ind = ind.astype('int')

    return (ind,fracs)


# ============================================================================================
def vertical_interp(original_depth, interpolated_depth):
    """ Input1: Original depth levels
        Input2: Target depth levels

        Output consists of index array 'ind' and fractions 'fracs', where
            z_out[i] = fracs[0,i]*z_in[ind[0,i]] + fracs[1,i]*z_in[ind[1,i]]
    """

    # linearly interpolate vertical levels first
    nz = np.shape(original_depth)[0]
    FI = interpolate.interp1d(original_depth, np.arange(0, nz), kind='linear')
    I = FI(interpolated_depth)
    ind = (np.floor(I), np.ceil(I))
    ind = np.asarray(ind)
    fracs = (ind[1, :] - I, I - ind[0, :])
    fracs = np.asarray(fracs)
    fracs[(fracs == 0) & (ind[0, :] == ind[1, :])] = 0.5
    ind = ind.astype('int')

    return (ind,fracs)

# ============================================================================================
def horizontal_interp_nonan(lon_in_1d, lat_in_1d, lon_out_2d, lat_out_2d, var_in_1d):
    """ Interpolates one-dimension data horizontally to a 2d numpy array.

        Method: triangular linear barycentryc interpolation, NOT using nans (i.e. find triangle with non-nan values)
                and nearest-neighbor interpolations for points not surrounded by 3 data points.

        lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate

        var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)

        lon\_out\_2d, lat\_out\_2d: 2d longitude and latitude of the target grid
    """

    # eliminate zeros as these will lead to unrealistic T & S values during interpolation
    x = lon_in_1d[var_in_1d != 0]
    y = lat_in_1d[var_in_1d != 0]
    var_in_1d = var_in_1d[var_in_1d != 0]

    if x.size != 0 :
        p_xx = interpolate.griddata((y, x), var_in_1d, (lat_out_2d, lon_out_2d), method='linear', fill_value=np.nan)
        p_ss = interpolate.griddata((y, x), var_in_1d, (lat_out_2d, lon_out_2d), method='nearest')
        p_xx[(np.isnan(p_xx))] = p_ss[(np.isnan(p_xx))]
    else:
        p_xx = np.zeros(np.shape(lat_out_2d))

    return p_xx

