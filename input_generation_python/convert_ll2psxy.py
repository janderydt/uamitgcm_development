# ============================================================================================
import numpy as np
from pyproj import Proj

# ============================================================================================
def uv_2psxy(U,V,lon,lat):

    ## Convert U and V from lat lon to ps
    ## assume dimensions U(n,x), V(n,x), lon(x), lat(x)

    delt = 1
    Cm = -110
    Ups = np.zeros(np.shape(U))
    Vps = np.zeros(np.shape(V))

    # expand arrays for vectorization
    lon = np.repeat(lon[np.newaxis,:], np.shape(U)[0], axis=0)
    lat = np.repeat(lat[np.newaxis,:], np.shape(U)[0], axis=0)

    # get end points in polar stereo from velocities (m)
    endlon = lon + U/(2 * np.pi * np.cos(lat*np.pi/180) * 6378137 / 360) * delt
    endlat = lat + V/(2 * np.pi * 6378137 / 360) * delt

    # convert end points to lat - lon
    p = Proj('+init=EPSG:3031')
    x, y = p(lon-Cm, lat)
    endx, endy = p(endlon-Cm, endlat)

    # calculate lat - lon displacements and convert to metres and then velocities
    Ups = (endx - x) / delt
    Vps = (endy - y) / delt

    return (Ups,Vps)