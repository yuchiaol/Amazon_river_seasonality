# ================================================================
# Yu-Chiao @ WHOI Dec 20, 2018
# Surface air temperature plots
# ================================================================

# ================================================================
# Import functions
# ================================================================
import argparse
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import isnan, radians
from IPython import get_ipython
import sys, os, ast
import xlrd
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from pylab import setp, genfromtxt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from datetime import datetime
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import StamenTerrain
import matplotlib.path as mpath
from scipy.signal import detrend

sys.path.append('/home/yliang/lib/python_tools/python_functions/whoi_projects/')
import whoi_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f
import ERA_interim_data_process_f, ORAS4_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/statistics')
import statistical_f, MCA_f, EOF_f
sys.path.append('/home/yliang/lib/python_tools/python_functions')
import plot_functions

# ================================================================
# Define functions
# ================================================================
def read_river_mask():

    print('read river mask')
    river_mask_tmp = np.loadtxt('/home/yliang/whoi_projects/Amazon_river/mask_new/obidosglobalnew.txt')
    river_mask_tmp = np.flipud(river_mask_tmp)

    river_mask = np.zeros((river_mask_tmp.shape))
    river_mask[:,0:361] = river_mask_tmp[:,360:].copy()
    river_mask[:,361:] = river_mask_tmp[:,:360].copy()
    river_mask[river_mask<=0.] = np.nan

    return river_mask

def read_gpcp_precipitation(year_st,year_ed,river_mask):

    print('read gpcp precipitation')
    year_ref = 1979
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GPCP/'
    filename = 'precip.mon.mean.nc'
    f = Dataset(dirname + filename, mode='r')
    lat = f.variables['lat'][:].data
    lon = f.variables['lon'][:].data
    var_tmp = f.variables['precip'][t0:t1,:,:].data
#   var_tmp[abs(var_tmp)>1.e20] = np.nan
    f.close()

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(np.linspace(0,360,721), np.linspace(-90,84,348))
    [lon_x, lat_y] = np.meshgrid(lon, lat)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
    mask_rg[mask_rg!=1] = np.nan

#    lon_amz1 = 300
#    lon_amz2 = 320
#    lat_amz1 = -10
#    lat_amz2 = 10

    area = data_process_f.area_calculate_nonuniform(lon, lat)

#    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)

    nt = var_tmp.shape[0]
    ts_gpcp = np.zeros((nt))
    for NT in range(nt):
        ts_gpcp[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg)/np.nansum(area*mask_rg)
#        ts_gpcp[NT] = np.nansum(var_tmp[NT,y1:y2,x1:x2]*area[y1:y2,x1:x2])/np.nansum(area[y1:y2,x1:x2])

#    plt.contourf(lon,lat,var_tmp[11,:,:]*mask_rg)
#    plt.colorbar()
#    plt.show()

    return ts_gpcp, area

def read_precl_precipitation(year_st,year_ed,river_mask):

    print('read precl precipitation')
    year_ref = 1948
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/PRECL/'
    filename = 'PRECL.mon.mean.1x1.nc'
    f = Dataset(dirname + filename, mode='r')
    lat = f.variables['lat'][:].data
    lat = np.flipud(lat)
    lon = f.variables['lon'][:].data
    var_tmp = f.variables['precip'][t0:t1,:,:].data
    for NT in range(var_tmp.shape[0]):
         var_tmp[NT,:,:] = np.flipud(np.squeeze(var_tmp[NT,:,:]))
#   var_tmp[abs(var_tmp)>1.e20] = np.nan
    f.close()

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(np.linspace(0,360,721), np.linspace(-90,84,348))
    [lon_x, lat_y] = np.meshgrid(lon, lat)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
    mask_rg[mask_rg!=1] = np.nan

    area = data_process_f.area_calculate_nonuniform(lon, lat)

    nt = var_tmp.shape[0]
    ts_out = np.zeros((nt))
    for NT in range(nt):
        ts_out[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg)/np.nansum(area*mask_rg)

    return ts_out

def read_gpcc_precipitation(year_st,year_ed,river_mask):

#    if year_ed > 2010:
#       year_end = 2010
#    else:
    year_end = year_ed

    print('read gpcc precipitation')
    year_ref = 1891
    t0 = (year_st-year_ref)*12
    t1 = (year_end-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GPCC/new/'
#    filename = 'precip.mon.total.1x1.v2018.nc'
    filename = 'precip.comb.v2018to2016-v6monitorafter.total.nc'
    f = Dataset(dirname + filename, mode='r')
    lat = f.variables['lat'][:].data
    lat = np.flipud(lat)
    lon = f.variables['lon'][:].data
    var_tmp = f.variables['precip'][t0:t1,:,:].data
    for NT in range(var_tmp.shape[0]):
         var_tmp[NT,:,:] = np.flipud(np.squeeze(var_tmp[NT,:,:]))
#   var_tmp[abs(var_tmp)>1.e20] = np.nan
    f.close()

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(np.linspace(0,360,721), np.linspace(-90,84,348))
    [lon_x, lat_y] = np.meshgrid(lon, lat)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
    mask_rg[mask_rg!=1] = np.nan

    area = data_process_f.area_calculate_nonuniform(lon, lat)

    nt = var_tmp.shape[0]
    ts_out = np.zeros((nt))
    for NT in range(nt):
        ts_out[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg)/np.nansum(area*mask_rg)

    return ts_out

def read_hybam_river(year_st,year_ed):

    print('read hybam river')
    dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
    filename = 'Amazon_river_obidos_1969_2015.xls'

    year_ref = 1968
    t0 = int((year_st-year_ref)*12)
    t1 = int((year_ed-year_ref)*12+12)

    dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/new/'
    filename = 'Amazon_new_1968_2019.nc'
    f = Dataset(dirname + filename, 'r')
    river_discharge = f.variables['river_discharge_interp'][t0:t1].data
    f.close()

    return river_discharge

def read_hybam_Orinoco_river(year_ed):

    print('read hybam Orinoco river')
    dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
    filename = 'Orinoco_hybam.xls'

    year_ref = 2003
    yearN = year_ed - year_ref + 1
    year_river = np.linspace(year_st,year_ed,yearN)

    wb = xlrd.open_workbook(dirname+filename)
    sheet = wb.sheet_by_index(0)

    river_discharge = np.zeros((yearN*12))

    for II in range(yearN*12):
        river_discharge[II] = sheet.cell_value(II+2,3)

    return river_discharge

def read_Dai_Trenberth_river_data(year_st, year_ed, river_index):
    dirname = '/stormtrack/data4/yliang/observations/streamflow/Dai_Trenberth/'
    filename = 'coastal-stns-Vol-monthly.updated-May2019.nc'

    yearN = year_ed - year_st + 1
    year_river_dai = np.linspace(year_st,year_ed,yearN)
    year_ref = 1900
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    f = Dataset(dirname + filename, mode='r')
    flow = f.variables['FLOW'][t0:t1,river_index].data
    flow[flow<0.] = np.nan
    river_name = f.variables['riv_name'][river_index].data
    print(river_name)
#    time = f.variables['time'][t0:t1].data
    f.close()

    return flow

def read_gecco2_salinity(year_st,year_ed):
    print('read GECCO2 salinity')

    year_ref = 1979

    yearN = int(year_ed - year_st + 1)
    nt = int(yearN*12)

#    lon_amz1 = 295
#    lon_amz2 = 311
#    lat_amz1 = -2
#    lat_amz2 = 15

    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GECCO2/regrided/'
    filename = 'GECCO2_S_regridded_1979_2016.nc'

    f = Dataset(dirname + filename, 'r')
    lat_tmp = f.variables['lat'][:].data.copy()
    lon = f.variables['lon'][:].data.copy()
    so_tmp = f.variables['salt_rg'][t0:t1,:,:].data.copy()
    f.close()

    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,11,-30,30,lon,lat_tmp)

    lat = lat_tmp[y1:y2+1].copy()
    so = so_tmp[:,y1:y2+1,:].copy()

    mask_so = np.nanmean(so, axis=0).copy()/np.nanmean(so, axis=0).copy()
#    mask_so[np.nanmean(so, axis=0)>38] = np.nan 

    area_so = data_process_f.area_calculate_nonuniform(lon, lat)
    area_so[np.isnan(so[1,:,:])==1] = np.nan
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)

# Remove climatology values
#    [so_rm, so_clim] = data_process_f.climatology_anomalies_calculate(so,1981,2010,year_st)

    ts_so = np.zeros((nt))
    for NT in range(nt):
        ts_so[NT] = np.nansum(so[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_so[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_so[y1:y2+1,x1:x2+1])

    return ts_so, lon, lat, mask_so

def read_oras4_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss):

    [so_oras, lon_oras, lat_oras, depth_oras] = ORAS4_data_process_f.read_variable(0,360,-30,30,0,year_st,year_ed,'so')

    print(depth_oras)

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(lon_sss, lat_sss)
    [lon_x, lat_y] = np.meshgrid(lon_oras, lat_oras)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), mask_sss.ravel(), (lat_y, lon_x), method='nearest')

    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    area_so = data_process_f.area_calculate_nonuniform(lon_oras, lat_oras)
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_oras,lat_oras)

#    mask_rg[np.isnan(mask_rg)==False] = 1.

    ts_so = np.zeros((so_oras.shape[0]))
    for NT in range(len(ts_so)):
        ts_so[NT] = np.nansum(so_oras[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

    return ts_so

def read_en4_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss):

    print('EN4 salinity')

    year_ref = 1979
    t0 = (year_st - year_ref)*12
    t1 = (year_ed - year_ref)*12 + 12  

    filename = '/stormtrack/data4/yliang/observations/EN4/EN4_2_1/EN.4.2.1.f.analysis.g10.197901_201812.nc'
    f = Dataset(filename, 'r')
    lon_en4 = f.variables['lon'][:].data
    lat_en4 = f.variables['lat'][:].data
    sss = f.variables['salinity'][t0:t1,0,:,:].data 
    f.close()

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(lon_sss, lat_sss)
    [lon_x, lat_y] = np.meshgrid(lon_en4, lat_en4)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), mask_sss.ravel(), (lat_y, lon_x), method='nearest')

    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    area_so = data_process_f.area_calculate_nonuniform(lon_en4, lat_en4)
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_en4,lat_en4)

    ts_so = np.zeros((sss.shape[0]))
    for NT in range(len(ts_so)):
        ts_so[NT] = np.nansum(sss[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

    return ts_so

def read_ecco_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss):

    year_ref = 1992
    yearN = int(year_ed - year_ref + 1)

    filename = '/stormtrack/data4/yliang/observations/ECCO/SALT.2011.nc'
    f = Dataset(filename, 'r')
    lon_tmp = f.variables['lon'][:,:].data
    lat_tmp = f.variables['lat'][:,:].data
    f.close()

    lon_ecco = lon_tmp[1,:].copy()
    lat_ecco = lat_tmp[:,1].copy()
    
    ny = len(lat_ecco)
    nx = len(lon_ecco)

    sss = np.zeros((yearN*12,ny,nx))
    NN = 0
    for NY in range(yearN):
        for NM in range(12):
            filename = '/stormtrack/data4/yliang/observations/ECCO/ECCO_v4r4/SALT_' + str(int(NY+year_ref)) + '_' + str(NM+1).zfill(2) + '.nc'
            f = Dataset(filename, 'r')
            sss[NN,:,:] = f.variables['SALT'][:,0,:,:].data.squeeze()
            f.close()
            NN = NN + 1

    lon_amz1 = -65
    lon_amz2 = -49
    lat_amz1 = -3
    lat_amz2 = 15

    area_so = data_process_f.area_calculate_nonuniform(lon_ecco, lat_ecco)
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_ecco,lat_ecco)

    mask_rg = np.nanmean(sss, axis=0).copy()/np.nanmean(sss, axis=0).copy()

    ts_so = np.zeros((sss.shape[0]))
    for NT in range(len(ts_so)):
        ts_so[NT] = np.nansum(sss[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

    return ts_so

def read_aquarius_salinity():

    dirname = '/stormtrack/data4/yliang/observations/Aquarius/'
    filename = 'sss201304.v5.0cap.nc'
    f = Dataset(dirname + filename, 'r')
    lon_sss = f.variables['lon'][:].data
    lat_sss = f.variables['lat'][:].data
    f.close()

    sss_tmp = np.zeros((36,len(lat_sss),len(lon_sss)))
    NN = 0
    for NY in range(3):
        for NM in range(12):
            filename = 'sss' + str(2012+NY) + str(NM+1).zfill(2) + '.v5.0cap.nc'
            f = Dataset(dirname + filename, 'r')
            sss_tmp[NN,:,:] = f.variables['sss_cap'][:,:].data.squeeze()
            NN = NN + 1

    sss_tmp[sss_tmp<-100] = np.nan

    mask_rg = np.nanmean(sss_tmp, axis=0).copy()/np.nanmean(sss_tmp, axis=0).copy()

    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    area_so = data_process_f.area_calculate_nonuniform(lon_sss, lat_sss)
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_sss,lat_sss)

    ts_so = np.zeros((sss_tmp.shape[0]))
    for NT in range(len(ts_so)):
        ts_so[NT] = np.nansum(sss_tmp[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

    return ts_so

def read_smos_sss_here():

    filename = '/stormtrack/data4/yliang/observations/SMOS/SMOS_SSS_2010_2016_V3.nc'
    f = Dataset(filename, 'r')
    sss = f.variables['sss_biasadj'][7:,::-1,:].data
    lon_sss = f.variables['lon'][:].data
    lat_sss = f.variables['lat'][::-1].data
    f.close()

    area_sss = data_process_f.area_calculate_nonuniform(lon_sss, lat_sss)

    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(-65,-50,-3,20,lon_sss,lat_sss)

    print(sss.shape)
    ts_sss_smos = np.zeros((sss.shape[0]))
    for NT in range(sss.shape[0]):
        ts_sss_smos[NT] = np.nansum(sss[NT,y1:y2+1,x1:x2+1]*area_sss[y1:y2+1,x1:x2+1])/np.nansum(area_sss[y1:y2+1,x1:x2+1]*sss[NT,y1:y2+1,x1:x2+1]/sss[NT,y1:y2+1,x1:x2+1])

    return ts_sss_smos

def read_gleam_evap_here(year_st,year_ed,river_mask):

    print('read gleam evaporation')
    year_ref = 1980
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GLEAM/new/'
    filename = 'E_1980_2018_GLEAM_v3.3a_MO.nc'
    f = Dataset(dirname + filename, 'r')
    evap = f.variables['E'][t0:t1,:,::-1].data
    lon_evap = f.variables['lon'][:].data
    lat_evap = f.variables['lat'][::-1].data
    f.close()

# Interpolate to precipitation grid
    [grid_x, grid_y] = np.meshgrid(np.linspace(-180,180,721), np.linspace(-90,84,348))
    [lon_x, lat_y] = np.meshgrid(lon_evap, lat_evap)
    mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
    mask_rg[mask_rg!=1] = np.nan

    mask_tmp = mask_rg.copy()*0.
    mask_tmp[:,0:720] = mask_rg[:,720:].copy()
    mask_tmp[:,720:] = mask_rg[:,0:720].copy()

    area = data_process_f.area_calculate_nonuniform(lon_evap, lat_evap)

    ts_evap = np.zeros((evap.shape[0]))
    for NT in range(evap.shape[0]):
        ts_evap[NT] = np.nansum((evap[NT,:,:].T)*area*mask_tmp)/np.nansum(area*mask_tmp)

    return ts_evap

def max_min_select(ts_monthly,year_N):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,:].max()
        ts_min[II+1] = ts_monthly[II+1,:].min()

    return ts_max, ts_min

def max_min_month_select(ts_monthly,year_N,t0_max,t1_max,t0_min,t1_min):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,t0_max:t1_max+1].mean()
        ts_min[II+1] = ts_monthly[II+1,t0_min:t1_min+1].mean()

    return ts_max, ts_min

def plot_box(ax,lon1,lon2,lat1,lat2, color_txt):

    ax.plot(np.linspace(lon1,lon1,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon2,lon2,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat1, lat1, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat2, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)

def plot_figure_here(ax1,lon,lat,map_diff,map_var_diff,map_ttest,clevel,clevel_line,cmap):

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    factor = 1.

    ax1.set_extent([-180, 180, 30, 90], ccrs.PlateCarree())
    map_2D = map_diff.copy()
    map_2D[0:43,:] = np.nan
    [map_ext, lon_ext] = data_process_f.extend_longitude(map_2D,lat,lon)
    im2      = ax1.contourf(lon_ext, lat, map_ext*factor, levels=clevel, extend='both',cmap=cmap, transform=ccrs.PlateCarree())
    map_2D = map_var_diff.copy()
    map_2D[0:43,:] = np.nan
    [map_ext, lon_ext] = data_process_f.extend_longitude(map_2D,lat,lon)
    im      = ax1.contour(lon_ext, lat, map_ext/100, levels=clevel_line, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.8)

    var_sig = map_ttest.copy()
    var_sig[0:43,:] = np.nan
    for II in range(len(lat[::2])):
        ax1.plot(lon, lat[II*2]*var_sig[II*2,:], 'ko', markersize=0.5, markerfacecoloralt='k', transform=ccrs.PlateCarree())
#    plot_box(ax1,50,120,50,70, 'c')
    ax1.coastlines('110m')
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax1.set_aspect('auto')

    return im2

def output_fields_here(filename,river_max,river_min,river3_max,river3_min,prcp_max,prcp_min,sss_max,sss_min,year_N):

    files = os.listdir(os.getcwd())
    for file in files:
        if file == filename:
           print('Delete ' + filename)
           os.system("rm -rf " + file)

    f = Dataset(filename, 'w',format='NETCDF4')
    f.description = 'temporary output fields'
    f.createDimension('dt', year_N)

    river_max_out = f.createVariable('river_max','f4',('dt'))
    river_min_out = f.createVariable('river_min','f4',('dt'))   
    river3_max_out = f.createVariable('river3_max','f4',('dt'))
    river3_min_out = f.createVariable('river3_min','f4',('dt'))
    prcp_max_out = f.createVariable('prcp_max','f4',('dt'))
    prcp_min_out = f.createVariable('prcp_min','f4',('dt'))   
    sss_max_out = f.createVariable('sss_max','f4',('dt'))
    sss_min_out = f.createVariable('sss_min','f4',('dt'))   

    river_max_out[:] = river_max[:]
    river_min_out[:] = river_min[:]
    river3_max_out[:] = river3_max[:]
    river3_min_out[:] = river3_min[:]
    prcp_max_out[:] = prcp_max[:]
    prcp_min_out[:] = prcp_min[:]
    sss_max_out[:] = sss_max[:]
    sss_min_out[:] = sss_min[:]
    f.close()

def read_fields_here(filename):

    f = Dataset(filename, 'r')
    river_max = f.variables['river_max'][:].data
    river_min = f.variables['river_min'][:].data
    river3_max = f.variables['river3_max'][:].data
    river3_min = f.variables['river3_min'][:].data
    prcp_max = f.variables['prcp_max'][:].data
    prcp_min = f.variables['prcp_min'][:].data
    sss_max = f.variables['sss_max'][:].data
    sss_min = f.variables['sss_min'][:].data
    f.close()

    return river_max,river_min,river3_max,river3_min,prcp_max,prcp_min,sss_max,sss_min

def trend_calculate_here(year, ts_temp):

    coeff = np.polyfit(np.linspace(1,len(year[6:-6]),len(year[6:-6])), ts_temp, 1)
    [t, prob] = stats.ttest_ind(np.linspace(1,len(year[6:-6]),len(year[6:-6]))*0, ts_temp)
    sig_label=''
    if prob<0.05: sig_label = '*'

    return coeff, sig_label


def trend_calculate_here_new(year, ts_temp):

    coeff = np.polyfit(np.linspace(1,len(year),len(year)), ts_temp, 1)
    [t, prob] = stats.ttest_ind(np.linspace(1,len(year),len(year))*0, ts_temp)
    sig_label=''
    if prob<0.05: sig_label = '*'

    return coeff, sig_label

# ================================================================
# main starts
# ================================================================
def main(inargs,year_st,year_ed,flag_detrend):
    """Run the program."""

    print('stop here!!!')

# ================================================================
# Describe this script and control flags
# ================================================================
if __name__ == '__main__':
    description='code template'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("run_flag", type=str, help="False: do not run main program, just plotting figure; True: run the main program")

    args = parser.parse_args()

    year_st = 1979
    year_ed = 2018

    flag_detrend = 0

    sig_level = 0.05 

    plt.close('all')

    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if ast.literal_eval(args.run_flag) == True:

# Read river mask
       river_mask = read_river_mask()

# Read GPCP precipitation data
       ts_gpcp, area = read_gpcp_precipitation(year_st,year_ed,river_mask)       

# Read PRECL precipitation data
       ts_precl = read_precl_precipitation(year_st,year_ed,river_mask)

# Read GPCC precipitation data
       ts_gpcc = read_gpcc_precipitation(year_st,year_ed,river_mask)

# Adjust unit
       for NT in range(len(ts_gpcc)):
           ts_gpcc[NT] = ts_gpcc[NT]/month_day[np.mod(NT,12)]       

# Read Hybam river data
       ts_river = read_hybam_river(year_st, year_ed)/100000.
#       print(ts_river)

# Read Hybam Orinoco river data
#       ts_river3 = read_hybam_Orinoco_river(year_ed)
#       print(ts_river3)

# Read Dai_Trenberth river data
       ts_river_dai = read_Dai_Trenberth_river_data(year_st, year_ed, 0)/100000.
#       print(ts_river_dai)
#       ts_river3_dai = read_Dai_Trenberth_river_data(year_st, year_ed, 0)
#       print(ts_river3_dai)

# Combine Orinoco river
#       ts_river3_dai[-144:] = ts_river3.copy()
#       print(ts_river3_dai)

# Read GECCO2 salinity data
       [ts_sss, lon_tmp, lat_sss, mask_tmp] = read_gecco2_salinity(year_st,2016)   

# adjust mask
       mask_sss = mask_tmp.copy()*np.nan
       mask_sss[:,-37:] = mask_tmp[:,0:37].copy()
       mask_sss[:,0:-37] = mask_tmp[:,37:].copy()
       lon_sss = lon_tmp.copy()
       lon_sss[-37:] = lon_tmp[0:37].copy()
       lon_sss[0:-37] = lon_tmp[37:].copy()   

# Read EN4 salinity data
       ts_en4_sss = read_en4_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss) 

# Read Oras4 salinity data
#       ts_oras4_sss = read_oras4_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss)

# Read ECCO salinity
       ts_ecco_sss = read_ecco_salinity(1992,2017,mask_sss,lon_sss,lat_sss)

# Read Aquarius SSS
       ts_aquarius_sss = read_aquarius_salinity()

#       print(np.corrcoef(ts_oras4_sss,ts_sss))

# Read SMOS sea surface salinity
       ts_sss_smos = read_smos_sss_here()


# read GLEAM evaporation data
       ts_evap = read_gleam_evap_here(1980,year_ed,river_mask)      

# Read SODA and ORAS5
       dirname = '/home/yliang/whoi_projects/Amazon_river/paper_prepare/supplementary_figure/'
       filename = 'soda_oras5_temp_output.nc'
       f = Dataset(dirname+filename, 'r')
       sss_oras5_monthly = f.variables['ts_monthly_oras5'][:,:].data
       sss_soda_monthly = f.variables['ts_monthly_soda'][:,:].data

# Perform 3-month running average
       N = 3
       ts_gpcp_rmean = ts_gpcp.copy()
       ts_gpcp_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_gpcp, np.ones((N,))/N, mode='valid')
       ts_precl_rmean = ts_precl.copy()
       ts_precl_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_precl, np.ones((N,))/N, mode='valid')
       ts_gpcc_rmean = ts_gpcc.copy()
       ts_gpcc_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_gpcc, np.ones((N,))/N, mode='valid')
       ts_river_rmean = (ts_river).copy()
       ts_river_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve((ts_river), np.ones((N,))/N, mode='valid')
       ts_river3_rmean = ts_river_dai.copy()
       ts_river3_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_river_dai, np.ones((N,))/N, mode='valid')
       ts_sss_rmean = ts_sss.copy()
       ts_sss_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_sss, np.ones((N,))/N, mode='valid')
       ts_sss_smos_rmean = ts_sss_smos.copy()
       ts_sss_smos_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_sss_smos, np.ones((N,))/N, mode='valid')
#       ts_oras4_sss_rmean = ts_oras4_sss.copy()
#       ts_oras4_sss_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_oras4_sss, np.ones((N,))/N, mode='valid')
       ts_en4_sss_rmean = ts_en4_sss.copy()
       ts_en4_sss_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_en4_sss, np.ones((N,))/N, mode='valid')
       ts_ecco_sss_rmean = ts_ecco_sss.copy()
       ts_ecco_sss_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_ecco_sss, np.ones((N,))/N, mode='valid')
       ts_evap_rmean = ts_evap.copy()
       ts_evap_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_evap, np.ones((N,))/N, mode='valid')

# Arrange to monthly data
       river_monthly = ts_river_rmean.reshape((year_N,12))
       river3_monthly = ts_river3_rmean.reshape((year_N,12))
       prcp_gpcp_monthly = ts_gpcp_rmean.reshape((year_N,12))
       prcp_precl_monthly = ts_precl_rmean.reshape((year_N,12))
       prcp_gpcc_monthly = ts_gpcc_rmean.reshape((year_N,12))
       sss_monthly = ts_sss_rmean.reshape((year_N-2,12))
#       sss_oras4_monthly = ts_oras4_sss_rmean.reshape((year_N,12))
       sss_en4_monthly = ts_en4_sss_rmean.reshape((year_N,12))
       sss_smos_monthly = ts_sss_smos_rmean.reshape((6,12))
       sss_ecco_monthly = ts_ecco_sss_rmean.reshape((year_ed-1992,12))
       sss_aquarius_monthly = ts_aquarius_sss.reshape(3,12)
       evap_monthly = ts_evap_rmean.reshape((year_N-1,12))

# Arrange max and min values
       [river_max, river_min] = max_min_select(river_monthly,year_N)
       print(river_max.shape) 
       print(year[1:-1])
       print(river_max[1:-1])  
       ts_temp = river_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam max', str(reg*10), sig)
       ts_temp = river_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam min', str(reg*10), sig)
       ts_temp = (river_max-river_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('hybam max-min', str(reg*10), sig)

       [river3_max, river3_min] = max_min_select(river3_monthly,year_N)
       print(river3_max.shape)
       print(year[1:-1])
       print(river3_max[1:-1])
       ts_temp = river3_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('dai max', str(reg*10), sig)
       ts_temp = river3_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('dai min', str(reg*10), sig)
       ts_temp = (river3_max-river3_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('dai max-min', str(reg*10), sig)

       [prcp_gpcp_max, prcp_gpcp_min] = max_min_select(prcp_gpcp_monthly,year_N)
       print(prcp_gpcp_max.shape)
       print(year[1:-1])
       print(prcp_gpcp_max[1:-1])
       ts_temp = prcp_gpcp_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcp max', str(reg*10), sig)
       ts_temp = prcp_gpcp_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcp min', str(reg*10), sig)
       ts_temp = (prcp_gpcp_max-prcp_gpcp_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcp max-min', str(reg*10), sig)

       [prcp_gpcc_max, prcp_gpcc_min] = max_min_select(prcp_gpcc_monthly,year_N)
       print(prcp_gpcc_max.shape)
       print(year[1:-1])
       print(prcp_gpcc_max[1:-1])
       ts_temp = prcp_gpcc_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcc max', str(reg*10), sig)
       ts_temp = prcp_gpcc_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcc min', str(reg*10), sig)
       ts_temp = (prcp_gpcc_max-prcp_gpcc_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('gpcc max-min', str(reg*10), sig)

       [prcp_precl_max, prcp_precl_min] = max_min_select(prcp_precl_monthly,year_N)
       print(prcp_precl_max.shape)
       print(year[1:-1])
       print(prcp_precl_max[1:-1])
       ts_temp = prcp_precl_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('precl max', str(reg*10), sig)
       ts_temp = prcp_precl_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('precl min', str(reg*10), sig)
       ts_temp = (prcp_precl_max-prcp_precl_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('precl max-min', str(reg*10), sig)

       [sss_max, sss_min] = max_min_select(sss_monthly,year_N-2)
       print(sss_max.shape)
       print(year[1:-3])
       print(sss_max[1:-1])
       ts_temp = sss_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-3], ts_temp[1:-1], sig_level)
       print('gecco2 max', str(reg*10), sig)
       ts_temp = sss_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-3], ts_temp[1:-1], sig_level)
       print('gecco2 min', str(reg*10), sig)
       ts_temp = (sss_max-sss_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-3], ts_temp[1:-1], sig_level)
       print('gecco2 max-min', str(reg*10), sig)

#       [sss_oras4_max, sss_oras4_min] = max_min_select(sss_oras4_monthly,year_N)
#       print(sss_oras4_max.shape)
#       print(year[1:-1])
#       print(sss_oras4_max[1:-1])
#       ts_temp = sss_oras4_max.copy()
#       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
#       print('oras4 max', str(reg*10), sig)
#       ts_temp = sss_oras4_min.copy()
#       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
#       print('oras4 min', str(reg*10), sig)
#       ts_temp = (sss_oras4_max-sss_oras4_min).copy()
#       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
#       print('oras4 max-min', str(reg*10), sig)

       [sss_oras5_max, sss_oras5_min] = max_min_select(sss_oras5_monthly,year_N)
       print(sss_oras5_max.shape)
       print(year[1:-1])
       print(sss_oras5_max[1:-1])
       ts_temp = sss_oras5_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('oras5 max', str(reg*10), sig)
       ts_temp = sss_oras5_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('oras5 min', str(reg*10), sig)
       ts_temp = (sss_oras5_max-sss_oras5_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('oras5 max-min', str(reg*10), sig)

       [sss_soda_max, sss_soda_min] = max_min_select(sss_soda_monthly,year_N-4)
       print(sss_soda_max.shape)
       print(year[2:-4])
       print(sss_soda_max[1:-1])
       ts_temp = sss_soda_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-4], ts_temp[1:-1], sig_level)
       print('soda max', str(reg*10), sig)
       ts_temp = sss_soda_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-4], ts_temp[1:-1], sig_level)
       print('soda min', str(reg*10), sig)
       ts_temp = (sss_soda_max-sss_soda_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-4], ts_temp[1:-1], sig_level)
       print('soda max-min', str(reg*10), sig)

       [sss_en4_max, sss_en4_min] = max_min_select(sss_en4_monthly,year_N)
       print(sss_en4_max.shape)
       print(year[1:-1])
       print(sss_en4_max[1:-1])
       ts_temp = sss_en4_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('en4 max', str(reg*10), sig)
       ts_temp = sss_en4_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('en4 min', str(reg*10), sig)
       ts_temp = (sss_en4_max-sss_en4_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('en4 max-min', str(reg*10), sig)

       [sss_ecco_max, sss_ecco_min] = max_min_select(sss_ecco_monthly,2017-1992+1)
       print(sss_ecco_max.shape)
       print(year[14:-2])
       print(sss_ecco_max[1:-1])
       ts_temp = sss_ecco_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[14:-2], ts_temp[1:-1], sig_level)
       print('ecco max', str(reg*10), sig)
       ts_temp = sss_ecco_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[14:-2], ts_temp[1:-1], sig_level)
       print('ecco min', str(reg*10), sig)
       ts_temp = (sss_ecco_max-sss_ecco_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[14:-2], ts_temp[1:-1], sig_level)
       print('ecco max-min', str(reg*10), sig)

       [sss_smos_max, sss_smos_min] = max_min_select(sss_monthly,6)
#       [sss_aquarius_max, sss_aquarius_min] = max_min_select(sss_aquarius_monthly,3)

       sss_aquarius_max = np.zeros((3))*np.nan
       sss_aquarius_min = np.zeros((3))*np.nan
       for II in range(3):
           sss_aquarius_max[II] = sss_aquarius_monthly[II,:].max()
           sss_aquarius_min[II] = sss_aquarius_monthly[II,:].min()

       [evap_max, evap_min] = max_min_select(evap_monthly,year_N-1)
       print(evap_max.shape)
       print(year[2:-1])
       print(evap_max[1:-1])
       ts_temp = evap_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('gleam max', str(reg*10), sig)
       ts_temp = evap_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('gleam min', str(reg*10), sig)
       ts_temp = (evap_max-evap_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('gleam max-min', str(reg*10), sig)

# Temperory output for JeffBlue
#       filename = 'precipitation_time_series_temp_output.nc'
#       files = os.listdir(os.getcwd())
#       for file in files:
#           if file == filename:
#               print('Delete ' + filename)
#               os.system("rm -rf " + file)

#       f = Dataset(filename, 'w',format='NETCDF4')
#       f.description = 'temporary output fields'
#       f.createDimension('dt1_long', len(ts_gpcc))
#       f.createDimension('dt1_short', len(prcp_gpcc_max))
#       f.createDimension('dt2_long', len(ts_gpcp))
#       f.createDimension('dt2_short', len(prcp_gpcp_max))

#       ts_gpcc_out = f.createVariable('ts_gpcc','f4',('dt1_long'))
#       ts_gpcc_max_out = f.createVariable('ts_gpcc_max','f4',('dt1_short'))
#       ts_gpcc_min_out = f.createVariable('ts_gpcc_min','f4',('dt1_short'))
#       ts_gpcp_out = f.createVariable('ts_gpcp','f4',('dt2_long'))
#       ts_gpcp_max_out = f.createVariable('ts_gpcp_max','f4',('dt2_short'))
#       ts_gpcp_min_out = f.createVariable('ts_gpcp_min','f4',('dt2_short'))
#       ts_precl_out = f.createVariable('ts_precl','f4',('dt2_long'))
#       ts_precl_max_out = f.createVariable('ts_precl_max','f4',('dt2_short'))
#       ts_precl_min_out = f.createVariable('ts_precl_min','f4',('dt2_short'))

#       ts_gpcc_out[:] = ts_gpcc[:].copy()
#       ts_gpcc_max_out[:] = prcp_gpcc_max[:].copy()
#       ts_gpcc_min_out[:] = prcp_gpcc_min[:].copy()
#       ts_gpcp_out[:] = ts_gpcp[:].copy()
#       ts_gpcp_max_out[:] = prcp_gpcp_max[:].copy()
#       ts_gpcp_min_out[:] = prcp_gpcp_min[:].copy()
#       ts_precl_out[:] = ts_precl[:].copy()
#       ts_precl_max_out[:] = prcp_precl_max[:].copy()
#       ts_precl_min_out[:] = prcp_precl_min[:].copy()

#       f.close()

# ================================================================
# Plot figures
# ================================================================
       fig = plt.figure()
       fig.set_size_inches(10, 10, forward=True)
       window = 21
       sig_level = 0.05

       plt.axes([0.13, 0.8, 0.35, 0.16])
       plt.plot(year,prcp_gpcp_max,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcp_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1) 
       plt.plot(year,prcp_gpcp_min,'r--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcp_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r--', linewidth=1)
       plt.plot(year,prcp_precl_max,'m-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_precl_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'm-', linewidth=1)
       plt.plot(year,prcp_precl_min,'m--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_precl_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'm--', linewidth=1)
       plt.plot(year,prcp_gpcc_max,'y-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcc_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y-', linewidth=1)
       plt.plot(year,prcp_gpcc_min,'y--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcc_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y--', linewidth=1)
       plt.ylabel('(mm/day)', fontsize=10)
       plt.axis([1979, 2018, 3, 10.5])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(a) Precipitation Max and Min', fontsize=10)

#       plt.axes([0.13, 0.7, 0.35, 0.05])
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],(prcp_gpcp_max[1:-1]+prcp_gpcc_max[1:-1]+prcp_precl_max[1:-1])/3.,window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],(prcp_gpcp_min[1:-1]+prcp_gpcc_min[1:-1]+prcp_precl_min[1:-1])/3.,window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       plt.plot([year[0], year[-1]], [0, 0], 'k-', linewidth=0.3)
#       plt.axis([1979, 2014, -2, 2])
#       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['','','','','','',''], fontsize=10)

       plt.axes([0.593, 0.8, 0.35, 0.16])
       plt.plot(year,prcp_gpcp_max-prcp_gpcp_min,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcp_max[1:-1]-prcp_gpcp_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1, label='GPCP')
       plt.plot(year,prcp_precl_max-prcp_precl_min,'m-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_precl_max[1:-1]-prcp_precl_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'm-', linewidth=1, label='PREC/L')
       plt.plot(year,prcp_gpcc_max-prcp_gpcc_min,'y-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],prcp_gpcc_max[1:-1]-prcp_gpcc_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y-', linewidth=1, label='GPCC')
       plt.legend(fontsize='small', ncol=3, loc='upper left')
       plt.ylabel('(mm/day)', fontsize=10)
       plt.axis([1979, 2018, 2, 7.5])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(b) Precipitation Seasonality (Max. minus Min.)', fontsize=10)

#       plt.axes([0.593, 0.7, 0.35, 0.05])
#       ts_tmp = (prcp_gpcp_max[1:-1]-prcp_gpcp_min[1:-1]) + (prcp_gpcc_max[1:-1]-prcp_gpcc_min[1:-1]) + (prcp_precl_max[1:-1]-prcp_precl_min[1:-1])
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],ts_tmp/3.,window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       plt.plot([year[0], year[-1]], [0, 0], 'k-', linewidth=0.3)
#       plt.axis([1979, 2014, -2, 2])
#       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['','','','','','',''], fontsize=10)

       plt.axes([0.13, 0.55, 0.35, 0.16])
       plt.plot(year,river_max,'b-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b-', linewidth=1)
       plt.plot(year,river_min,'b--', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b--', linewidth=1)
       plt.plot(year,river3_max,'c-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river3_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'c-', linewidth=1)
       plt.plot(year,river3_min,'c--', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river3_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'c--', linewidth=1)
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
       plt.axis([1979, 2018, 0.5, 3])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(c) Amazon River Discharge Max. and Min.', fontsize=10)

#       plt.axes([0.13, 0.38, 0.35, 0.05])
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],river_max[1:-1],window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],river_min[1:-1],window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       plt.plot([year[0], year[-1]], [0, 0], 'k-', linewidth=0.3)
#       plt.axis([1979, 2014, -1, 1])
#       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['','','','','','',''], fontsize=10)

       plt.axes([0.593, 0.55, 0.35, 0.16])
       plt.plot(year,river_max-river_min, 'b-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river_max[1:-1]-river_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'b-', linewidth=1, label='HYBAM')
       plt.plot(year,river3_max-river3_min, 'c-', linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],river3_max[1:-1]-river3_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'c-', linewidth=0.5, label='Dai-Trenberth')
       plt.legend(loc='upper left', ncol=1, fontsize='small')
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
#       plt.axis([1979, 2014, 0.7, 2])
       plt.xlim(1979,2018)
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(d) Amazon River Discharge Seasonality (Max. minus Min. x10$^5$)', fontsize=10)

#       plt.axes([0.593, 0.38, 0.35, 0.05])
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],river_max[1:-1]-river_min[1:-1],window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       plt.plot([year[0], year[-1]], [0, 0], 'k-', linewidth=0.3)
#       plt.axis([1979, 2014, -1, 1])
#       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['','','','','','',''], fontsize=10)

       tt_ecco = np.linspace(1992,2017,2017-1992+1)
       ax1 = fig.add_axes([0.13, 0.3, 0.35, 0.16])
       ax1.plot(year[:-2],sss_max,'b-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-3],sss_max[1:-1],1)
       ax1.plot(year[1:-3], corrcoef[0]*year[1:-3]+corrcoef[1], 'b-', linewidth=1)
       ax1.plot(year[:-2],sss_min,'b--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-3],sss_min[1:-1],1)
       ax1.plot(year[1:-3], corrcoef[0]*year[1:-3]+corrcoef[1], 'b--', linewidth=1)

       ax1.plot(tt_ecco,sss_ecco_max,'g-',linewidth=0.5)
       corrcoef = np.polyfit(tt_ecco[1:-1],sss_ecco_max[1:-1],1)
       ax1.plot(tt_ecco[1:-1], corrcoef[0]*tt_ecco[1:-1]+corrcoef[1], 'g-', linewidth=1)
       ax1.plot(tt_ecco,sss_ecco_min,'g--',linewidth=0.5)
       corrcoef = np.polyfit(tt_ecco[1:-1],sss_ecco_min[1:-1],1)
       ax1.plot(tt_ecco[1:-1], corrcoef[0]*tt_ecco[1:-1]+corrcoef[1], 'g--', linewidth=1)

       ax1.set_ylabel('(PSU)', color='k', fontsize=10)
       plt.xlim(1979,2018)

       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(e) 5-meter Sea Salinity Max. and Min.', fontsize=10)

       ax1 = fig.add_axes([0.13, 0.05, 0.35, 0.16])
       ax1.plot(year,sss_oras5_max,'k-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_max[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'k-', linewidth=1)
       ax1.plot(year,sss_oras5_min,'k--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_min[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'k--', linewidth=1)

       ax1.plot(year,sss_en4_max,'y-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_en4_max[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y-', linewidth=1)
       ax1.plot(year,sss_en4_min,'y--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_en4_min[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y--', linewidth=1)

       ax1.plot(year[1:-3],sss_soda_max,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-4],sss_soda_max[1:-1],1)
       ax1.plot(year[2:-4], corrcoef[0]*year[2:-4]+corrcoef[1], 'r-', linewidth=1)
       ax1.plot(year[1:-3],sss_soda_min,'r--',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-4],sss_soda_min[1:-1],1)
       ax1.plot(year[2:-4], corrcoef[0]*year[2:-4]+corrcoef[1], 'r--', linewidth=1)

       ax1.set_ylabel('(PSU)', color='k', fontsize=10)

#       plt.axis([1979, 2014, 33.3, 36])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(g) 5-meter Sea Salinity Max. and Min.', fontsize=10)
       plt.xlim(1979,2018)

#       plt.axes([0.13, 0.06, 0.35, 0.05])
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],sss_max[1:-1],window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       [y_slope, y_sig] = statistical_f.running_linear_trends(year[1:-1],sss_min[1:-1],window,sig_level)
#       y_sig[y_sig==0] = np.nan
#       plt.bar(year[1:-1],y_slope*20)
#       for NT in range(len(year[1:-1])):
#           plt.plot([year[1+NT], year[1+NT]], [y_slope[NT]*y_sig[NT]*20, y_slope[NT]*y_sig[NT]*20], 'k*', markersize=1.6)
#       plt.plot([year[0], year[-1]], [0, 0], 'k-', linewidth=0.3)
#       plt.axis([1979, 2014, -0.3, 0.3])
#       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['','','','','','',''], fontsize=10)

       ax1 = fig.add_axes([0.593, 0.3, 0.35, 0.16])
       ax1.plot(year[:-2],sss_max-sss_min,'b-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-3],sss_max[1:-1]-sss_min[1:-1],1)
       ax1.plot(year[1:-3], corrcoef[0]*year[1:-3]+corrcoef[1], 'b-', linewidth=1, label='GECCO2')
       plt.xlim(1979,2018)

       ax1.tick_params('y', colors='b')
       ax1.set_ylabel('(PSU)', color='b', fontsize=10)
       ax1.legend(loc='upper left', ncol=3, fontsize='small')
       ax2 = ax1.twinx()
       ax2.plot(tt_ecco, sss_ecco_max-sss_ecco_min, 'g-', linewidth=0.5)
       corrcoef = np.polyfit(tt_ecco[1:-1],sss_ecco_max[1:-1]-sss_ecco_min[1:-1],1)
       ax2.plot(tt_ecco[1:-1], corrcoef[0]*tt_ecco[1:-1]+corrcoef[1], 'g-', linewidth=1, label='ECCO4')
       ax2.tick_params('y', colors='g')
       ax2.set_ylabel('', color='g', fontsize=8)
       ax2.legend(loc='lower right', ncol=1, fontsize='small')
#       plt.axis([1979, 2014, 0.5, 3.5])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(f) 5-meter Sea Salnity Seasonality (Max. minus Min.)', fontsize=10)
       plt.xlim(1979,2018)


       ax1 = fig.add_axes([0.593, 0.05, 0.35, 0.16])
       ax1.plot(year,sss_oras5_max-sss_oras5_min,'k-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_max[1:-1]-sss_oras5_min[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'k-', linewidth=1, label='ORAS5')

       ax1.plot(year,sss_en4_max-sss_en4_min,'y-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_en4_max[1:-1]-sss_en4_min[1:-1],1)
       ax1.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'y-', linewidth=1, label='EN4')

       ax1.plot(year[1:-3],sss_soda_max-sss_soda_min,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-4],sss_soda_max[1:-1]-sss_soda_min[1:-1],1)
       ax1.plot(year[2:-4], corrcoef[0]*year[2:-4]+corrcoef[1], 'r-', linewidth=1, label='SODA3.1.1')

       ax1.legend(loc='lower right', ncol=3, fontsize='small')
       ax1.set_ylabel('(PSU)', color='k', fontsize=10)

       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015],['1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=10)
       plt.title('(h) 5-meter Sea Salnity Seasonality (Max. minus Min.)', fontsize=10)
       plt.xlim(1979,2018)


       plt.savefig('Fig2_time_series_all_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

    else:

       print('hahahah')


