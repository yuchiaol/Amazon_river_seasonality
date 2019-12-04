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

    dirname = '/stormtrack/data4/yliang/observations/GPCC/'
#    filename = 'GPCC_V6_190101-201012_360X720.nc'
    filename = 'precip.mon.total.1x1.v2018.nc'
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

    yearN = year_ed - year_st + 1
    year_river = np.linspace(year_st,year_ed,yearN)
    year_ref = 1969
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    wb = xlrd.open_workbook(dirname+filename)
    sheet = wb.sheet_by_index(0)

    river_discharge = np.zeros(((2015-1969+1)*12))

    for II in range((2015-1969+1)*12):
        river_discharge[II] = sheet.cell_value(II+2,3)

    return river_discharge[t0:t1]

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
    filename = 'coastal-stns-Vol-monthly.updated-Aug2014.nc'

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
    filename = 'GECCO2_S_regridded.nc'

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

    filename = '/stormtrack/data4/yliang/observations/EN4/EN4_2_0_197901_201612.nc'
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

    filename = '/stormtrack/data4/yliang/observations/ECCO/SALT.1994.nc'
    f = Dataset(filename, 'r')
    lon_tmp = f.variables['lon'][:,:].data
    lat_tmp = f.variables['lat'][:,:].data
    f.close()

    lon_ecco = lon_tmp[1,:].copy()
    lat_ecco = lat_tmp[:,1].copy()
    
    ny = len(lat_ecco)
    nx = len(lon_ecco)

    sss = np.zeros((yearN*12,ny,nx))
    for NY in range(yearN):
        filename = '/stormtrack/data4/yliang/observations/ECCO/SALT.' + str(int(NY+year_ref)) + '.nc'
        f = Dataset(filename, 'r')
        sss[12*NY:12*NY+12,:,:] = f.variables['SALT'][:,0,:,:].data

    lon_amz1 = -65
    lon_amz2 = -49
    lat_amz1 = -3
    lat_amz2 = 10

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
    lat_amz2 = 20

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
    year_ref = 1979
    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GLEAM/'
    filename = 'E_monthly_GLEAM_1980_2017.nc'
    f = Dataset(dirname + filename, 'r')
    evap = f.variables['E'][t0:t1,:,:].data
    lon_evap = f.variables['lon'][:].data
    lat_evap = f.variables['lat'][:].data
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
        ts_evap[NT] = np.nansum(evap[NT,:,:]*area*mask_tmp)/np.nansum(area*mask_tmp)

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
    year_ed = 2014

    flag_detrend = 0

    sig_level = 0.05 

    plt.close('all')

    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if ast.literal_eval(args.run_flag) == True:

# ================================================================
# Read SODA SSS
# ================================================================
       dirname = '/stormtrack/data4/yliang/observations/SODA/'
       filename = 'soda3.3.1_mn_ocean_reg_1980.nc'
       f = Dataset(dirname + filename, mode='r')
       lat = f.variables['latitude'][:].data
       lon = f.variables['longitude'][:].data
       var_tmp = f.variables['salt'][:,0,:,:].data.squeeze()
       var_tmp[abs(var_tmp)>1.e19] = np.nan
       f.close()

       mask_rg = var_tmp[0,:,:]/var_tmp[0,:,:]
       area = data_process_f.area_calculate_nonuniform(lon, lat)

       ny = len(lat)
       nx = len(lon)

       var_test = np.zeros(((year_N-1)*12,ny,nx))
       for NY in range(year_N-1):
               filename = 'soda3.3.1_mn_ocean_reg_' + str(int(1980+NY)) + '.nc'
               print(filename)
               f = Dataset(dirname + filename, mode='r')
               var_tmp = f.variables['salt'][:,0,:,:].data.squeeze()
               var_tmp[abs(var_tmp)>1.e19] = np.nan
               f.close()
               var_test[0+NY*12:12+NY*12,:,:] = var_tmp.copy()

# Calculate time series
       lon_amz1 = 295
       lon_amz2 = 311
       lat_amz1 = -3
       lat_amz2 = 15

       [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)

       ts_test = np.zeros((year_N*12-12))
       for NT in range(year_N*12-12):
           ts_test[NT] = np.nansum(var_test[NT,y1:y2+1,x1:x2+1]*area[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

# Perform 3-month running average
       N = 3
       ts_test_rmean = ts_test.copy()
       ts_test_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_test, np.ones((N,))/N, mode='valid')

# Arrange to monthly data
       ts_monthly_soda = ts_test_rmean.reshape((year_N-1,12))

# Arrange max and min values
       [sss_soda_max, sss_soda_min] = max_min_select(ts_monthly_soda,year_N-1)
       print(sss_soda_max.shape)
       print(year[1:-1])
       print(sss_soda_max[1:-1])
       ts_temp = sss_soda_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('SODA max', str(reg*10), sig)
       ts_temp = sss_soda_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('SODA min', str(reg*10), sig)
       ts_temp = (sss_soda_max-sss_soda_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[2:-1], ts_temp[1:-1], sig_level)
       print('SODA max-min', str(reg*10), sig)

# ================================================================
# Read ORAS5 data
# ================================================================
       dirname = '/stormtrack/data4/yliang/observations/ORAS5/oras5/ORCA025/vosaline/opa3/'
       filename = 'vosaline_ORAS5_1m_197901_grid_T_02.nc'
       f = Dataset(dirname + filename, mode='r')
       nav_lat = f.variables['nav_lat'][:,:].data
       nav_lon = f.variables['nav_lon'][:,:].data
       var_tmp = f.variables['vosaline'][:,0,:,:].data.squeeze()
       var_tmp[var_tmp>1.e30] = np.nan
       f.close()

       lon_rg = np.linspace(-180,180,360)
       lat_rg = np.linspace(-77,89,166)
       [lon_x, lat_y] = np.meshgrid(lon_rg, lat_rg)
       var_tmp_rg = griddata((nav_lat.ravel(), nav_lon.ravel()), var_tmp.ravel(), (lat_y, lon_x), method='nearest')

       mask_rg = var_tmp_rg/var_tmp_rg
       area = data_process_f.area_calculate_nonuniform(lon_rg, lat_rg)

       ny = len(lat_rg)
       nx = len(lon_rg)
        
       var_test = np.zeros((year_N*12,ny,nx))
       NN = 0
       for NY in range(year_N):
           for NM in range(12):
               filename = 'vosaline_ORAS5_1m_' + str(int(NY+1979)) + str(NM+1).zfill(2) + '_grid_T_02.nc'
               print(filename) 
               f = Dataset(dirname + filename, mode='r')
               var_tmp = f.variables['vosaline'][:,0,:,:].data.squeeze()
               var_tmp[var_tmp>1.e30] = np.nan
               f.close()
               var_tmp_rg = griddata((nav_lat.ravel(), nav_lon.ravel()), var_tmp.ravel(), (lat_y, lon_x), method='nearest')
               var_test[NN,:,:] = var_tmp_rg.copy()
               NN = NN + 1

# Calculate time series
       lon_amz1 = -65
       lon_amz2 = -49
       lat_amz1 = -3
       lat_amz2 = 15

       [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_rg,lat_rg)

       ts_test = np.zeros((year_N*12))
       for NT in range(year_N*12):
           ts_test[NT] = np.nansum(var_test[NT,y1:y2+1,x1:x2+1]*area[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])/np.nansum(area[y1:y2+1,x1:x2+1]*mask_rg[y1:y2+1,x1:x2+1])

# Perform 3-month running average
       N = 3
       ts_test_rmean = ts_test.copy()
       ts_test_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_test, np.ones((N,))/N, mode='valid')

# Arrange to monthly data
       ts_monthly_oras5 = ts_test_rmean.reshape((year_N,12))

# Arrange max and min values
       [sss_oras5_max, sss_oras5_min] = max_min_select(ts_monthly_oras5,year_N)
       print(sss_oras5_max.shape) 
       print(year[1:-1])
       print(sss_oras5_max[1:-1])  
       ts_temp = sss_oras5_max.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('ORAS5 max', str(reg*10), sig)
       ts_temp = sss_oras5_min.copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('ORAS5 min', str(reg*10), sig)
       ts_temp = (sss_oras5_max-sss_oras5_min).copy()
       [reg, sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(year[1:-1], ts_temp[1:-1], sig_level)
       print('ORAS5 max-min', str(reg*10), sig)

# ================================================================
# Output file
# ================================================================
       filename = 'soda_oras5_temp_output.nc'
       print(filename)
       files = os.listdir(os.getcwd())
       for file in files:
           if file == filename:
              print('Delete ' + filename)
              os.system("rm -rf " + file)

       f = Dataset(filename, 'w',format='NETCDF4')
       f.description = 'temporary output fields'
       f.createDimension('d1', ts_monthly_oras5.shape[0])
       f.createDimension('d2', ts_monthly_soda.shape[0])
       f.createDimension('d3', 12)

       ts_monthly_oras5_out = f.createVariable('ts_monthly_oras5','f4',('d1','d3'))
       ts_monthly_soda_out = f.createVariable('ts_monthly_soda','f4',('d2','d3'))

       ts_monthly_oras5_out[:,:] = ts_monthly_oras5[:,:].copy()
       ts_monthly_soda_out[:,:] = ts_monthly_soda[:,:].copy() 

       f.close()

# ================================================================
# Plot figures
# ================================================================
if True:
       fig = plt.figure()
       fig.set_size_inches(10, 10, forward=True)
       window = 21
       sig_level = 0.05

       plt.axes([0.13, 0.7, 0.35, 0.2])
       plt.plot(year,sss_oras5_max,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1)
       plt.plot(year,sss_oras5_min,'r--',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r--', linewidth=1)
       plt.ylabel('(PSU)', fontsize=10)
       plt.axis([1979, 2014, 31, 35])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(a) ORAS5 5-meter Salinity Max and Min', fontsize=10)

       plt.axes([0.593, 0.7, 0.35, 0.2])
       plt.plot(year,sss_oras5_max-sss_oras5_min,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[1:-1],sss_oras5_max[1:-1]-sss_oras5_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1, label='ORAS5')
       plt.legend(fontsize=10, ncol=3, loc='lower left')
       plt.ylabel('(PSU)', fontsize=10)
       plt.axis([1979, 2014, 1, 4])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(b) ORAS5 5-meter Salinity Seasonality (Max. minus Min.)', fontsize=10)

       plt.axes([0.13, 0.4, 0.35, 0.2])
       plt.plot(year[1:],sss_soda_max,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-1],sss_soda_max[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1)
       plt.plot(year[1:],sss_soda_min,'r--',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-1],sss_soda_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r--', linewidth=1)
       plt.ylabel('(PSU)', fontsize=10)
       plt.axis([1979, 2014, 32, 36])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(c) SODA 5-meter Salinity Max and Min', fontsize=10)

       plt.axes([0.593, 0.4, 0.35, 0.2])
       plt.plot(year[1:],sss_soda_max-sss_soda_min,'r-',linewidth=0.5)
       corrcoef = np.polyfit(year[2:-1],sss_soda_max[1:-1]-sss_soda_min[1:-1],1)
       plt.plot(year[1:-1], corrcoef[0]*year[1:-1]+corrcoef[1], 'r-', linewidth=1, label='SODA')
       plt.legend(fontsize=10, ncol=3, loc='lower left')
       plt.ylabel('(PSU)', fontsize=10)
       plt.axis([1979, 2014, 0, 3])
       plt.xticks([1980, 1985, 1990, 1995, 2000, 2005, 2010],['1980','1985','1990','1995','2000','2005','2010'], fontsize=10)
       plt.title('(d) SODA 5-meter Salinity Seasonality (Max. minus Min.)', fontsize=10)

       plt.savefig('sss_time_series_all_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

       print('hahahah')


