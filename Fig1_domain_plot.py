flag_run = 1
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
import ERA_interim_data_process_f
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

    print(river_mask.shape)

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

    return ts_gpcp

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

    if year_ed > 2010:
       year_end = 2010
    else:
       year_end = year_ed.copy()

    print('read gpcc precipitation')
    year_ref = 1901
    t0 = (year_st-year_ref)*12
    t1 = (year_end-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GPCC/'
    filename = 'GPCC_V6_190101-201012_360X720.nc'
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

    lon_amz1 = 295
    lon_amz2 = 310
    lat_amz1 = -3
    lat_amz2 = 15

#    lon_amz1 = 300
#    lon_amz2 = 318
#    lat_amz1 = -2
#    lat_amz2 = 11

    t0 = (year_st-year_ref)*12
    t1 = (year_ed-year_ref)*12+12

    dirname = '/stormtrack/data4/yliang/observations/GECCO2/regrided/'
    filename = 'GECCO2_S_regridded.nc'

    f = Dataset(dirname + filename, 'r')
    lat = f.variables['lat'][:].data.copy()
    lon = f.variables['lon'][:].data.copy()
    so = f.variables['salt_rg'][t0:t1,:,:].data.copy()
    f.close()

    so_mask = so[111,:,:].copy()*0.
    so_mask_tmp = so_mask.copy()
    so_mask_tmp[np.isnan(so_mask_tmp)==False] = 1

    area_so = data_process_f.area_calculate_nonuniform(lon, lat)
    area_so[np.isnan(so[1,:,:])==1] = np.nan
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)

# Remove climatology values
#    [so_rm, so_clim] = data_process_f.climatology_anomalies_calculate(so,1981,2010,year_st)

    ts_so = np.zeros((nt))
    for NT in range(nt):
        ts_so[NT] = np.nansum(so[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1])

    so_mask[y1-1:y2+1,x1:x2+1] = 1.
    so_mask = so_mask_tmp*so_mask
    so_mask[np.isnan(so_mask)==True] = 0.

#    plt.contourf(lon,lat,so_mask)
#    plt.colorbar()
#    plt.show()

#    sys.exit()

    return ts_so, lat, lon, so_mask

def read_ecco_salinity(year_st,year_ed):

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
    lat_amz2 = 15

    area_so = data_process_f.area_calculate_nonuniform(lon_ecco, lat_ecco)
    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon_ecco,lat_ecco)

    print(lon_ecco)

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

    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(-65,-50,-2,15,lon_sss,lat_sss)

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

def bar_plot_here(ax, tt, ts_input, color_txt):

    for NM in range(12):
        std_value = ts_input[:,NM].std()
        mean_value = np.nanmean(ts_input[:,NM])
        ax.plot([tt[NM], tt[NM]], [-std_value, std_value]+mean_value, color_txt + '-', linewidth=0.5)
        ax.plot([tt[NM]-0.1, tt[NM]+0.1], [-std_value, -std_value]+mean_value, color_txt + '-', linewidth=0.5)
        ax.plot([tt[NM]-0.1, tt[NM]+0.1], [std_value, std_value]+mean_value, color_txt + '-', linewidth=0.5)

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

# ================================================================
# main starts
# ================================================================
if flag_run == 1:

    plt.close('all')
    year_st = 1979
    year_ed = 2014


    year_N = year_ed - year_st + 1
    year = np.linspace(year_st,year_ed,year_N)

    month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if True:

# Read river mask
       river_mask = read_river_mask()

#       plt.figure()
#       plt.contourf(river_mask)
#       plt.show()

# Read GPCP precipitation data
       ts_gpcp = read_gpcp_precipitation(year_st,year_ed,river_mask)       

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
       [ts_sss, lat_sss, lon_sss, so_mask] = read_gecco2_salinity(year_st,year_ed)   

# Read ECCO salinity data
       ts_sss_ecco = read_ecco_salinity(year_st,year_ed)

# Read SMOS sea surface salinity
       ts_sss_smos = read_smos_sss_here()

# read GLEAM evaporation data
#       ts_evap = read_gleam_evap_here(1980,year_ed,river_mask)      

# Read Aquarius sss
       ts_aquarius_sss = read_aquarius_salinity()

# Read EN4
       [ts_xxx, lat_sss, lon_tmp, mask_tmp] = read_gecco2_salinity(year_st,year_ed)

       mask_sss = mask_tmp.copy()*np.nan
       mask_sss[:,-37:] = mask_tmp[:,0:37].copy()
       mask_sss[:,0:-37] = mask_tmp[:,37:].copy()
       lon_sss = lon_tmp.copy()
       lon_sss[-37:] = lon_tmp[0:37].copy()
       lon_sss[0:-37] = lon_tmp[37:].copy()

       ts_en4_sss = read_en4_salinity(year_st,year_ed,mask_sss,lon_sss,lat_sss)

# Read SODA and ORAS5
       dirname = '/home/yliang/whoi_projects/Amazon_river/paper_prepare/supplementary_figure/'
       filename = 'soda_oras5_temp_output.nc'
       f = Dataset(dirname+filename, 'r')
       sss_oras5_monthly = f.variables['ts_monthly_oras5'][:,:].data
       sss_soda_monthly = f.variables['ts_monthly_soda'][:,:].data
       f.close()

# Arrange to monthly data
       river_monthly = ts_river.reshape((year_N,12))
       river3_monthly = ts_river_dai.reshape((year_N,12))
       prcp_gpcp_monthly = ts_gpcp.reshape((year_N,12))
       prcp_precl_monthly = ts_precl.reshape((year_N,12))
       prcp_gpcc_monthly = ts_gpcc.reshape((year_N-(year_ed-2010),12))
       sss_monthly = ts_sss.reshape((year_N,12))
       sss_ecco_monthly = ts_sss_ecco.reshape(((year_ed-1992+1),12))
       sss_smos_monthly = ts_sss_smos.reshape((6,12))
#       evap_monthly = ts_evap.reshape((year_N-1,12))
       sss_aquarius_monthly = ts_aquarius_sss.reshape((3,12))
       sss_en4_monthly = ts_en4_sss.reshape((year_N,12))

       print(river_monthly.shape)
 
# ================================================================
# Read trmm and cmap
# ================================================================
       dirname = '/home/yliang/whoi_projects/Amazon_river/observations/trmm_cmap/'
       filename = 'ts_trmm_cmap_temp.nc'
       f = Dataset(dirname+filename, 'r')
       ts_cmap = f.variables['ts_cmap'][:].data
       ts_trmm = f.variables['ts_trmm'][:].data
       f.close()

       cmap_monthly = ts_cmap.reshape((int(len(ts_cmap)/12),12))
       trmm_monthly = ts_trmm.reshape((int(len(ts_trmm)/12),12))

# ================================================================
# Read GECCO salinity, U, and V
# ================================================================
       dirname = '/stormtrack/data4/yliang/observations/GECCO2/regrided/'
       filename = 'GECCO2_S_regridded_1979_2016.nc'
       f = Dataset(dirname+filename, 'r')
       lon_salt = f.variables['lon'][:].data
       lat_salt = f.variables['lat'][:].data
       salt_rg = np.nanmean(f.variables['salt_rg'][:,:,:].data, axis=0)
       f.close()
#       salt_rg[np.isnan(salt_rg)] = 0.

       filename = 'GECCO2_U_regridded.nc'
       f = Dataset(dirname+filename, 'r')
       lon_u = f.variables['lon'][:].data
       lat_u = f.variables['lat'][:].data
       u_rg = np.nanmean(f.variables['U_rg'][:,:,:].data, axis=0)
       f.close()

       filename = 'GECCO2_V_regridded.nc'
       f = Dataset(dirname+filename, 'r')
       lon_v = f.variables['lon'][:].data
       lat_v = f.variables['lat'][:].data
       v_rg = np.nanmean(f.variables['V_rg'][:,:,:].data, axis=0)
       f.close()

# ================================================================
# Plot figures
# ================================================================
if True:
       font_size = 10

       plt.close('all')

       fig = plt.figure()
       fig.set_size_inches(10, 7, forward=True)

       month_txt = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
       tt = np.linspace(1,12,12)
       scale_in = 2
       v_scale = 5

# Figure 1(a)
       lon_obidos = -55.6753
       lat_obidos = -1.9225
       lon_bolivar = -63.5361
       lat_bolivar = 8.1536
       ax1 = plt.axes([0.08, 0.14, 0.42, 0.76], projection=ccrs.PlateCarree())
       ax1.set_extent([-80, -45, -20, 20])
       im = ax1.contourf(lon_salt-360,lat_salt,salt_rg, levels=np.linspace(34.5,37,15), projection=ccrs.PlateCarree(), cmap='viridis', extend='both')
       ax1.contour(lon_salt-360,lat_salt,salt_rg, [35.399, 35.4], colors='m', linewidths=0.6, projection=ccrs.PlateCarree())
       Q = ax1.quiver(lon_u[::scale_in]-360, lat_u[::scale_in], u_rg[::scale_in,::scale_in], v_rg[::scale_in,::scale_in], transform=ccrs.PlateCarree(), scale=v_scale, headwidth=7, headlength=7, width=0.004)
       ax1.quiverkey(Q, 0.46, 0.49, 0.5, str(0.5) + ' m/s', labelpos='N', coordinates='figure', transform=ccrs.PlateCarree())
       ax1.stock_img()
       ax1.coastlines('110m')
       ax1.add_feature(cfeature.RIVERS)
       ax1.plot(lon_obidos, lat_obidos, 'm*', transform=ccrs.PlateCarree())
       ax1.text(-69, -2, 'Obidos Gauge Station', fontsize=font_size, transform=ccrs.PlateCarree(),color='m')
       ax1.plot(lon_bolivar, lat_bolivar, 'b*', transform=ccrs.PlateCarree())
       ax1.text(-76, 6, 'Ciudad Bolivar Gauge Station', fontsize=font_size, transform=ccrs.PlateCarree(), color='b')
       mask_plot = river_mask.copy()
       mask_plot[np.isnan(river_mask)==True] = 0.
       ax1.contour(np.linspace(0,360,721), np.linspace(-90,84,348), mask_plot, [0, 1], colors='k', linewidths=1, transform=ccrs.PlateCarree())
#       ax1.contour(lon_sss, lat_sss-1, so_mask, [0, 1], colors='k', linewidths=1, transform=ccrs.PlateCarree())
       ax1.plot([-49, -49], [-0, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.plot([-65, -49], [15, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.plot([-65, -65], [10, 15], 'r-', linewidth=1, transform=ccrs.PlateCarree())
       ax1.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
       ax1.set_yticks([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20], crs=ccrs.PlateCarree())
       lon_formatter = LongitudeFormatter(zero_direction_label=True)
       lat_formatter = LatitudeFormatter()
       ax1.xaxis.set_major_formatter(lon_formatter)
       ax1.yaxis.set_major_formatter(lat_formatter)
       ax1.set_aspect('auto')
       plt.title('(a) Geography of Amazon River Basin and Amazon Plume Region', fontsize=font_size)

       cbaxes = fig.add_axes([0.08, 0.08, 0.42, 0.01])
       cbar = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[34.5, 35, 35.5, 36, 36.5, 37])
       cbar.set_clim(34.5,37)
       cbar.set_label('(PSU)',fontsize=font_size-1)
       setp(cbar.ax.yaxis.get_ticklabels(), fontsize=font_size-1)


# Figure 1(b)
       ax1 = plt.axes([0.58, 0.7, 0.35, 0.2])
       ax1.plot(tt, np.nanmean(prcp_gpcp_monthly, axis=0), 'ro-', label='GPCP', markersize=3)
       bar_plot_here(ax1, tt, prcp_gpcp_monthly, 'r')
       ax1.plot(tt, np.nanmean(prcp_precl_monthly, axis=0), 'mo-', label='PREC/L', markersize=3)
       bar_plot_here(ax1, tt, prcp_precl_monthly, 'm')
       ax1.plot(tt, np.nanmean(prcp_gpcc_monthly, axis=0), 'yo-', label='GPCC', markersize=3)
       bar_plot_here(ax1, tt, prcp_gpcc_monthly, 'y') 
       ax1.plot(tt, np.nanmean(trmm_monthly, axis=0), 'go-', label='TRMM', markersize=3)
       bar_plot_here(ax1, tt, trmm_monthly, 'g')
       ax1.plot(tt, np.nanmean(cmap_monthly, axis=0), 'co-', label='CMAP', markersize=3)
       bar_plot_here(ax1, tt, cmap_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(mm/day)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='upper right', ncol=3, fontsize=7)
       plt.title('(b) Precipitation Climatology', fontsize=font_size)

# Figure 1(c)
       ax1 = plt.axes([0.58, 0.42, 0.35, 0.2])
       ax1.plot(tt, np.nanmean(river3_monthly, axis=0), 'co-', label='Dai-Trenberth', markersize=3)
       ax1.plot(tt, np.nanmean(river_monthly, axis=0), 'bo-', label='HYBAM', markersize=3)
       bar_plot_here(ax1, tt, river_monthly, 'b')
       bar_plot_here(ax1, tt, river3_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(m$^3$/s x10$^5$)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='upper right', ncol=1, fontsize=7)
       plt.title('(c) Amazon River Discharge Climatology', fontsize=font_size)       

# Figure 1(d)
       ax1 = plt.axes([0.58, 0.14, 0.35, 0.2])
#       ax1.plot(tt, np.nanmean(sss_aquarius_monthly, axis=0), 'ko-', label='Aquarius', markersize=3)
#       bar_plot_here(ax1, tt, sss_aquarius_monthly, 'k')
#       ax1.plot(tt, np.nanmean(sss_smos_monthly, axis=0), 'bo-', label='SMOS', markersize=3)
#       bar_plot_here(ax1, tt, sss_smos_monthly, 'b')
       ax1.plot(tt, np.nanmean(sss_ecco_monthly, axis=0), 'go-', label='ECCO', markersize=3)
       bar_plot_here(ax1, tt, sss_ecco_monthly, 'g')
       ax1.plot(tt, np.nanmean(sss_monthly, axis=0), 'bo-', label='GECCO2', markersize=3)
       bar_plot_here(ax1, tt, sss_monthly, 'b')
       ax1.plot(tt, np.nanmean(sss_oras5_monthly, axis=0), 'yo-', label='ORAS5', markersize=3)
       bar_plot_here(ax1, tt, sss_oras5_monthly, 'y')
       ax1.plot(tt, np.nanmean(sss_soda_monthly, axis=0), 'mo-', label='SODA3.1.1', markersize=3)
       bar_plot_here(ax1, tt, sss_soda_monthly, 'm')
       ax1.plot(tt, np.nanmean(sss_en4_monthly, axis=0), 'co-', label='EN4', markersize=3)
       bar_plot_here(ax1, tt, sss_en4_monthly, 'c')
       ax1.grid('on')
       plt.ylabel('(PSU)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='lower right', ncol=1, fontsize=6)
       plt.title('(d) 5-meter Ocean Salinity Climatology', fontsize=font_size)
       plt.savefig('geographic_tmp_plot.jpg', format='jpeg', dpi=200)

       fig = plt.figure()
       ax1 = plt.axes([0.3, 0.3, 0.5, 0.5]) 
       ax1.plot(tt, np.nanmean(sss_aquarius_monthly, axis=0), 'ko-', label='Aquarius', markersize=3)
       bar_plot_here(ax1, tt, sss_aquarius_monthly, 'k')
       ax1.plot(tt, np.nanmean(sss_smos_monthly, axis=0), 'bo-', label='SMOS', markersize=3)
       bar_plot_here(ax1, tt, sss_smos_monthly, 'b')
#       ax1.plot(tt, np.nanmean(sss_ecco_monthly, axis=0), 'co-', label='ECCO', markersize=3)
#       bar_plot_here(ax1, tt, sss_ecco_monthly, 'c')
#       ax1.plot(tt, np.nanmean(sss_monthly, axis=0), 'go-', label='GECCO2', markersize=3)
#       bar_plot_here(ax1, tt, sss_monthly, 'g')
       ax1.grid('on')
       plt.ylabel('(PSU)', fontsize=10)
       plt.xticks(tt, month_txt, fontsize=9)
       plt.legend(loc='lower right', ncol=1, fontsize=6)
       plt.title('5-meter Ocean Salinity Climatology', fontsize=font_size)
       plt.savefig('suppl_geographic_tmp_plot.jpg', format='jpeg', dpi=200)

       plt.show()

#    else:

#       print('hahahah')

# ================================================================
# Describe this script and control flags
# ================================================================


