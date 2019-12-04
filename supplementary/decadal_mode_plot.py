flag_run = 1
# ================================================================
# Yu-Chiao @ WHOI August 21, 2018
# Investigate snow variability and trend
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
#from mpl_toolkits.basemap import Basemap
from IPython import get_ipython
import sys, os
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from pylab import setp, genfromtxt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from datetime import datetime
from sklearn import linear_model
import xlrd
from scipy.io import loadmat
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

sys.path.append('/home/yliang/lib/python_tools/python_functions/whoi_projects/')
import whoi_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f, ORAS4_data_process_f, HadSST_data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/statistics')
import statistical_f, MCA_f, EOF_f
sys.path.append('/home/yliang/lib/python_tools/python_functions')
import plot_functions

# ================================================================
# Define functions
# ================================================================

# ================================================================
# Set parameters
# ================================================================
year_st = 1979
year_ed = 2015
yearN = year_ed - year_st + 1
year = np.linspace(year_st,year_ed,yearN)
year_ref = 1958
t0 = (year_st-year_ref)*12
t1 = (year_ed-year_ref)*12+12

lon_amz1 = 309
lon_amz2 = 310
lat_amz1 = -2
lat_amz2 = 2

sig_level = 0.05
q_fdr = sig_level*2.

layer = 0

if flag_run == 1:
# ================================================================
# Read AMO index
# ================================================================
   year_st = 1948
   year_ed = 2013
   year_ref = 1856
   t0 = (year_st-year_ref)
   t1 = (year_ed-year_ref)+1
   dirname = '/stormtrack/data4/yliang/Indices/'
   filename = 'amon.us.long.data'
   f = open(dirname+filename, "r")
   data_read = genfromtxt(dirname + filename, dtype='float', skip_header=0)[t0:t1,1:]
   f.close()
   amo_index = np.nanmean(data_read, axis=1)
   amo_index = mlab.detrend_linear(amo_index)
   N = 11
   amo_rmean = amo_index.copy()*np.nan
   amo_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(amo_index, np.ones((N,))/N, mode='valid')

# ================================================================
# Read PDO index
# ================================================================
   year_st = 1948
   year_ed = 2013
   year_ref = 1900
   t0 = (year_st-year_ref)
   t1 = (year_ed-year_ref)+1
   dirname = '/stormtrack/data4/yliang/Indices/'
   filename = 'PDO.latest.txt'
   f = open(dirname+filename, "r")
   data_read = genfromtxt(dirname + filename, dtype='float', skip_header=0)[t0:t1,1:]
   f.close()
   pdo_index = np.nanmean(data_read, axis=1)
   pdo_index = mlab.detrend_linear(pdo_index)
   N = 11
   pdo_rmean = pdo_index.copy()*np.nan
   pdo_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(pdo_index, np.ones((N,))/N, mode='valid')

# ================================================================
# Read IPO index
# ================================================================
   year_st = 1948
   year_ed = 2013
   year_ref = 1854
   t0 = (year_st-year_ref)
   t1 = (year_ed-year_ref)+1
   dirname = '/stormtrack/data4/yliang/Indices/'
   filename = 'tpi.timeseries.ersstv5.data'
   f = open(dirname+filename, "r")
   data_read = genfromtxt(dirname + filename, dtype='float', skip_header=0)[t0:t1,1:]
   f.close()
   ipo_index = np.nanmean(data_read, axis=1)
   ipo_index = mlab.detrend_linear(ipo_index)
   N = 11
   ipo_rmean = ipo_index.copy()*np.nan
   ipo_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ipo_index, np.ones((N,))/N, mode='valid')

# ================================================================
# main begins
# ================================================================
#  Read river mask
   x = loadmat('/home/yliang/whoi_projects/Amazon_river/observations/precipitation/river_basin_mask.mat')
   river_tmp = x['A'].T
   river_mask = np.zeros((river_tmp.shape))
   river_mask[:,0:180] = river_tmp[:,180:]
   river_mask[:,180:] = river_tmp[:,0:180]
   del river_tmp

# ================================================================
# Read GPCP precipitation
# ================================================================
   year_st = 1979
   year_ed = 2015 
   yearN = year_ed - year_st + 1
   year_gpcp = np.linspace(year_st,year_ed,yearN)
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
   [grid_x, grid_y] = np.meshgrid(np.linspace(1,360,360), np.linspace(-90,90,180))
   [lon_x, lat_y] = np.meshgrid(lon, lat)
   mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest') 
   mask_rg[mask_rg!=1] = np.nan

#   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(270,330,-20,0,lon,lat)

   area = data_process_f.area_calculate(lon, lat)

   nt = var_tmp.shape[0]
   ts_gpcp = np.zeros((nt))

   for NT in range(nt):
       ts_gpcp[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg) /np.nansum(area*mask_rg)

# ================================================================
# Read PREC/L precipitation
# ================================================================
   year_st = 1948
   year_ed = 2014
   yearN = year_ed - year_st + 1
   year_precl = np.linspace(year_st,year_ed,yearN)
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
   [grid_x, grid_y] = np.meshgrid(np.linspace(1,360,360), np.linspace(-90,90,180))
   [lon_x, lat_y] = np.meshgrid(lon, lat)
   mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
   mask_rg[mask_rg!=1] = np.nan

#   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(270,330,-20,0,lon,lat)

   area = data_process_f.area_calculate(lon, lat)

   nt = var_tmp.shape[0]
   ts_precl = np.zeros((nt))

   for NT in range(nt):
       ts_precl[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg)/np.nansum(area*mask_rg)

# ================================================================
# Read GPCC precipitation
# ================================================================
   year_st = 1948
   year_ed = 2010
   yearN = year_ed - year_st + 1
   year_gpcc = np.linspace(year_st,year_ed,yearN)
   year_ref = 1901
   t0 = (year_st-year_ref)*12
   t1 = (year_ed-year_ref)*12+12

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
   [grid_x, grid_y] = np.meshgrid(np.linspace(1,360,360), np.linspace(-90,90,180))
   [lon_x, lat_y] = np.meshgrid(lon, lat)
   mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (lat_y, lon_x), method='nearest')
   mask_rg[mask_rg!=1] = np.nan

#   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(270,330,-20,0,lon,lat)

   area = data_process_f.area_calculate(lon, lat)

   nt = var_tmp.shape[0]
   ts_gpcc = np.zeros((nt))

   for NT in range(nt):
       ts_gpcc[NT] = np.nansum(var_tmp[NT,:,:]*area*mask_rg)/np.nansum(area*mask_rg)

# ================================================================
# Read Hybam river discharge data
# ================================================================
   dirname = '/home/yliang/whoi_projects/Amazon_river/observations/streamflow/'
   filename = 'Amazon_river_obidos_1969_2015.xls'

   year_st = 1969
   year_ed = 2015
   yearN = year_ed - year_st + 1
   year_river = np.linspace(year_st,year_ed,yearN)
   year_ref = 1969
   t0 = (year_st-year_ref)*12
   t1 = (year_ed-year_ref)*12+12

   wb = xlrd.open_workbook(dirname+filename)
   sheet = wb.sheet_by_index(0)

   river_discharge = np.zeros((yearN*12))

   for II in range(yearN*12):
       river_discharge[II] = sheet.cell_value(II+2,3)

# ================================================================
# Read Dai and Trenberth dataset
# ================================================================
   dirname = '/stormtrack/data4/yliang/observations/streamflow/Dai_Trenberth/'
   filename = 'coastal-stns-Vol-monthly.updated-Aug2014.nc'

   year_st = 1948
   year_ed = 2013
   yearN = year_ed - year_st + 1
   year_river_dai = np.linspace(year_st,year_ed,yearN)
   year_ref = 1900
   t0 = (year_st-year_ref)*12
   t1 = (year_ed-year_ref)*12+12

   f = Dataset(dirname + filename, mode='r')
   flow = f.variables['FLOW'][t0:t1,0].data
   flow[flow<0.] = np.nan
   river_name = f.variables['riv_name'][0].data
   time = f.variables['time'][t0:t1].data
   f.close()

# ================================================================
# Perform 3-month running average
# ================================================================
   N = 3
   ts_gpcp_rmean = ts_gpcp.copy()*np.nan
   ts_gpcp_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_gpcp, np.ones((N,))/N, mode='valid')
   ts_precl_rmean = ts_precl.copy()*np.nan
   ts_precl_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_precl, np.ones((N,))/N, mode='valid')
   ts_gpcc_rmean = ts_gpcc.copy()*np.nan
   ts_gpcc_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_gpcc, np.ones((N,))/N, mode='valid')
   river_discharge_rmean = river_discharge.copy()*np.nan
   river_discharge_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(river_discharge, np.ones((N,))/N, mode='valid')
   flow_rmean = flow.copy()*np.nan
   flow_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(flow, np.ones((N,))/N, mode='valid')

# ================================================================
# Arrange time series
# ================================================================
   river_monthly = river_discharge_rmean.reshape((len(year_river),12))
   river_dai_monthly = flow_rmean.reshape((len(year_river_dai),12))
   gpcp_monthly = ts_gpcp_rmean.reshape((len(year_gpcp),12))   
   gpcc_monthly = ts_gpcc_rmean.reshape((len(year_gpcc),12))
   precl_monthly = ts_precl_rmean.reshape((len(year_precl),12))

   river_max = np.zeros((len(year_river)))
   river_min = np.zeros((len(year_river)))

   river_dai_max = np.zeros((len(year_river_dai)))
   river_dai_min = np.zeros((len(year_river_dai)))

   gpcp_max = np.zeros((len(year_gpcp)))
   gpcp_min = np.zeros((len(year_gpcp)))

   gpcc_max = np.zeros((len(year_gpcc)))
   gpcc_min = np.zeros((len(year_gpcc)))

   precl_max = np.zeros((len(year_precl)))
   precl_min = np.zeros((len(year_precl)))

   for II in range(len(year_river)):
       river_max[II] = river_monthly[II,:].max()
       river_min[II] = river_monthly[II,:].min()

   for II in range(len(year_river_dai)):
       river_dai_max[II] = np.nanmax(river_dai_monthly[II,:])
       river_dai_min[II] = np.nanmin(river_dai_monthly[II,:])

   for II in range(len(year_gpcp)):
       gpcp_max[II] = gpcp_monthly[II,:].max()
       gpcp_min[II] = gpcp_monthly[II,:].min()

   for II in range(len(year_gpcc)):
       gpcc_max[II] = gpcc_monthly[II,:].max()
       gpcc_min[II] = gpcc_monthly[II,:].min()

   for II in range(len(year_precl)):
       precl_max[II] = precl_monthly[II,:].max()
       precl_min[II] = precl_monthly[II,:].min()

   river_clim = np.nanmean(river_monthly, axis=0)
   river_dai_clim = np.nanmean(river_dai_monthly, axis=0)
   gpcp_clim = np.nanmean(gpcp_monthly, axis=0)
   gpcc_clim = np.nanmean(gpcc_monthly, axis=0)
   precl_clim = np.nanmean(precl_monthly, axis=0)

# ================================================================
# Read field
# ================================================================
if True:

   font_size = 11

   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 7, forward=True)

   factor = 1./10000.

   N = 11

   ax1 = fig.add_axes([0.25, 0.6, 0.45, 0.35])
   ax1.plot(year_river_dai, (river_dai_max-river_dai_min)*factor, 'b-', linewidth=0.5, alpha=0.5)
   ax1.plot(year_river_dai[5:-5], (np.convolve(river_dai_max-river_dai_min, np.ones((N,))/N, mode='valid'))*factor, 'b-', linewidth=2)
   ax1.set_ylabel('discharge (m$^3$/s x 10000)', color='b', fontsize=font_size)
   ax1.tick_params('y', colors='b')
#   ax1.axis([1950, 2016, 70000*factor, 180000*factor])
#   plt.xticks([1975,1980,1985,1990,1995,2000,2005,2010,2015],['1975','1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=font_size, rotation=30)
   plt.xticks([1950,1960,1970,1980,1990,2000,2010],['1950','1960','1970','1980','1990','2000','2010'], fontsize=font_size, rotation=30)
   ax2 = ax1.twinx()
   ax2.plot(year_river_dai, amo_index, 'm-', linewidth=0.5)
   ax2.plot(year_river_dai, amo_rmean, 'm-', linewidth=2)
   ax2.set_ylabel('AMV Index', color='m', fontsize=font_size)
   ax2.tick_params('y', colors='m')
   ax2.axis([1950, 2016, -0.5, 0.5])
   plt.title('(a) AMV Index', fontsize=font_size)

# PDO
   ax1 = fig.add_axes([0.25, 0.15, 0.45, 0.35])
   ax1.plot(year_river_dai, (river_dai_max-river_dai_min)*factor, 'b-', linewidth=0.5, alpha=0.5)
   ax1.plot(year_river_dai[5:-5], (np.convolve(river_dai_max-river_dai_min, np.ones((N,))/N, mode='valid'))*factor, 'b-', linewidth=2)
   ax1.set_ylabel('discharge (m$^3$/s x 10000)', color='b', fontsize=font_size)
   ax1.tick_params('y', colors='b')
#   ax1.axis([1950, 2016, -0.5, 0.5])
#   plt.xticks([1975,1980,1985,1990,1995,2000,2005,2010,2015],['1975','1980','1985','1990','1995','2000','2005','2010','2015'], fontsize=font_size, rotation=30)
   plt.xticks([1950,1960,1970,1980,1990,2000,2010],['1950','1960','1970','1980','1990','2000','2010'], fontsize=font_size, rotation=30)
   ax2 = ax1.twinx()
   ax2.plot(year_river_dai, -pdo_index, 'm-', linewidth=0.5)
   ax2.plot(year_river_dai, -pdo_rmean, 'm-', linewidth=2)
   ax2.plot(year_river_dai, -ipo_index, 'y-', linewidth=0.5)
   ax2.plot(year_river_dai, -ipo_rmean, 'y-', linewidth=2)
   ax2.set_ylabel('PDV & IPV Index x -1', color='m', fontsize=font_size)
   ax2.tick_params('y', colors='m')
   ax2.axis([1950, 2016, -2, 2])
   plt.title('(b) PDV & IPV Indices', fontsize=font_size)

   plt.savefig('decadal_index_temp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

