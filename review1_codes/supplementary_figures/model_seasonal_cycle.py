run_flag = 0

import numpy as np
import matplotlib.pyplot as plt
#import xarray as xr
import cartopy.crs as ccrs
import ast, os, sys
from pylab import setp
from netCDF4 import Dataset
from scipy.io import loadmat
from scipy.interpolate import griddata
from scipy import stats
from numpy import unravel_index

sys.path.append('/home/yliang/lib/python_tools/python_functions/data_process/')
import data_process_f
sys.path.append('/home/yliang/lib/python_tools/python_functions/statistics')
import statistical_f

def read_here1(dirname, filename):

    print(filename)
    f = Dataset(dirname + filename, mode='r')
    ts_max = f.variables['ts_max'][:].data
    ts_min = f.variables['ts_min'][:].data
    f.close()

    return ts_max, ts_min

def max_min_select(ts_monthly,year_N):

    print('max min select')
    ts_max = np.zeros((year_N))*np.nan
    ts_min = np.zeros((year_N))*np.nan
    for II in range(year_N-2):
        ts_max[II+1] = ts_monthly[II+1,:].max()
        ts_min[II+1] = ts_monthly[II+1,:].min()

    return ts_max, ts_min

def read_river_mask():

    print('read river mask')
    river_mask_tmp = np.loadtxt('/home/yliang/whoi_projects/Amazon_river/mask_new/obidosglobalnew.txt')
    river_mask_tmp = np.flipud(river_mask_tmp)

    river_mask = np.zeros((river_mask_tmp.shape))
    river_mask[:,0:361] = river_mask_tmp[:,360:].copy()
    river_mask[:,361:] = river_mask_tmp[:,:360].copy()
    river_mask[river_mask<=0.] = np.nan

    return river_mask

def sss_calculate_here(TLONG,TLAT,sss_test):

# Regrid
    lon = np.linspace(0,360, 360)
    lat = np.linspace(-79,89,180)
    [gridx_out, gridy_out] = np.meshgrid(lon, lat)
    ny = 180
    nx = 360
    [gridx_in, gridy_in] = [TLONG, TLAT]
    nt = sss_test.shape[0]

    sss_test_interp = np.zeros((nt,ny,nx))
    for NT in range(nt):
        print('Ocean Salinity')
        print(NT)
        sss_test_interp[NT,:,:] = griddata((gridy_in.ravel(), gridx_in.ravel()), np.squeeze(sss_test[NT,:,:]).ravel(), (gridy_out, gridx_out), method='nearest')

# Select region
    lon_amz1 = 295
    lon_amz2 = 311
    lat_amz1 = -3
    lat_amz2 = 15

    mask_sss = (sss_test_interp[11,:,:]/sss_test_interp[11,:,:]).copy()

    [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(lon_amz1,lon_amz2,lat_amz1,lat_amz2,lon,lat)
    area_so = data_process_f.area_calculate_nonuniform(lon, lat)
    ts_sss_test = np.zeros((nt))
    for NT in range(nt):
        ts_sss_test[NT] = np.nansum(sss_test_interp[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])

    return ts_sss_test

def rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,year_ref):

    ts_rain_test = np.zeros((year_N*12))
    ts_runoff_test = np.zeros((year_N*12))
    NN = 0
    for NY in range(year_N):
        for NM in range(12):
            filename = 'P_R_' + str(year_ref+NY+30) + '-' + str(NM+1).zfill(2) + '.nc'
            print(filename)
            f = Dataset(dirname + filename, mode='r')
            P_in = f.variables['RAIN'][:,:,:].data.squeeze() + f.variables['SNOW'][:,:,:].data.squeeze()
            P_in[P_in>1.e11] = np.nan
            R_in = f.variables['QCHANR'][:,176,248].data.squeeze()
#            R_in[R_in>1.e11] = np.nan
            f.close()
            ts_rain_test[NN] = np.nansum(P_in[:,:]*mask_rg*area_in)/np.nansum(mask_rg*area_in)
#            ts_runoff_test[NN] = np.nansum(R_in[:,:]*mask_rg*area_in)#/np.nansum(mask_rg*area_in)
            ts_runoff_test[NN] = R_in
            NN = NN + 1

    return ts_rain_test, ts_runoff_test

if run_flag == 1:
# ================================================================
# Basic setting
# ================================================================
   sig_level = 0.05

# ================================================================
# Read river forcing
# ================================================================
   year_N = 64
   tt = np.linspace(1948,2011,year_N)
   N = 3

# ================================================================
# Read simulation results
# ================================================================
#   ind_yr = 31 # 1979
   ind_yr = 30 # 1978  
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/'
   filename = 'GIAF_control_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control1 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   TLONG = f.variables['TLONG'][:,:].copy()
   TLAT = f.variables['TLAT'][:,:].copy()
   f.close()
   sss_control1[sss_control1>1.e30] = np.nan

   filename = 'GIAF_control_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control2 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   f.close()
   sss_control2[sss_control2>1.e30] = np.nan

   filename = 'GIAF_control_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control3 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   f.close()
   sss_control3[sss_control3>1.e30] = np.nan

   filename = 'GIAF_control_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control4 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   f.close()
   sss_control4[sss_control4>1.e30] = np.nan

   filename = 'GIAF_control_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control5 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   f.close()
   sss_control5[sss_control5>1.e30] = np.nan

   filename = 'GIAF_exp_1p75_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_1 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp3_1[sss_exp3_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_2 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp3_2[sss_exp3_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_3 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp3_3[sss_exp3_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_4 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp3_4[sss_exp3_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_5 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp3_5[sss_exp3_5>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_1 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp2_1[sss_exp2_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_2 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp2_2[sss_exp2_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_3 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp2_3[sss_exp2_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_4 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp2_4[sss_exp2_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_5 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_exp2_5[sss_exp2_5>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_1 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_spin_1[sss_spin_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_2 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_spin_2[sss_spin_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_3 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_spin_3[sss_spin_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_4 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_spin_4[sss_spin_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_5 = f.variables['SALT'][ind_yr*12:,0,:,:].data.copy()
   sss_spin_5[sss_spin_5>1.e30] = np.nan
   f.close()

   ts_sss_control1 = sss_calculate_here(TLONG,TLAT,sss_control1)
   ts_sss_control2 = sss_calculate_here(TLONG,TLAT,sss_control2)
   ts_sss_control3 = sss_calculate_here(TLONG,TLAT,sss_control3)
   ts_sss_control4 = sss_calculate_here(TLONG,TLAT,sss_control4)
   ts_sss_control5 = sss_calculate_here(TLONG,TLAT,sss_control5)

#   ts_sss_exp3_1 = sss_calculate_here(TLONG,TLAT,sss_exp3_1)
#   ts_sss_exp3_2 = sss_calculate_here(TLONG,TLAT,sss_exp3_2)
#   ts_sss_exp3_3 = sss_calculate_here(TLONG,TLAT,sss_exp3_3)
#   ts_sss_exp3_4 = sss_calculate_here(TLONG,TLAT,sss_exp3_4)
#   ts_sss_exp3_5 = sss_calculate_here(TLONG,TLAT,sss_exp3_5)

#   ts_sss_exp2_1 = sss_calculate_here(TLONG,TLAT,sss_exp2_1)
#   ts_sss_exp2_2 = sss_calculate_here(TLONG,TLAT,sss_exp2_2)
#   ts_sss_exp2_3 = sss_calculate_here(TLONG,TLAT,sss_exp2_3)
#   ts_sss_exp2_4 = sss_calculate_here(TLONG,TLAT,sss_exp2_4)
#   ts_sss_exp2_5 = sss_calculate_here(TLONG,TLAT,sss_exp2_5)

#   ts_sss_spin_1 = sss_calculate_here(TLONG,TLAT,sss_spin_1)
#   ts_sss_spin_2 = sss_calculate_here(TLONG,TLAT,sss_spin_2)
#   ts_sss_spin_3 = sss_calculate_here(TLONG,TLAT,sss_spin_3)
#   ts_sss_spin_4 = sss_calculate_here(TLONG,TLAT,sss_spin_4)
#   ts_sss_spin_5 = sss_calculate_here(TLONG,TLAT,sss_spin_5)

# ================================================================
# Read precipitation and river runoff
# ================================================================
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle1/'
   filename = 'P_R_1970-02.nc'
   f = Dataset(dirname + filename, mode='r')
   lon_in = f.variables['lon'][:].data
   lat_in = f.variables['lat'][:].data
   R_temp = f.variables['QCHANR'][:,:,:].data.squeeze()
   f.close()

   area_in = data_process_f.area_calculate_nonuniform(lon_in, lat_in)

# Read amazon river basin
# Interpolate to model grid
   river_mask = read_river_mask()
   lon_mask = np.linspace(0,360,721)
   lat_mask = np.linspace(-90,84,348)
   [grid_x, grid_y] = np.meshgrid(lon_mask, lat_mask)
   [grid_x_in, grid_y_in] = np.meshgrid(lon_in, lat_in) 
   mask_rg = griddata((grid_y.ravel(), grid_x.ravel()), np.squeeze(river_mask).ravel(), (grid_y_in, grid_x_in), method='nearest')
#   mask_rg[mask_rg!=1] = np.nan

   year_N = int(32)

   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle1/'
   [ts_rain_control1, ts_runoff_control1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle2/'
   [ts_rain_control2, ts_runoff_control2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle3/'
   [ts_rain_control3, ts_runoff_control3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle4/'
   [ts_rain_control4, ts_runoff_control4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle5/'
   [ts_rain_control5, ts_runoff_control5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle1/'

#   [ts_rain_exp2_1, ts_runoff_exp2_1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle2/'
#   [ts_rain_exp2_2, ts_runoff_exp2_2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle3/'
#   [ts_rain_exp2_3, ts_runoff_exp2_3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle4/'
#   [ts_rain_exp2_4, ts_runoff_exp2_4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle5/'
#   [ts_rain_exp2_5, ts_runoff_exp2_5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle1/'
#   [ts_rain_exp3_1, ts_runoff_exp3_1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle2/'
#   [ts_rain_exp3_2, ts_runoff_exp3_2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle3/'
#   [ts_rain_exp3_3, ts_runoff_exp3_3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle4/'
#   [ts_rain_exp3_4, ts_runoff_exp3_4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
#   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle5/'
#   [ts_rain_exp3_5, ts_runoff_exp3_5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)

# ================================================================
# Plot figures
# ================================================================
if True:

   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   text_txt = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
   ttt = np.linspace(1,12,12)

   plt.axes([0.57, 0.37, 0.4, 0.35])
   ts_test = np.nanmean(ts_sss_control1.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bs-', label='cycle1')
   ts_test = np.nanmean(ts_sss_control2.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'b^-', label='cycle2')
   ts_test = np.nanmean(ts_sss_control3.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'b*-', label='cycle3')
   ts_test = np.nanmean(ts_sss_control4.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bv-', label='cycle4')
   ts_test = np.nanmean(ts_sss_control5.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bD-', label='cycle5')

   plt.ylabel('(PSU)')
   plt.legend()
   plt.xticks(ttt, text_txt, fontsize=10)
   plt.title('(b) 5-m ocean salinity seasonal cycle')

   plt.axes([0.07, 0.37, 0.4, 0.35])
   ts_test = np.nanmean(ts_runoff_control1.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bs-', label='cycle1')
   ts_test = np.nanmean(ts_runoff_control2.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'b^-', label='cycle2')
   ts_test = np.nanmean(ts_runoff_control3.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'b*-', label='cycle3')
   ts_test = np.nanmean(ts_runoff_control4.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bv-', label='cycle4')
   ts_test = np.nanmean(ts_runoff_control5.reshape((32,12)), axis=0)
   plt.plot(ttt,ts_test,'bD-', label='cycle5')

   plt.ylabel('(m$^{3}$/s)')
   plt.legend()
   plt.xticks(ttt, text_txt, fontsize=10)
   plt.title('(a) river discharge seasonal cycle')

   plt.savefig('model_seasonal_cycle_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()

   tt2 = np.linspace(1978,2009,32)

   plt.axes([0.1, 0.57, 0.4, 0.35])
   plt.plot(tt2[:], rain_range_control1*86400, 'k-', label='control')
   plt.plot(tt2[:], rain_range_exp2_1*86400, 'm-', label='x1.5')
   plt.plot(tt2[:], rain_range_exp3_1*86400, 'r-', label='x1.75')
   plt.legend(ncol=3)
   plt.title('(a) Amazon precipitation seasonality forcings')
   plt.ylabel('(mm/s)')
   plt.axis([1979, 2007, 2, 10])

   plt.axes([0.6, 0.57, 0.37, 0.35])
   runoff_control_trend = np.array([runoff_control1_trend, runoff_control2_trend, runoff_control3_trend, runoff_control4_trend, runoff_control5_trend])*10
#   runoff_control_range = runoff_control_trend.std()
   runoff_exp2_trend = np.array([runoff_exp2_1_trend, runoff_exp2_2_trend, runoff_exp2_3_trend, runoff_exp2_4_trend, runoff_exp2_5_trend])*10
#   runoff_exp2_range = runoff_exp2_trend.std()
   runoff_exp3_trend = np.array([runoff_exp3_1_trend, runoff_exp3_2_trend, runoff_exp3_3_trend, runoff_exp3_4_trend, runoff_exp3_5_trend])*10
#   runoff_exp3_range = runoff_exp3_trend.std()
   plt.plot([1, 1.5, 1.75], np.array([runoff_control_trend.mean(), runoff_exp2_trend.mean(), runoff_exp3_trend.mean()]), 'b-', markerfacecolor='none', markersize=6)
   plt.plot([1, 1.5, 1.75], np.array([runoff_control_trend.mean(), runoff_exp2_trend.mean(), runoff_exp3_trend.mean()]), 'bo', markerfacecolor='none', markersize=6, label='average')
   plt.plot([1, 1,], [runoff_control1_trend*10, runoff_control1_trend*10], 'bs', markersize=6, markerfacecolor='none', label='cycle1')
   plt.plot([1, 1,], [runoff_control2_trend*10, runoff_control2_trend*10], 'b^', markersize=6, markerfacecolor='none', label='cycle2')
   plt.plot([1, 1,], [runoff_control3_trend*10, runoff_control3_trend*10], 'b*', markersize=6, markerfacecolor='none', label='cycle3')
   plt.plot([1, 1,], [runoff_control4_trend*10, runoff_control4_trend*10], 'bv', markersize=6, markerfacecolor='none', label='cycle4')
   plt.plot([1, 1,], [runoff_control5_trend*10, runoff_control5_trend*10], 'bD', markersize=6, markerfacecolor='none', label='cycle5')

   plt.plot([1.5, 1.5], [runoff_exp2_1_trend*10, runoff_exp2_1_trend*10], 'bs', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [runoff_exp2_2_trend*10, runoff_exp2_2_trend*10], 'b^', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [runoff_exp2_3_trend*10, runoff_exp2_3_trend*10], 'b*', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [runoff_exp2_4_trend*10, runoff_exp2_4_trend*10], 'bv', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [runoff_exp2_5_trend*10, runoff_exp2_5_trend*10], 'bD', markersize=6, markerfacecolor='none')

   plt.plot([1.75, 1.75], [runoff_exp3_1_trend*10, runoff_exp3_1_trend*10], 'bs', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [runoff_exp3_2_trend*10, runoff_exp3_2_trend*10], 'b^', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [runoff_exp3_3_trend*10, runoff_exp3_3_trend*10], 'b*', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [runoff_exp3_4_trend*10, runoff_exp3_4_trend*10], 'bv', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [runoff_exp3_5_trend*10, runoff_exp3_5_trend*10], 'bD', markersize=6, markerfacecolor='none')

   plt.axis([0.75, 2, 14000, 30000])
   plt.xticks([1, 1.5, 1.75], ['Control', 'x1.5', 'x1.75'])
   plt.ylabel('(m$^3$/s per decade)')
   plt.grid('on')
   plt.legend(loc='lower right')
   plt.title('(b) Amazon discharge seasonality sensitivity tests')

   plt.axes([0.1, 0.14, 0.4, 0.35])
   plt.plot(tt, river_max_control-river_min_control, 'k-', label='control')
#   plt.plot(tt, river_max_exp1-river_min_exp1, 'b-', label='x0.5')
   plt.plot(tt, river_max_exp2-river_min_exp2, 'm-', label='x1.5')
   plt.plot(tt, river_max_exp3-river_min_exp3, 'r-', label='x1.75')
   plt.legend(ncol=4)
   plt.title('(c) Amazon runoff seasonality forcings')
   plt.ylabel('(m$^3$/s)')
   plt.axis([1979, 2007, 50000, 375000])
    
   plt.axes([0.6, 0.14, 0.37, 0.35])
   plt.plot([1, 1.5, 1.75], [sss_control_mean*10, sss_exp2_mean*10, sss_exp3_mean*10], 'b-', markerfacecolor='none')
   plt.plot([1, 1.5, 1.75], [sss_control_mean*10, sss_exp2_mean*10, sss_exp3_mean*10], 'bo', markerfacecolor='none', label='average')

   plt.plot([1, 1,], [sss_control1_trend*10, sss_control1_trend*10], 'bs', markersize=6, markerfacecolor='none', label='cycle1')
   plt.plot([1, 1,], [sss_control2_trend*10, sss_control2_trend*10], 'b^', markersize=6, markerfacecolor='none', label='cycle2')
   plt.plot([1, 1,], [sss_control3_trend*10, sss_control3_trend*10], 'b*', markersize=6, markerfacecolor='none', label='cycle3')
   plt.plot([1, 1,], [sss_control4_trend*10, sss_control4_trend*10], 'bv', markersize=6, markerfacecolor='none', label='cycle4')
   plt.plot([1, 1,], [sss_control5_trend*10, sss_control5_trend*10], 'bD', markersize=6, markerfacecolor='none', label='cycle5')

   plt.plot([1.5, 1.5], [sss_exp2_1_trend*10, sss_exp2_1_trend*10], 'bs', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [sss_exp2_2_trend*10, sss_exp2_2_trend*10], 'b^', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [sss_exp2_3_trend*10, sss_exp2_3_trend*10], 'b*', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [sss_exp2_4_trend*10, sss_exp2_4_trend*10], 'bv', markersize=6, markerfacecolor='none')
   plt.plot([1.5, 1.5], [sss_exp2_5_trend*10, sss_exp2_5_trend*10], 'bD', markersize=6, markerfacecolor='none')

   plt.plot([1.75, 1.75], [sss_exp3_1_trend*10, sss_exp3_1_trend*10], 'bs', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [sss_exp3_2_trend*10, sss_exp3_2_trend*10], 'b^', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [sss_exp3_3_trend*10, sss_exp3_3_trend*10], 'b*', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [sss_exp3_4_trend*10, sss_exp3_4_trend*10], 'bv', markersize=6, markerfacecolor='none')
   plt.plot([1.75, 1.75], [sss_exp3_5_trend*10, sss_exp3_5_trend*10], 'bD', markersize=6, markerfacecolor='none')

#   plt.plot([1, 1], [sss_core2_trends[-2]*10, sss_core2_trends[-2]*10], 'm*', label='CORE-II NCAR')
   plt.plot([1, 1], [sss_spin_mean*10, sss_spin_mean*10], 'k*', label='spin-up')
   plt.axis([0.75, 2, 0.005, 0.035])
#   plt.axis([0.75, 2, 0.01, 0.04])
   plt.xticks([1, 1.5, 1.75], ['Control', 'x1.5', 'x1.75'])
   plt.ylabel('(PSU per decade)')
   plt.legend(loc='lower right')
   plt.grid('on')
   plt.title('(d) Ocean salinity seasonality sensitivity tests')

   plt.savefig('core2_tmp_plot.jpg', format='jpeg', dpi=200)
   plt.show()

   sys.exit()

   plt.axes([0.59, 0.15, 0.4, 0.35])
   tt = np.linspace(1,len(core2_case_txt)+1,len(core2_case_txt)+1)
   plt.bar(tt, sss_core2_trends, width=0.5, color='m', alpha=0.5)
   plt.plot([tt[-1], tt[-1]], [sss_core2_trends[-1], sss_core2_trends[-1]+np.nanstd(sss_core2_trends[:-1])], 'k-')
   plt.plot([tt[-1], tt[-1]], [sss_core2_trends[-1], sss_core2_trends[-1]-np.nanstd(sss_core2_trends[:-1])], 'k-')
   for II in range(len(sss_core2_sig)):
       if sss_core2_sig[II] == 1:
          plt.plot([tt[II], tt[II]], [sss_core2_trends[II], sss_core2_trends[II]], 'k*')
   plt.bar([tt[-1]+1, tt[-1]+1], [sss_control_trend, sss_control_trend], width=0.5, color='m', alpha=0.5)
   plt.bar([tt[-1]+2, tt[-1]+2], [sss_exp3_trend, sss_exp3_trend], width=0.5, color='m', alpha=0.5)

   plt.ylabel('(PSU/year)')
   plt.xticks(np.linspace(1,len(core2_case_txt)+3,len(core2_case_txt)+3), ['BERGEN','CERFAC','CMCC','CNRM','FSU','GFDL-GOLD','GFDL-MOM','ICTP','INMOM','MRI-A','MRI-F','NCAR','Mean', 'Control', 'x1.75'], rotation=60) 
   plt.title('(d) CORE-II Models')

   sys.exit()

