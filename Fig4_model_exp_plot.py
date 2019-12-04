run_flag = 0

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
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
        ts_sss_test[NT] = np.nansum(sss_test_interp[NT,y1:y2+1,x1:x2+1]*area_so[y1:y2+1,x1:x2+1])/np.nansum(area_so[y1:y2+1,x1:x2+1]*mask_sss[y1:y2+1,x1:x2+1])

    ts_rmean = ts_sss_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_sss_test.copy(), np.ones((N,))/N, mode='valid')
    sss_monthly = ts_rmean.reshape((int(nt/12),12))
    [sss_max_test, sss_min_test] = max_min_select(sss_monthly,int(nt/12))

    [sss_test_trend, sss_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt/12,nt/12)[0:-2], sss_max_test[1:-1]-sss_min_test[1:-1], sig_level)

    return sss_test_trend, sss_test_sig, pvalue 

def rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,year_ref):

    ts_rain_test = np.zeros((year_N*12))
    ts_runoff_test = np.zeros((year_N*12))
    NN = 0
    for NY in range(year_N):
        for NM in range(12):
            filename = 'P_R_' + str(year_ref+NY+31) + '-' + str(NM+1).zfill(2) + '.nc'
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

    ts_rmean = ts_rain_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_rain_test.copy(), np.ones((N,))/N, mode='valid')
    rain_monthly = ts_rmean.reshape((year_N,12))
    [rain_max_test, rain_min_test] = max_min_select(rain_monthly,year_N)

    [rain_test_trend, rain_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,year_N,year_N)[0:-2], rain_max_test[1:-1]-rain_min_test[1:-1], sig_level)

    ts_rmean = ts_runoff_test.copy()
    ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_test.copy(), np.ones((N,))/N, mode='valid')
    runoff_monthly = ts_rmean.reshape((year_N,12))
    [runoff_max_test, runoff_min_test] = max_min_select(runoff_monthly,year_N)

    [runoff_test_trend, runoff_test_sig, pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,year_N,year_N)[0:-2], runoff_max_test[1:-1]-runoff_min_test[1:-1], sig_level)

    return rain_test_trend, rain_test_sig, runoff_test_trend, runoff_test_sig, rain_max_test-rain_min_test

if run_flag == 1:
# ================================================================
# Basic setting
# ================================================================
   sig_level = 0.05

   dirname = '/home/yliang/whoi_projects/Amazon_river/model/CORE-II/'
   argo_case_txt = ['argo_ciso_tmp.nc','argo_iprc_tmp.nc','argo_isas_tmp.nc','argo_jamstec_tmp.nc','argo_scripps_tmp.nc']
   cesm2_case_txt = ['cesm2_control_tmp.nc']
#   core2_case_txt = ['core2_awi_tmp.nc','core2_bergen_tmp.nc','core2_cerfacs_tmp.nc','core2_cmcc_tmp.nc','core2_cnrm_tmp.nc','core2_fsu_tmp.nc','core2_gfdlg_tmp.nc','core2_gfdlm_tmp.nc','core2_ictp_tmp.nc','core2_inmom_tmp.nc','core2_mria_tmp.nc','core2_mrif_tmp.nc','core2_ncar_tmp.nc','core2_nocs_tmp.nc']
   core2_case_txt = ['core2_bergen_tmp.nc','core2_cerfacs_tmp.nc','core2_cmcc_tmp.nc','core2_cnrm_tmp.nc','core2_fsu_tmp.nc','core2_gfdlg_tmp.nc','core2_gfdlm_tmp.nc','core2_ictp_tmp.nc','core2_inmom_tmp.nc','core2_mria_tmp.nc','core2_mrif_tmp.nc','core2_ncar_tmp.nc']


   sss_argo_trends = np.zeros((len(argo_case_txt)+1))
   sss_argo_sig = np.zeros((len(argo_case_txt)+1))
   sss_cesm2_trends = np.zeros((len(cesm2_case_txt)+1))
   sss_core2_trends = np.zeros((len(core2_case_txt)+1))
   sss_core2_sig = np.zeros((len(core2_case_txt)+1))

   ts_sss_core2 = np.zeros((len(core2_case_txt)+1,60))

   run_trend_core2 = np.zeros((len(core2_case_txt), 58)) 
   run_sig_core2 = np.zeros((len(core2_case_txt), 58))

# ================================================================
# Read Argo data
# ================================================================
   for II in range(len(argo_case_txt)):
       [ts_max, ts_min] = read_here1(dirname, argo_case_txt[II])
       nt = len(ts_max)
#       sss_argo_trends[II] = np.polyfit(np.linspace(1,nt,nt)[1:-1],ts_max[1:-1]-ts_min[1:-1],1)[0]
#       sss_argo_sig[II] = stats.ttest_ind(np.linspace(1,nt,nt)[1:-1]*0.,ts_max[1:-1]-ts_min[1:-1])[1]
       [sss_argo_trends[II], sss_argo_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[1:-1], ts_max[1:-1]-ts_min[1:-1], sig_level)
       if II == 0: ts_argo_ciso = ts_max - ts_min
       if II == 1: ts_argo_iprc = ts_max - ts_min
       if II == 2: ts_argo_isas = ts_max - ts_min
       if II == 3: ts_argo_jamstec = ts_max - ts_min
       if II == 4: ts_argo_scripps = ts_max - ts_min
   sss_argo_trends[-1] = sss_argo_trends[:-1].mean()

   for II in range(len(core2_case_txt)):
       [ts_max, ts_min] = read_here1(dirname, core2_case_txt[II])
       nt = len(ts_max[-31:])
#       nt = len(ts_max[:])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_max[-31:-1], sig_level)
       print('sss core2 max', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_min[-31:-1], sig_level)
       print('sss core2 min', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       [sss_core2_trends[II], sss_core2_sig[II], pvalue] = statistical_f.linear_regression_1D_t_test_with_sample_size_adjusted(np.linspace(1,nt,nt)[:-1], ts_max[-31:-1]-ts_min[-31:-1], sig_level)
       print('sss core2 max-min', str(sss_core2_trends[II]*10), sss_core2_sig[II])

       ts_sss_core2[II,:] = ts_max - ts_min

#       [y_slope, y_sig] = statistical_f.running_linear_trends(np.linspace(1,nt,nt)[1:-1], ts_max[1:-1]-ts_min[1:-1],21,sig_level)
#       run_trend_core2[II,:] = y_slope.copy()
#       run_sig_core2[II,:] = y_sig.copy()
                
   sss_core2_trends[-1] = np.nanmean(sss_core2_trends[:-1])
   ts_sss_core2[-1,:] = np.nanmean(ts_sss_core2[0:-1,:], axis=0)

# ================================================================
# Read river forcing
# ================================================================
   year_N = 64
   tt = np.linspace(1948,2011,year_N)
   N = 3

# Read control run
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/forcing/runoff/'
   filename = 'runoff.daitren.iaf.20120419.nc'
   f = Dataset(dirname + filename, mode='r')
   runoff_control = f.variables['runoff'][:,:,:].data
   runoff_area = f.variables['area'][:,:].data
   lat_max = unravel_index(np.nanmean(runoff_control[:,:,:], axis=0).argmax(), np.nanmean(runoff_control[:,:,:], axis=0).shape)[0]
   lon_max = unravel_index(np.nanmean(runoff_control[:,:,:], axis=0).argmax(), np.nanmean(runoff_control[:,:,:], axis=0).shape)[1]
   runoff_control[runoff_control==0.] = np.nan
   ts_runoff_tmp = runoff_control[:,lat_max,lon_max]*runoff_area[lat_max,lon_max].copy()/1000.
   f.close()

   ts_rmean = ts_runoff_tmp.copy()
   ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_tmp.copy(), np.ones((N,))/N, mode='valid')
   river_monthly = ts_rmean.reshape((year_N,12))
   [river_max_control, river_min_control] = max_min_select(river_monthly,year_N)

# Read x0.5 case
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/forcing/runoff/'
   filename = 'runoff_daitren_iaf_control_seasonx0p5.nc'
   f = Dataset(dirname + filename, mode='r')
   runoff_exp1 = f.variables['runoff'][:,:,:].data
   runoff_exp1[runoff_exp1==0.] = np.nan
   ts_runoff_tmp = runoff_exp1[:,lat_max,lon_max]*runoff_area[lat_max,lon_max].copy()/1000.
   f.close()

   ts_rmean = ts_runoff_tmp.copy()
   ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_tmp.copy(), np.ones((N,))/N, mode='valid')
   river_monthly = ts_rmean.reshape((year_N,12))
   [river_max_exp1, river_min_exp1] = max_min_select(river_monthly,year_N)

# Read x1.5 case
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/forcing/runoff/'
   filename = 'runoff_daitren_iaf_control_seasonx1p5.nc'
   f = Dataset(dirname + filename, mode='r')
   runoff_exp2 = f.variables['runoff'][:,:,:].data
   runoff_exp2[runoff_exp2==0.] = np.nan
   ts_runoff_tmp = runoff_exp2[:,lat_max,lon_max]*runoff_area[lat_max,lon_max].copy()/1000.
   f.close()

   ts_rmean = ts_runoff_tmp.copy()
   ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_tmp.copy(), np.ones((N,))/N, mode='valid')
   river_monthly = ts_rmean.reshape((year_N,12))
   [river_max_exp2, river_min_exp2] = max_min_select(river_monthly,year_N)

# Read x1.75 case
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/forcing/runoff/'
   filename = 'runoff_daitren_iaf_control_seasonx1p75.nc'
   f = Dataset(dirname + filename, mode='r')
   runoff_exp3 = f.variables['runoff'][:,:,:].data
   runoff_exp3[runoff_exp3==0.] = np.nan
   ts_runoff_tmp = runoff_exp3[:,lat_max,lon_max]*runoff_area[lat_max,lon_max].copy()/1000.
   f.close()

   ts_rmean = ts_runoff_tmp.copy()
   ts_rmean[int((N-1)/2):-int((N-1)/2)] = np.convolve(ts_runoff_tmp.copy(), np.ones((N,))/N, mode='valid')
   river_monthly = ts_rmean.reshape((year_N,12))
   [river_max_exp3, river_min_exp3] = max_min_select(river_monthly,year_N)

#   plt.plot(tt, river_max_control-river_min_control, 'bo-')
#   plt.plot(tt, river_max_exp3-river_min_exp3, 'ro-')
#   plt.show()

# ================================================================
# Read simulation results
# ================================================================
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/'
   filename = 'GIAF_control_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control1 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   TLONG = f.variables['TLONG'][:,:].copy()
   TLAT = f.variables['TLAT'][:,:].copy()
   f.close()
   sss_control1[sss_control1>1.e30] = np.nan

   filename = 'GIAF_control_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control2 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   f.close()
   sss_control2[sss_control2>1.e30] = np.nan

   filename = 'GIAF_control_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control3 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   f.close()
   sss_control3[sss_control3>1.e30] = np.nan

   filename = 'GIAF_control_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control4 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   f.close()
   sss_control4[sss_control4>1.e30] = np.nan

   filename = 'GIAF_control_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_control5 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   f.close()
   sss_control5[sss_control5>1.e30] = np.nan

   filename = 'GIAF_exp_1p75_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_1 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp3_1[sss_exp3_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_2 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp3_2[sss_exp3_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_3 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp3_3[sss_exp3_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_4 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp3_4[sss_exp3_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p75_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp3_5 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp3_5[sss_exp3_5>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_1 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp2_1[sss_exp2_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_2 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp2_2[sss_exp2_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_3 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp2_3[sss_exp2_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_4 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp2_4[sss_exp2_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_exp_1p5_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_exp2_5 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_exp2_5[sss_exp2_5>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle1_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_1 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_spin_1[sss_spin_1>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle2_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_2 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_spin_2[sss_spin_2>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle3_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_3 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_spin_3[sss_spin_3>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle4_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_4 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_spin_4[sss_spin_4>1.e30] = np.nan
   f.close()

   filename = 'GIAF_control_spinup_cycle5_select_var.nc'
   f = Dataset(dirname + filename, mode='r')
   sss_spin_5 = f.variables['SALT'][31*12:,0,:,:].data.copy()
   sss_spin_5[sss_spin_5>1.e30] = np.nan
   f.close()

   [sss_control1_trend, sss_control1_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_control1)
   [sss_control2_trend, sss_control2_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_control2)
   [sss_control3_trend, sss_control3_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_control3)
   [sss_control4_trend, sss_control4_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_control4)
   [sss_control5_trend, sss_control5_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_control5)

   [sss_exp3_1_trend, sss_exp3_1_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp3_1)
   [sss_exp3_2_trend, sss_exp3_2_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp3_2)
   [sss_exp3_3_trend, sss_exp3_3_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp3_3)
   [sss_exp3_4_trend, sss_exp3_4_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp3_4)
   [sss_exp3_5_trend, sss_exp3_5_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp3_5)

   [sss_exp2_1_trend, sss_exp2_1_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp2_1)
   [sss_exp2_2_trend, sss_exp2_2_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp2_2)
   [sss_exp2_3_trend, sss_exp2_3_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp2_3)
   [sss_exp2_4_trend, sss_exp2_4_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp2_4)
   [sss_exp2_5_trend, sss_exp2_5_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_exp2_5)

   [sss_spin_1_trend, sss_spin_1_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_spin_1)
   [sss_spin_2_trend, sss_spin_2_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_spin_2)
   [sss_spin_3_trend, sss_spin_3_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_spin_3)
   [sss_spin_4_trend, sss_spin_4_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_spin_4)
   [sss_spin_5_trend, sss_spin_5_sig, pvalue] = sss_calculate_here(TLONG,TLAT,sss_spin_5)

#   sss_control_std = np.max([sss_control1_trend, sss_control2_trend]) - np.min([sss_control1_trend, sss_control2_trend])
   sss_control_std = np.nanstd([sss_control1_trend, sss_control2_trend, sss_control3_trend, sss_control4_trend, sss_control5_trend])
   sss_control_mean = np.nanmean([sss_control1_trend, sss_control2_trend, sss_control3_trend, sss_control4_trend, sss_control5_trend])
   sss_spin_mean = np.nanmean([sss_spin_1_trend, sss_spin_2_trend, sss_spin_3_trend, sss_spin_4_trend, sss_spin_5_trend])

#   sss_exp3_std = np.max([sss_exp3_1_trend, sss_exp3_2_trend]) - np.min([sss_exp3_1_trend, sss_exp3_2_trend])
   sss_exp3_mean = np.nanmean([sss_exp3_1_trend, sss_exp3_2_trend, sss_exp3_3_trend, sss_exp3_4_trend, sss_exp3_5_trend])
   sss_exp2_mean = np.nanmean([sss_exp2_1_trend, sss_exp2_2_trend, sss_exp2_3_trend, sss_exp2_4_trend, sss_exp2_5_trend])

# ================================================================
# Read precipitation and river runoff
# ================================================================
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle1/'
   filename = 'P_R_1970-02.nc'
   f = Dataset(dirname + filename, mode='r')
   lon_in = f.variables['lon'][:].data
   lat_in = f.variables['lat'][:].data
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

   year_N = int(31)

   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle1/'
   [rain_control1_trend, rain_control1_sig, runoff_control1_trend, runoff_control1_sig, rain_range_control1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle2/'
   [rain_control2_trend, rain_control2_sig, runoff_control2_trend, runoff_control2_sig, rain_range_control2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle3/'
   [rain_control3_trend, rain_control3_sig, runoff_control3_trend, runoff_control3_sig, rain_range_control3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle4/'
   [rain_control4_trend, rain_control4_sig, runoff_control4_trend, runoff_control4_sig, rain_range_control4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/control/cycle5/'
   [rain_control5_trend, rain_control5_sig, runoff_control5_trend, runoff_control5_sig, rain_range_control5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle1/'
   [rain_exp2_1_trend, rain_exp2_1_sig, runoff_exp2_1_trend, runoff_exp2_1_sig, rain_range_exp2_1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle2/'
   [rain_exp2_2_trend, rain_exp2_2_sig, runoff_exp2_2_trend, runoff_exp2_2_sig, rain_range_exp2_2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle3/'
   [rain_exp2_3_trend, rain_exp2_3_sig, runoff_exp2_3_trend, runoff_exp2_3_sig, rain_range_exp2_3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle4/'
   [rain_exp2_4_trend, rain_exp2_4_sig, runoff_exp2_4_trend, runoff_exp2_4_sig, rain_range_exp2_4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p5/cycle5/'
   [rain_exp2_5_trend, rain_exp2_5_sig, runoff_exp2_5_trend, runoff_exp2_5_sig, rain_range_exp2_5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle1/'
   [rain_exp3_1_trend, rain_exp3_1_sig, runoff_exp3_1_trend, runoff_exp3_1_sig, rain_range_exp3_1] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,1948)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle2/'
   [rain_exp3_2_trend, rain_exp3_2_sig, runoff_exp3_2_trend, runoff_exp3_2_sig, rain_range_exp3_2] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2010)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle3/'
   [rain_exp3_3_trend, rain_exp3_3_sig, runoff_exp3_3_trend, runoff_exp3_3_sig, rain_range_exp3_3] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2072)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle4/'
   [rain_exp3_4_trend, rain_exp3_4_sig, runoff_exp3_4_trend, runoff_exp3_4_sig, rain_range_exp3_4] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2134)
   dirname = '/stormtrack/data4/yliang/modelling/CESM2/CRU_V7_clm2/extracted/px1p75/cycle5/'
   [rain_exp3_5_trend, rain_exp3_5_sig, runoff_exp3_5_trend, runoff_exp3_5_sig, rain_range_exp3_5] = rain_runoff_calculate_here(year_N,dirname,mask_rg,area_in,2196)

# ================================================================
# Plot figures
# ================================================================
if True:

   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   tt2 = np.linspace(1979,2009,31)

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

   plt.plot([1, 1], [sss_core2_trends[-2]*10, sss_core2_trends[-2]*10], 'm*', label='CORE-II NCAR')
   plt.plot([1, 1], [sss_spin_mean*10, sss_spin_mean*10], 'k*', label='spin-up')
   plt.axis([0.75, 2, 0.025, 0.05])
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

