import numpy as np
import os #for change directory
import scipy.io as sio #for read/write matlab file
import scipy.ndimage.filters as sfilter
import scipy.ndimage
import scipy.misc

import ptycho_recon.ptycho as pty
import ptycho_recon.utility_function as utils
import ptycho_recon.probe as probe

#import h5py

from numpy import *
import time

##################################### read data ########################################
print("load data")
#data_name = 'data_ws2_wse2_80keV_09_roi2.mat'
data_name = 'data_mos2_sample.mat'
currentdir = os.getcwd() #current directory
#os.chdir(data_dir) #change to data directory
data = sio.loadmat(data_name)
dp = data['dp']*1.0

sx = 0 if not 'sx' in data else int(squeeze(data['sx']))
sy = 0 if not 'sy' in data else int(squeeze(data['sy']))
#if not 'sx' in data else int(squeeze(data['sx']))
#sy = 0 if not data.has_key('sy') else int(squeeze(data['sy']))
######################################### Parameters ##################################### 
ADU_background_cutoff = 20
ADU_electronCount_ratio = 149.0

N_roi = 128

##############################################################################################
N_dp = dp.shape[0]
voltage = squeeze(data['voltage']) #kev
alpha_max = squeeze(data['alpha_max']) #mrad
df = squeeze(data['df'])
cs = squeeze(data['cs']) #angstrom
#scanStepSize_x = squeeze(data['scanStepSize_x']) #angstrom
#scanStepSize_y = squeeze(data['scanStepSize_y'])
scanStepSize_x = 0.21 #angstrom
scanStepSize_y = 0.21

dk = squeeze(data['dk'])

print("dk =", dk)
dx = 1.0/dk/N_roi
print("dx =", dx)

################################## data processing ##########################################
rot_angle_d = 30

print("processing data")
result_dir_extra = "/preprocessCBED"
dp, result_dir_extra = utils.transpose_cbed(dp, result_dir_extra)
dp, result_dir_extra = utils.background_removal(dp, ADU_background_cutoff, result_dir_extra)

print("recon data size:", dp.shape)
################################## make initial probe function ##################################
probe_init = probe.STEMprobe()
probe_init.df = df
probe_init.cs = cs
probe_init.alpha_max = alpha_max
probe_init.voltage = voltage

################################## calculate scan positions ####################################
#calculate scan positions
N_scan_y = dp.shape[2]
N_scan_x = dp.shape[3]

ppX, ppY, result_dir_extra = utils.calculate_scan_positions(N_scan_x, N_scan_y, scanStepSize_x, scanStepSize_y, rot_angle_d, result_dir_extra)

Ny_max = max([abs(round(np.min(ppY)/dx)-floor(N_roi/2.0)), abs(round(np.max(ppY)/dx)+ceil(N_roi/2.0))])*2+1
Nx_max = max([abs(round(np.min(ppX)/dx)-floor(N_roi/2.0)), abs(round(np.max(ppX)/dx)+ceil(N_roi/2.0))])*2+1
N_image = int(max([Ny_max,Nx_max]))+20
print("Image size:", N_image)
'''
######### reshape dp and scan positions #########
N_scan_y = dp.shape[2]
N_scan_x = dp.shape[3]
N_scan_tot = N_scan_y * N_scan_x

dp_temp = zeros((N_scan_tot,dp.shape[0],dp.shape[1]))
ppX = ppX.reshape(N_scan_tot)
ppY = ppY.reshape(N_scan_tot)

for i in range(N_scan_y):
    for j in range(N_scan_x):
        index = i*N_scan_x + j
        dp_temp[index,:,:] = sqrt(dp[:,:,i,j])
'''
########################################### reconstruction #####################################
print('Reconstruction')
reconObject = pty.ptycho(dp, dk, probe_init, ppX, ppY)
reconObject.paraDict['N_image'] = N_image
reconObject.paraDict['N_roi'] = N_roi
reconObject.paraDict['Niter'] = 50
reconObject.paraDict['Niter_update_probe'] = 0
reconObject.paraDict['Niter_save'] = 5

reconObject.paraDict['rotationAngle']  = rot_angle_d

reconObject.paraDict['printID'] = 'MoS2'

############## position correction ##############
reconObject.paraDict['Niter_update_position'] = 30

############## mixed-states recon ##############
reconObject.paraDict['N_probe'] = 2
reconObject.paraDict['Niter_update_states'] = 10

result_dir = currentdir + "/mos2"
result_dir_extra = reconObject.initialize(result_dir)

start_time = time.time()
reconObject.recon() #start reconstruction
total_time = time.time() - start_time

timeLeftMin, timeLeftSec = divmod(total_time, 60)
timeLeftHour, timeLeftMin = divmod(timeLeftMin, 60)
print('Total recon time: %02d:%02d:%02d' %(timeLeftHour,timeLeftMin,timeLeftSec))


