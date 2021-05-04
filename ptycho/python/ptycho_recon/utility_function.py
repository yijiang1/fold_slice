import numpy as np
from numpy import *
from scipy import ndimage
import scipy.ndimage
import pyfftw
from numpy.fft import *
import pyfftw

import zipfile as zp
import os
import warnings

def readraw(filename):
    scanx = int(filename.rstrip('.raw').split('_')[-1].lstrip('abcdefghijklmnopqrstuvwxyz'))
    scany = int(filename.rstrip('.raw').split('_')[-2].lstrip('abcdefghijklmnopqrstuvwxyz'))
    contents = np.fromfile(filename, dtype = 'float32')
    data_arr = np.reshape(contents, (125, 125, scany, scanx), order = 'C')
    return data_arr
    
##############################################################################################
def shift(input,dx,px,py):
    N_image = input.shape[0]

    dk= 1.0/(dx*N_image)
    kx = np.linspace(-np.floor(N_image/2.0),np.ceil(N_image/2.0)-1,N_image)
    [kX,kY] = np.meshgrid(kx,kx)
    kX = kX*dk; kY = kY*dk;

    f = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(input)))
    f = f*np.exp(-2*np.pi*1j*px*kX)*np.exp(-2*np.pi*1j*py*kY)
    f = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(f)))

    return f

##############################################################################################
def convertADUtoElectronCount(input, ADU_electronCount_ratio, directory = ""):
    print("converting ADU to # of electrons:", ADU_electronCount_ratio)
    output = input/ADU_electronCount_ratio
    directory = directory + "_ADUtoElectron" + str(ADU_electronCount_ratio)
    return output, directory

##############################################################################################
def partition_scan(N_scan_x, N_scan_y, partion_x, partion_y):
    scan_partition = np.zeros((N_scan_y, N_scan_x), dtype=np.int)
    delta_x = int(np.ceil(N_scan_x / partion_x))
    delta_y = int(np.ceil(N_scan_y / partion_y))
    value = 0
    for i in xrange(partion_y):
        for j in xrange(partion_x):
             index_x_lb = delta_x*j
             index_x_ub = min(delta_x*(j+1), N_scan_x)
             index_y_lb = delta_y*i
             index_y_ub = min(delta_y*(i+1), N_scan_y)
             scan_partition[index_y_lb:index_y_ub,index_x_lb:index_x_ub] = value
             value = value + 1
    return scan_partition

##############################################################################################
def shift_cbed(input, px, py, directory = ""):
    N_cbed = input.shape[0]
    N_tot = N_cbed*N_cbed
    dx = 1.0
    dk= 1.0/(dx*N_cbed)
    kx = np.linspace(-np.floor(N_cbed/2.0),np.ceil(N_cbed/2.0)-1,N_cbed)
    kx = ifftshift(kx)
    [kX,kY] = np.meshgrid(kx,kx)
    kX = kX*dk; kY = kY*dk;
    output = zeros(input.shape)
    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            f = np.fft.fft2(input[:,:,i,j]);
            f = f*exp(-2*pi*1j*px*kX)*exp(-2*pi*1j*py*kY)
            output[:,:,i,j] = abs(np.fft.ifft2(f)); #fix normalization
    directory = directory + "_shift_sx" + str(px)+"_sy" + str(py)
    return output, directory

##############################################################################################
def lowpassfilter(input, dx,cutoff):
    N = input.shape[0]
    dk= 1.0/(dx*N)
    kx = np.linspace(-np.floor(N/2.0),np.ceil(N/2.0)-1,N)
    [kX,kY] = np.meshgrid(kx,kx)
    kR = np.sqrt(kX**2 + kY**2)

    a = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(input)))
    a[kR>(N/2*cutoff)] = 0
    output = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(a)))

    return output

##############################################################################################
def lowpassfilter_alpha(input, cutoff, dk_x, dk_y, alpha_max, voltage):
    N = input.shape[0]
    print("applying low pass filter to image")
    print("mask cutoff =",cutoff ,'alpha')
    kx = linspace(-floor(N/2.0),ceil(N/2.0)-1, N)
    [kX,kY] = meshgrid(kx,kx)
    kX = kX*dk_x; kY = kY*dk_y;
    kR = np.sqrt(kX**2+ kY**2)
    
    wavelength = 12.398/np.sqrt((2*511.0 + voltage) * voltage)  #angstrom

    k_cutoff = cutoff * alpha_max *1e-3 / wavelength
    
    f = fftshift(fft2(ifftshift(input)))
    f[kR > k_cutoff] = 0
    output = real(fftshift(ifft2(ifftshift(f))))

    return output
##############################################################################################
def propagtor_function(N, dk_x, dk_y,dz, wavelength):
    kx = np.linspace(-np.floor(N/2.0),np.ceil(N/2.0)-1,N)
    [kX,kY] = np.meshgrid(kx,kx)
    kX = kX*dk_x; kY = kY*dk_y;

    kR = np.sqrt(kX**2 + kY**2)
    cutoff = (np.ceil(N/2.0)-1)*min(dk_x,dk_y)*2/3

    P = zeros((dz.size,N,N), dtype=np.complex128)
    for i in range(dz.size):
        temp = np.exp(-1j*np.pi*wavelength*kR**2*dz[i])
        temp[kR>cutoff] = 0  #apply a low pass filter to keep 2/3 of maximum spatial frequency
        P[i,:,:] = ifftshift(temp)

    '''
    x = np.linspace(-np.floor(N/2.0),np.ceil(N/2.0)-1,N)*dx
    [X,Y] = np.meshgrid(x,x)
    R = np.sqrt(X**2 + Y**2)
    p = 1/(1j*wavelength*dz)*np.exp(1j*np.pi/(wavelength*dz)*R**2)
    '''
    return P

##############################################################################################
def recenter(input):
    N = input.shape[0]
    center_index = N // 2
    [yy,xx] = where(abs(input) == np.max(abs(input)))

    output = roll(input,int(-(yy[0]-center_index)), axis = 0)
    output = roll(output,int(-(xx[0]-center_index)), axis = 1)

    return output

##############################################################################################
def upsample_cbed_ff(dp, dk_x, dk_y, resizeFactor = 2, directory = ""):
    print("upsample cbed using free float ptychography")
    print("old dk_x=", dk_x, "old dk_y=", dk_y)

    #resize data
    output = np.zeros((int(dp.shape[0]*resizeFactor), int(dp.shape[1]*resizeFactor), dp.shape[2], dp.shape[3]))
    output[0:-1:resizeFactor,0:-1:resizeFactor,:,:] = dp
    mask = np.ones((int(dp.shape[0]*resizeFactor), int(dp.shape[1]*resizeFactor)))
    print(mask.shape)
    mask[0:-1:2,0:-1:2] = 0
    dk_x_r = dk_x/resizeFactor
    dk_y_r = dk_y/resizeFactor
    print("new dk_x=", dk_x_r, "new dk_y=", dk_y_r)
    directory = directory + "_upsampleCBED" + str(resizeFactor)
    #output[output<0] = 0
    return output, mask, dk_x_r, dk_y_r, directory


##############################################################################################
def resize_cbed(dp, resizeFactor, dk_x, dk_y, directory = "", order = 1):
    print("resize cbed")
    print("old dk_x=", dk_x, "old dk_y=", dk_y)

    #resize data
    output = np.zeros((int(dp.shape[0]*resizeFactor),int(dp.shape[1]*resizeFactor),dp.shape[2],dp.shape[3]))
    for i in range(0,dp.shape[2]):
        for j in range(0,dp.shape[3]):
            scipy.ndimage.interpolation.zoom(dp[:,:,i,j],[resizeFactor,resizeFactor],output[:,:,i,j], order)
    dk_x_r = dk_x/resizeFactor
    dk_y_r = dk_y/resizeFactor
    print("new dk_x=", dk_x_r, "new dk_y=", dk_y_r)
    directory = directory + "_resizeCBED" + str(resizeFactor)
    #output[output<0] = 0
    return output, dk_x_r, dk_y_r, directory

##############################################################################################
def resample_cbed(dp, Npix, dk_x, dk_y, directory = ""):
    print("resample cbed using every ", str(Npix), 'pixels...')
    print("old dk_x=", dk_x, "old dk_y=", dk_y)
    output = dp[0:-1:Npix,0:-1:Npix,:,:]
    if Npix>1: directory = directory + "_resampleCBED" + str(Npix)
    dk_x_new = dk_x*Npix; dk_y_new = dk_y*Npix
    print("new dk_x=", dk_x_new, "new dk_y=", dk_y_new)
    return output, dk_x_new, dk_y_new, directory
    
##############################################################################################
def crop_cbed(dp, N_dp_x_new, N_dp_y_new, directory = ""):
    print("crop cbed to ", str(N_dp_y_new), 'x', str(N_dp_x_new))
    cen_x = floor(dp.shape[1]/2.0)
    cen_y = floor(dp.shape[0]/2.0)
    index_x_lb = (cen_x - floor(N_dp_x_new/2.0)).astype(np.int)
    index_x_ub = (cen_x + ceil(N_dp_x_new/2.0)).astype(np.int)
    index_y_lb = (cen_y - floor(N_dp_y_new/2.0)).astype(np.int)
    index_y_ub = (cen_y + ceil(N_dp_y_new/2.0)).astype(np.int)
    
    #crop data
    output = np.zeros((N_dp_x_new, N_dp_y_new,dp.shape[2],dp.shape[3]))
    for i in range(0,dp.shape[2]):
        for j in range(0,dp.shape[3]):
            output[:,:,i,j] = dp[index_y_lb:index_y_ub,index_x_lb:index_x_ub,i,j]
    directory = directory + "_crop_Ndpx" + str(N_dp_x_new) + "_Ndpy" + str(N_dp_y_new)
    return output, directory

##############################################################################################
def pad_cbed(input, N_pad_x, N_pad_y, value = 0, directory = ""):
    print("pad " + str(value) + "s to cbed patterns...")
    Ny = input.shape[0]; Nx = input.shape[1];
    pad_pre_y = int(np.ceil((N_pad_y - Ny) / 2.0))
    pad_post_y = int(np.floor((N_pad_y - Ny) / 2.0))
    pad_pre_x = int(np.ceil((N_pad_x - Nx) / 2.0))
    pad_post_x = int(np.floor((N_pad_x - Nx) / 2.0))
  
    output = np.pad(input, ((pad_pre_y, pad_post_y), (pad_pre_x, pad_post_x), (0,0), (0,0)), 'constant', constant_values=value)
    directory = directory + "_padCBED" + str(value) + "_" + str(N_pad_x)
    return output, directory

##############################################################################################
def apply_circular_mask(input, radius, offset_x=0, offset_y=0, directory = ""):
    N_dp = input.shape[0]
    print("applying circullar mask to cbed")
    print("mask radius =",radius)
    print("mask offset_x =",offset_x, "mask offset_y =",offset_y)

    x = np.linspace(-np.floor(N_dp/2.0),np.ceil(N_dp/2.0)-1,N_dp)
    [X,Y] = np.meshgrid(x,x)
    mask_disk = np.sqrt(X**2+ Y**2)
    
    output = input.copy()

    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            temp = input[:,:,i,j].copy()
            temp = np.roll(temp, offset_y, axis=0)
            temp = np.roll(temp, offset_x, axis=1)
            if radius>0:
               temp[mask_disk>radius] = 0
            output[:,:,i,j] = temp.copy()
    if radius>0: directory = directory + "_lowPassFilter" + str(radius)
    if offset_x!= 0: directory = directory + "_sx" + str(offset_x)
    if offset_y!= 0: directory = directory + "_sy" + str(offset_y)
    
    return output, directory, mask_disk

##############################################################################################
def apply_circular_mask_alpha(dp, cutoff, dk_x, dk_y, alpha_max, voltage, offset_x=0, offset_y=0, directory = ""):
    if cutoff>0:
        N_dp = dp.shape[0]
        print("applying circullar mask to cbed")
        print("mask cutoff =",cutoff ,'alpha')
        print("mask offset_x =",offset_x, "mask offset_y =",offset_y)
        kx = linspace(-floor(N_dp/2.0),ceil(N_dp/2.0)-1, N_dp)

        [kX,kY] = meshgrid(kx,kx)
        kX = kX*dk_x; kY = kY*dk_y;
        kR = np.sqrt(kX**2+ kY**2)
        
        wavelength = 12.398/np.sqrt((2*511.0 + voltage) * voltage)  #angstrom

        k_cutoff = cutoff * alpha_max *1e-3 / wavelength
    output = dp.copy()
    for i in range(dp.shape[2]):
        for j in range(dp.shape[3]):
            temp = dp[:,:,i,j].copy()
            temp = np.roll(temp, offset_y, axis=0)
            temp = np.roll(temp, offset_x, axis=1)
            #np.roll(dp[:,:,i,j], offset_y, axis=0)
            #np.roll(dp[:,:,i,j], offset_x, axis=1)
            if cutoff > 0:
               #dp[kR > k_cutoff,i,j]= 0
               temp[kR > k_cutoff] = 0
            output[:,:,i,j] = temp.copy()
    directory = directory + "_cutoff" + str(cutoff)+"alpha"
    if offset_x!= 0: directory = directory + "_sx" + str(offset_x)
    if offset_y!= 0: directory = directory + "_sy" + str(offset_y)
    
    return output, directory
##############################################################################################
def add_poisson_noise(input, current, readOutTime, directory = ""):
    N_dp_tot = input.shape[0] * input.shape[1]
    print("applying poisson noise to cbed")
    print("beam current =", current, 'pA')
    Nc_avg = current*1e-12*readOutTime/(1.6e-19)/N_dp_tot
    print("average count per pixel = ", Nc_avg) 
    print("snr = ", sqrt(Nc_avg)) 
    
    output = input.copy()
    snr = zeros((input.shape[2],input.shape[3]))
    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            cbed_noise = input[:,:,i,j].copy()
            cbed_noise = cbed_noise/np.sum(input[:,:,i,j])*(N_dp_tot*Nc_avg*1.0)
            cbed_noise = random.poisson(cbed_noise)
            cbed_noise = cbed_noise*np.sum(input[:,:,i,j])/(N_dp_tot*Nc_avg*1.0)
            output[:,:,i,j] = cbed_noise
            snr[i,j] = np.mean(input[:,:,i,j])/np.std(input[:,:,i,j] - cbed_noise)
    directory = directory + "_poissonNoise" + str(current) + "pA"
    return output, snr, directory

##############################################################################################
def average_cbed(input, windowSize, scanStepSize_x, scanStepSize_y, directory = ""):
    print("average cbed patterns... window size =",windowSize)
    
    if windowSize==1:
        output = input
    else:
        output = zeros((input.shape[0],input.shape[1],input.shape[2]//windowSize,input.shape[3]//windowSize))
        for i in range(output.shape[2]):
            for j in range(output.shape[3]):
                temp = np.sum(input[:,:,i*windowSize:(i+1)*windowSize,j*windowSize:(j+1)*windowSize], axis=(2,3))/windowSize**2
                output[:,:,i,j] = temp.copy()
        directory = directory + "_averageCBED"+str(windowSize)
        scanStepSize_x = scanStepSize_x * windowSize
        scanStepSize_y = scanStepSize_y * windowSize
    
    return output, scanStepSize_x, scanStepSize_y, directory

##############################################################################################
def resample_scan(input, windowSize, scanStepSize_x, scanStepSize_y, directory = ""):
    print("resample scans ... window size =",windowSize)
     
    if windowSize==1:
        output = input
    else:
        output = input[:,:,::windowSize,::windowSize]
        
        scanStepSize_x = scanStepSize_x * windowSize
        scanStepSize_y = scanStepSize_y * windowSize

        directory = directory + "_resampleScan"+str(windowSize)

    return output, scanStepSize_x, scanStepSize_y, directory

##############################################################################################
def transpose_cbed(input, directory = ""):
    print("transpose cbed patterns...")
    output = zeros((input.shape[1],input.shape[0],input.shape[2],input.shape[3]))

    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            output[:,:,i,j] = input[:,:,i,j].T
    directory = directory + "_transpose" 
    return output, directory

##############################################################################################
def rot90_cbed(input, k, directory = ""):
    print("rotate cbed patterns by", 90*k, "degrees...")
    output = zeros((input.shape[1],input.shape[0],input.shape[2],input.shape[3]))
    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            output[:,:,i,j] = rot90(input[:,:,i,j], k)
    directory = directory + "_rotate90_"+str(k) 
    return output, directory

##############################################################################################
def transpose_scan_positions(input, directory = ""):
    print("transpose scan positions...")
    output = zeros((input.shape[0],input.shape[1],input.shape[3],input.shape[2]))

    for i in range(output.shape[2]):
        for j in range(output.shape[3]):
            output[:,:,i,j] = input[:,:,j,i]

    directory = directory + "_transposeScanPos"
    return output, directory

##############################################################################################
def flip_scan_positions(input, flipType, directory = ""):
    assert flipType in ['lr','ud'], "Flip type %s not known!" % flipType
    if flipType=="ud":
        print("flip scan positions up and down (y-axis, third dimension)")
        output = input[:,:,::-1,:]
    if flipType=="lr":
        print("flip scan positions left and right (x-axis, forth dimension)")
        output = input[:,:,:,::-1]

    directory = directory + "_flipScanPos_" + flipType
    return output, directory
##############################################################################################
def normalize_wave_function(input, dx):
    output = input.copy()
    c = sqrt( 1.0/ ( np.sum(np.abs(input)**2 * dx**2 ) ))
    output = output * c
    return output

##############################################################################################
def normalize_cbed(input, dk, directory = ""):
    print("normalizing cbed patterns: sum(dp) = 1")
    
    for i in range(input.shape[2]):
        for j in range(input.shape[3]):
            c = 1.0/(sum(abs(input[:,:,i,j])))
            input[:,:,i,j] = input[:,:,i,j] * c
    directory = directory + "_normCBED"
    return input, directory
    
##############################################################################################
def background_removal(input, bg_level, directory = ""):
    print("removing background... level=", bg_level)
    output = input.copy()
    output[output<=bg_level] = 0
    directory = directory + "_bgRemoval"+str(bg_level)
    return output, directory

##############################################################################################
def background_subtraction(input, bg_level, directory = ""):
    print("subtracting background... level=", bg_level)

    output = input - bg_level
    output[output<0] = 0

    directory = directory + "_bgSubtraction"+str(bg_level)
    return output, directory

##############################################################################################
def calculate_scan_positions(N_scan_x, N_scan_y, scanStepSize_x, scanStepSize_y, rot_angle_d = 0, directory = "", randomOffset = 0, ppX = 0, ppY = 0):
    print("calculating scan positions...")
    print("N_scan_x =", N_scan_x, "scanStepSize_x =", scanStepSize_x)
    print("N_scan_y =", N_scan_y, "scanStepSize_y =", scanStepSize_y)
    print("rot_angle =", rot_angle_d)
    rot_angle = rot_angle_d*pi/180.0

    ppx = linspace(-floor(N_scan_x/2.0),ceil(N_scan_x/2.0)-1,N_scan_x)*scanStepSize_x
    ppy = linspace(-floor(N_scan_y/2.0),ceil(N_scan_y/2.0)-1,N_scan_y)*scanStepSize_y
    [ppX0, ppY0] = meshgrid(ppx,ppy)
    if not isscalar(ppX): ppX0 = ppX
    if not isscalar(ppY): ppY0 = ppY

    if randomOffset > 0:
        ppX0 = ppX0 + (np.random.rand(ppX0.shape[0], ppX0.shape[1])*2-1)*scanStepSize_x*randomOffset
        ppY0 = ppY0 + (np.random.rand(ppY0.shape[0], ppY0.shape[1])*2-1)*scanStepSize_y*randomOffset

    ppY_rot = ppX0*-sin(rot_angle) + ppY0*cos(rot_angle)
    ppX_rot = ppX0*cos(rot_angle) + ppY0*sin(rot_angle)

    directory = directory + "/scanStepSize" + str(np.around(scanStepSize_x,4)) + "_rotAngle"+str(rot_angle_d)
    if randomOffset>0: directory = directory + "_randomOffset" + str(randomOffset)
    if not isscalar(ppX) or not isscalar(ppY): directory = directory + "_externalCoord"

    return ppX_rot, ppY_rot, directory

##############################################################################################
def gaussian(N, sigma):
    x = linspace(-floor(N/2.0),ceil(N/2.0)-1, N)

    [X,Y] = meshgrid(x,x)
    g = exp(-(X**2+Y**2)/(2*sigma**2)) + np.zeros((N,N), dtype=np.complex128)
    return g

##############################################################################################
def guess_bad_scan(I, threshold):
    print('Determining bad scans')

    I_pad = np.lib.pad(I, (1, 1), 'edge')
    # calculate standard deviation in a 3 x 3 window
    averageI2 = scipy.ndimage.filters.uniform_filter(I_pad ** 2)
    averageI = scipy.ndimage.filters.uniform_filter(I_pad)
    std = np.sqrt(abs(averageI2 - averageI**2))[1:-1, 1:-1]

    medianI = scipy.ndimage.filters.median_filter(I_pad, 2)[1:-1, 1:-1]

    return abs(I - medianI) > std * threshold

