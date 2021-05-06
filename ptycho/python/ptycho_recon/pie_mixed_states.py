import scipy.io as sio #for read/write matlab file
from scipy import ndimage
import utility_function as u
import utility_function_recon as ur
import pyfftw
from numpy import *
from numpy.fft import *
import numpy as np
import time

########## mixed states ##########
def reconPIE_mixed_state(dp, paraDict):
    N_probe = paraDict['N_probe']; N_object = paraDict['N_object']
    N_roi = paraDict['N_roi']
    Niter = paraDict['Niter'];
    N_tot = N_roi*N_roi
    if 'randomSeed' in paraDict: random.seed(mod(int(paraDict['randomSeed']),4294967295))

    auxiFunc = ur.auxiliary_function(paraDict)

    ########## initialize exit wave ########## 
    psi = zeros((N_probe, N_roi, N_roi), dtype=np.complex128)
    psi_old = zeros((N_probe, N_roi, N_roi), dtype=np.complex128)
    delta_psi = zeros((N_probe, N_roi, N_roi), dtype=np.complex128)
    
    ########## initialize CBED ##########
    dp_tot = sum(dp)
    dp_avg = np.sum(dp * dp, axis = 0)/paraDict['N_scan']
    cbed_region_mag = zeros((N_probe, paraDict['N_dp'], paraDict['N_dp']))

    ########## initialize transmission function ##########
    O = auxiFunc.initializeObject()

    ########## initialize probe function ##########
    if 'previous_probe' in paraDict:
        probes = paraDict['previous_probe'].copy()
    elif 'probes0' in paraDict:
        probes = paraDict['probes0'].copy()
    else:
        probes = np.zeros((N_probe, N_roi, N_roi), dtype=np.complex128)
        if paraDict['normalizeInitialProbe']:
            #print('normalize initial probe intensity to match CBED')
            probe = paraDict['probe0'] * sqrt(np.sum(dp_avg) / np.sum(abs(paraDict['probe0'])**2) / N_tot)
        else:
            probe = paraDict['probe0']
        for i in range(N_probe):
            probes[i,:,:] = probe / (i+1)

    probes_shifted = np.zeros((N_probe,N_roi,N_roi), dtype=np.complex128)
    probes_old = np.zeros((N_probe,N_roi,N_roi), dtype=np.complex128)

    ########## miscellaneous ##########
    startNiter = 0
    timeLeft = 0
    start_time = time.time()
    time_counter = 1

    #Use single object and probe before mixed states update
    N_object_recon = 1
    N_probe_recon = 1

    dp_error = zeros(paraDict['N_scan']) #difference between data and recon wave
    s, dp_error_old = auxiFunc.initializeDataError()

    if 'previousIteration' in paraDict: 
        startNiter = paraDict['previousIteration']
        if startNiter >= paraDict['Niter_update_states']:
            N_object_recon = N_object; N_probe_recon = N_probe

    ################################# prepare result dictionary ################################# 
    resultDir = {'object':O}
    resultDir['psi'] = psi
    resultDir['probes'] = probes
    resultDir['probe0'] = paraDict['probe0'].copy()
    if 'probes0' in paraDict: resultDir['probes0'] = paraDict['probes0'].copy()
    resultDir['dx_x'] = 1.0/(paraDict['dk_x']*N_roi); 
    resultDir['dx_y'] = 1.0/(paraDict['dk_y']*N_roi);
    resultDir['dk_x'] = paraDict['dk_x']; resultDir['dk_y'] = paraDict['dk_y']
    resultDir['ppX'] = paraDict['ppX']; resultDir['ppY'] = paraDict['ppY']
    resultDir['dp_avg'] = dp_avg
    resultDir['s'] = s
    resultDir['badPixels'] = paraDict['badPixels']
    
    resultDir['dp_error'] = dp_error
    if 'filter_r_probe' in paraDict:
        resultDir['filter_r_probe'] = paraDict['filter_r_probe']
    if 'filter_r_psi' in paraDict:
        resultDir['filter_r_psi'] = paraDict['filter_r_psi']
    if 'filter_f_probe' in paraDict:
        resultDir['filter_f_probe'] = paraDict['filter_f_probe']
    if 'filter_f_psi' in paraDict:
        resultDir['filter_f_psi'] = paraDict['filter_f_psi']
    if 'probe0_info' in paraDict:
        resultDir['probe0_info'] = paraDict['probe0_info']

    ################################## main recon loop #################################    
    for k in range(startNiter, Niter):
        if mod(k, paraDict['Niter_print'])==0: auxiFunc.printStatus(timeLeft, k)

        if k == paraDict['Niter_update_probe']: 
            print('start probe update')
            start_time = time.time()
            time_counter = 1
        if k == paraDict['Niter_update_position']: 
            print('start position correction')
            start_time = time.time()
            time_counter = 1
        if k == paraDict['Niter_update_states']:
            print('start mixed states update')
            N_object_recon = N_object; N_probe_recon = N_probe
            start_time = time.time()
            time_counter = 1
            
            probes_temp = gramschmidt(probes.reshape(N_probe, N_roi*N_roi))
            probes[:,:,:] = probes_temp.reshape(N_probe, N_roi, N_roi)
           
        update_order = random.permutation(paraDict['N_scan']) #random order
 
        for i in update_order:
            O_old = auxiFunc.getObjectROI(O, i)
            for p in range(N_probe_recon):
                probes_shifted[p,:,:] = auxiFunc.shiftProb(probes[p,:,:], i, 'toScanPosition')

                #overlap projection
                psi[p,:,:] = O_old * probes_shifted[p,:,:]

                #Fourier projection
                psi_old[p,:,:] = psi[p,:,:]
                psi[p,:,:], cbed_region_mag[p,:,:] = auxiFunc.FFTpsi(psi[p,:,:])

            psi_f_mag_tot = sqrt(np.sum(cbed_region_mag**2, axis = 0))
            dp_error[i] = np.sum((psi_f_mag_tot - dp[i,:,:])**2)

            #Fourier projection
            for p in range(N_probe_recon):
                psi[p,:,:] = auxiFunc.updateFourierIntensity(psi[p,:,:], dp[i,:,:], psi_f_mag_tot)
                if 'filter_r_psi' in paraDict: psi[p,:,:] = psi[p,:,:] * paraDict['filter_r_psi']

            delta_psi = psi - psi_old
            probe_old = probes_shifted[:,:,:]
            O_update = auxiFunc.calculateMixedStatesUpdate(probe_old, delta_psi, 'o')
            auxiFunc.updateObj(O, O_update, i)

            if k>=paraDict['Niter_update_probe']:
                O_tot_max = amax(abs(O_old)**2)
                for p in range(N_probe_recon):
                    probe_update = auxiFunc.calculateUpdate(O_old, delta_psi[p,:,:], 'p')
                    probes_shifted[p,:,:] += probe_update
                    probes[p,:,:] = auxiFunc.shiftProb(probes_shifted[p,:,:], i, 'toOrigin', checkFilter=True)

            if k>=paraDict['Niter_update_position']:
                auxiFunc.gradPositionCorrection(probe_old[0,:,:], O_old, i, delta_psi[0,:,:])

        ############################### orthogonalise probe #######################################
        if k >= paraDict['Niter_update_states']:
            probes_temp = gramschmidt(probes.reshape(N_probe, N_roi*N_roi))
            probes[:,:,:] = probes_temp.reshape(N_probe, N_roi, N_roi)
           
        ############################### calcuate data error #######################################
        s[k] = np.sum(dp_error)/dp_tot
        dp_error_old = dp_error.copy()

        ############################### save results #######################################
        if mod(k+1, paraDict['Niter_save'])==0: #save results
            saveName = paraDict['saveName'] + '_Niter'+str(k+1) + '.mat'
            sio.savemat(saveName,resultDir)
    
        timeLeft = (time.time()-start_time)/time_counter * (Niter-k-1)
        time_counter += 1

def proj(u, v):
    return u * np.vdot(u,v) / np.vdot(u,u)  

def gramschmidt(V):
    U = np.copy(V)
    for i in range(1, V.shape[0]):
        for j in range(i):
            U[i,:] -= proj(U[j,:], V[i,:])
    return U

