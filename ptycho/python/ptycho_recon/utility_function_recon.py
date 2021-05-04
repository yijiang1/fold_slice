import numpy as np
from numpy import *
from scipy import ndimage
import scipy.ndimage
import pyfftw
from numpy.fft import *
import multiprocessing

class auxiliary_function:
      def __init__(self, paraDict):
          self.paraDict = paraDict          
          self.createFourierCoord(self.paraDict['N_roi'])
          self.initializeFFTW(self.paraDict['N_roi'])
          self.initializeScanPositions()
          self.initializeDiffractionPatterns()
          self.alpha = self.paraDict['alpha']
          self.beta = self.paraDict['beta']

      def calculateUpdate(self, psi, delta_psi, type):
          psi_mag = abs(psi)**2
          if type=='o': #calculate update for object
              w = self.alpha
          elif type=='p': #calculate update for probe/psi
              w = self.beta
          else: raise RuntimeError('Invalid input!')

          return w * conj(psi)/ amax(psi_mag) * delta_psi

      def calculateMixedStatesUpdate(self, psi, delta_psi, type):
          psi_tot = sum(abs(psi)**2,axis=0)
          #psi_mag = abs(psi)**2
          if type=='o': #calculate update for object
              w = self.alpha
          elif type=='p': #calculate update for probe/psi
              w = self.beta
          else: raise RuntimeError('Invalid input!')

          return w * sum(conj(psi)*delta_psi, axis = 0)/amax(psi_tot)

      #################### object ####################
      def initializeObject(self):
          if 'previous_obj' in self.paraDict:
              O = self.paraDict['previous_obj']
          else:
              N_image = self.paraDict['N_image']
              if self.paraDict['uniformInitialObject']:
                  O = np.ones((N_image, N_image), dtype=np.complex128)
              else:
                  O = np.random.rand(N_image,N_image) + 1j*np.random.rand(N_image,N_image)
              O = O/abs(O)
          return O

      def updateObj(self, O, objUpdate, ind_dp):
          O[self.ind_y_lb_s[ind_dp]:self.ind_y_ub_s[ind_dp],self.ind_x_lb_s[ind_dp]:self.ind_x_ub_s[ind_dp]] += objUpdate

      def getObjectROI(self, O, ind_dp):
          return O[self.ind_y_lb_s[ind_dp]:self.ind_y_ub_s[ind_dp],self.ind_x_lb_s[ind_dp]:self.ind_x_ub_s[ind_dp]]

      #################### probe ####################
      def shiftProb(self, probe, ind_dp, direction, checkFilter=False):
          if direction=='toScanPosition': #from origin to scan position
              px = self.px_f[ind_dp]
              py = self.py_f[ind_dp]
          elif direction=='toOrigin': #from scan position back to origin
              px = -self.px_f[ind_dp]
              py = -self.py_f[ind_dp]
          else: raise RuntimeError('Invalid input!')

          self.r[:,:] = probe
          self.fft_forward.update_arrays(self.r, self.f)
          self.fft_forward.execute()
          self.f = self.f*exp(-2*pi*1j*px*self.kX)*exp(-2*pi*1j*py*self.kY)
          if checkFilter and 'filter_f_probe' in self.paraDict: 
              self.f = self.f * self.paraDict['filter_f_probe']
          self.fft_inverse.update_arrays(self.f, self.r)
          self.fft_inverse.execute();
          return self.r / self.N_tot #fix normalization

      def orthoProbe(self, probes):
          probes_temp = gramschmidt(probes.reshape(paraDict['N_probe'], paraDict['N_roi']**2))
          probes[:,:,:] = probes_temp.reshape(paraDict['N_probe'], paraDict['N_roi'], paraDict['N_roi'])
          #sort probes based on power
          power = sum(abs(probes)**2, axis=(1,2))
          power_ind = argsort(-power)
          probes[:,:,:] = probes[power_ind,:,:]
          return probes

      def processProbe(self, probe):
          if 'filter_r_probe' in self.paraDict: probe = probe * paraDict['filter_r_probe']
          if 'probe_profile' in self.paraDict:
              probe_mag_sum = sum(abs(probe))
              probe = probe / abs(probe) * self.paraDict['probe_profile']
              probe = probe / sum(abs(probe)) * probe_mag_sum
          return probe

      def FFTpsi(self, psi):
          self.r[:,:] = psi
          self.fft_forward.update_arrays(self.r, self.f)
          self.fft_forward.execute()
          self.f = fftshift(self.f)
          psi_f_cbed_region_mag = abs(self.f[self.ind_dp_lb:self.ind_dp_ub, self.ind_dp_lb:self.ind_dp_ub])
          psi_f = self.f
          return psi_f, psi_f_cbed_region_mag

      def updateFourierIntensity(self, psi_f, dp, denominator):
          #psi_f: wave function in Fourier space
          self.f[:,:] = psi_f;
          f_cbed = self.f[self.ind_dp_lb:self.ind_dp_ub, self.ind_dp_lb:self.ind_dp_ub]
          
          f_cbed[self.paraDict['badPixels']==0] = f_cbed[self.paraDict['badPixels']==0]/(denominator[self.paraDict['badPixels']==0]+1e-16) * dp[self.paraDict['badPixels']==0]
          self.f[self.ind_dp_lb:self.ind_dp_ub, self.ind_dp_lb:self.ind_dp_ub] = f_cbed
          self.f = ifftshift(self.f)
          if 'filter_f_psi' in self.paraDict: self.f = self.f * self.paraDict['filter_f_psi']
          self.fft_inverse.update_arrays(self.f,self.r)
          self.fft_inverse.execute();          
          return self.r / self.N_tot

      #################### position correction ####################
      def gradPositionCorrection(self, probe, O, ind_dp, delta_psi):
          dx_O, dy_O = self.getObjectGradient(O)
          
          dx_OP = dx_O*probe
          shift_x = np.sum(real(conj(dx_OP)*delta_psi))/np.sum(abs(dx_OP)**2)
          
          dy_OP = dy_O*probe
          shift_y = np.sum(real(conj(dy_OP)*delta_psi))/np.sum(abs(dy_OP)**2)

          #update position
          #print(shift_y)
          self.paraDict['ppY'][ind_dp] = self.paraDict['ppY'][ind_dp] + shift_y*self.dx_y 
          self.paraDict['ppX'][ind_dp] = self.paraDict['ppX'][ind_dp] + shift_x*self.dx_x

          #position = pi + pf = integer + fraction
          py_i = np.round(self.paraDict['ppY'][ind_dp] / self.dx_y)
          self.py_f[ind_dp] = self.paraDict['ppY'][ind_dp] - py_i * self.dx_y 
          px_i = np.round(self.paraDict['ppX'][ind_dp] / self.dx_x)
          self.px_f[ind_dp] = self.paraDict['ppX'][ind_dp] - px_i * self.dx_x

          #calculate ROI indices in the whole fov
          self.ind_x_lb_s[ind_dp] = (px_i - floor(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_x_ub_s[ind_dp] = (px_i + ceil(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_y_lb_s[ind_dp] = (py_i - floor(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_y_ub_s[ind_dp] = (py_i + ceil(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)

      #################### initialization ####################
      def initializeDataError(self):
          if 's' in self.paraDict:
              if self.paraDict['Niter'] > self.paraDict['s'].size:
                  s = self.pad(paraDict['s'].flatten(),(0, self.paraDict['Niter']-self.paraDict['s'].size),'constant')
              else:
                  s = self.paraDict['s'].flatten()
          else:
              s = zeros(self.paraDict['Niter']) #data error
          if 'dp_error_old' in self.paraDict:
              dp_error_old  = self.paraDict['dp_error_old'].flatten()
          else:
              dp_error_old = full(self.paraDict['N_scan'], np.inf)
          return s, dp_error_old

      def initializeScanPositions(self):
          self.center_index_image = int(self.paraDict['N_image']/2)
          self.dx_x = 1.0/(self.paraDict['dk_x'] * self.paraDict['N_roi'])
          self.dx_y = 1.0/(self.paraDict['dk_y'] * self.paraDict['N_roi'])

          #position = pi + pf = integer + fraction
          py_i = np.round(self.paraDict['ppY'] / self.dx_y)
          self.py_f = self.paraDict['ppY'] - py_i * self.dx_y 
          px_i = np.round(self.paraDict['ppX'] / self.dx_x)
          self.px_f = self.paraDict['ppX'] - px_i * self.dx_x

          #calculate ROI indices in entire image
          self.ind_x_lb_s = (px_i - floor(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_x_ub_s = (px_i + ceil(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_y_lb_s = (py_i - floor(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)
          self.ind_y_ub_s = (py_i + ceil(self.paraDict['N_roi']/2.0) + self.center_index_image).astype(np.int)

      def initializeDiffractionPatterns(self):
          center_index_roi = int(self.paraDict['N_roi']/2)
          #calculate dp indices in ROI image
          self.ind_dp_lb = int(-floor(self.paraDict['N_dp']/2.0) + center_index_roi)
          self.ind_dp_ub = int(ceil(self.paraDict['N_dp']/2.0) + center_index_roi)

      def createFourierCoord(self, N):
          kx = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)
          kx = ifftshift(kx)
          [self.kX, self.kY] = meshgrid(kx,kx)
          self.kX = self.kX * self.paraDict['dk_x']
          self.kY = self.kY * self.paraDict['dk_y']

      def printStatus(self, timeLeft, iter, extraMessage=''):
          timeLeftMin, timeLeftSec = divmod(timeLeft, 60)
          timeLeftHour, timeLeftMin = divmod(timeLeftMin, 60)
          print(self.paraDict['printID'] + extraMessage + '-Iter:%d Time remain: %02d:%02d:%02d' %(iter,timeLeftHour,timeLeftMin,timeLeftSec))

      #################### FFT ####################
      def initializeFFTW(self, N):
          self.f = pyfftw.empty_aligned((N,N),dtype='complex128',n=16) 
          self.r = pyfftw.empty_aligned((N,N),dtype='complex128',n=16)
          self.N_tot = N*N
          self.fft_forward = pyfftw.FFTW(self.r, self.f, axes=(0,1))
          self.fft_inverse = pyfftw.FFTW(self.f, self.r, direction='FFTW_BACKWARD', axes=(0,1))

      def shift(self, func, px, py):
          """Shift function via FFT"""
          self.r[:,:] = ifftshift(func)
          self.fft_forward.update_arrays(self.r, self.f)
          self.fft_forward.execute()
          self.f = self.f*exp(-2*pi*1j*px*self.kX)*exp(-2*pi*1j*py*self.kY)
          self.fft_inverse.update_arrays(self.f, self.r)
          self.fft_inverse.execute();
          return fftshift(self.r) / self.N_tot #fix normalization

      ####################  ####################
      def proj(u, v):
          return u * np.vdot(u,v) / np.vdot(u,u)  

      def gramschmidt(V):
          U = np.copy(V)
          for i in range(1, V.shape[0]):
              for j in range(i):
                  U[i,:] -= proj(U[j,:], V[i,:])
          return U

      def getObjectGradient(self, O):
          Ny, Nx = O.shape
          kx = fftshift(linspace(0,Nx-1,Nx)*1.0/Nx-0.5)
          ky = fftshift(linspace(0,Ny-1,Ny)*1.0/Ny-0.5)
          [kX, kY] = meshgrid(kx,ky)

          O_fx = fft(O,axis=1)
          O_fy = fft(O,axis=0)

          O_dx = ifft(O_fx*kX*2j*pi,axis=1)
          O_dy = ifft(O_fy*kY*2j*pi,axis=0)

          return O_dx, O_dy

