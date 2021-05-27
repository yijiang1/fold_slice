from numpy import *
from numpy.fft import *
import numpy as np
import filters as filters
import scipy.io as sio #for read/write matlab file
from scipy import ndimage
import os #for change directory

import utility_function as u
import pie_mixed_states
from probe import STEMprobe

##############################################################################################
class ptycho:
      """Ptychography reconstruction"""
      def __init__(self, dp, dk, initialProbe, ppX, ppY):
          self.initialProbe = initialProbe
          self.paraDict = {'dk':dk}

          ########## reshape scan positions and diffraction patterns ##########
          N_scan_y = dp.shape[2]
          N_scan_x = dp.shape[3]
          N_scan_tot = N_scan_y * N_scan_x

          self.dp = zeros((N_scan_tot,dp.shape[0],dp.shape[1]))
          self.paraDict['ppX'] = ppX.reshape(N_scan_tot)
          self.paraDict['ppY'] = ppY.reshape(N_scan_tot)

          for i in range(N_scan_y):
              for j in range(N_scan_x):
                  index = i*N_scan_x + j
                  self.dp[index,:,:] = sqrt(dp[:,:,i,j])

          self.paraDict['dk_y'] = dk
          self.paraDict['dk_x'] = dk
          
          self.paraDict['N_dp'] = self.dp.shape[1]
          self.paraDict['N_scan'] = self.dp.shape[0]
          self.paraDict['badPixels'] = zeros((self.dp.shape[1],self.dp.shape[1]))
          self.paraDict['Niter'] = 200
          self.paraDict['Niter_save'] = 50
          self.paraDict['Niter_print'] = 1
          self.paraDict['beta'] = 1.0
          self.paraDict['alpha'] = 0.1
          
          self.paraDict['Niter_update_probe'] = 10
          self.paraDict['uniformInitialObject'] = True
          self.paraDict['normalizeInitialProbe'] = True

          self.paraDict['filter_r_type_psi'] = 'none'
          self.paraDict['filter_r_type_probe'] = 'none'
          self.paraDict['filter_f_type_psi'] = 'cbed'
          self.paraDict['filter_f_type_probe'] = 'none'
         
          self.paraDict['saveData'] = False
          self.paraDict['loadData'] = False
          
          self.paraDict['reconID'] = 0
          self.paraDict['printID'] = ''
          
          #mixed-states
          self.paraDict['N_probe'] = 1
          self.paraDict['N_object'] = 1
      
      def recon(self):
          print("begin ptychographic reconstruction")
          pie_mixed_states.reconPIE_mixed_state(self.dp, self.paraDict)

      def initialize(self, result_dir):
           ########## initial probe ##########
           #create initial probe
           self.initialProbe.dx = 1.0/(self.paraDict['dk']*self.paraDict['N_roi'])
           #print self.initialProbe.dx
           self.initialProbe.Nside = self.paraDict['N_roi']

           if not 'probe0' in self.paraDict: self.paraDict['probe0'] = self.initialProbe.generateProbe()
           #self.initialProbe.printParameters()

           ########## save data ##########
           if self.paraDict['saveData']:
               print('saving dp_recon...')
               if 'dataDir' in self.paraDict:
                   if not os.path.exists(self.paraDict['dataDir']): os.makedirs(self.paraDict['dataDir'])
                   os.chdir(self.paraDict['dataDir'])
               else:
                   if not os.path.exists(result_dir): os.makedirs(result_dir)
                   os.chdir(result_dir)
               sio.savemat('dp_recon',{'dp_recon':self.dp,'dk':self.paraDict['dk']})

           ########## filters ##########
           a = self.generateFilters('', self.paraDict['filter_r_type_psi'], 'r', 'psi')
           a = self.generateFilters('', self.paraDict['filter_r_type_probe'], 'r', 'probe')
           a = self.generateFilters('', self.paraDict['filter_f_type_psi'], 'f', 'psi')
           a = self.generateFilters('', self.paraDict['filter_f_type_probe'], 'f', 'probe')

           self.paraDict['saveName'] =  "recon"

           ########## create result dir ##########
           if not os.path.exists(result_dir): os.makedirs(result_dir)
           os.chdir(result_dir)
           return result_dir

      def generateFilters(self, result_dir, filterType, space, waveFunction):
          if self.paraDict['filter_' + space + '_type_' + waveFunction] == 'none':
              return result_dir
          #print("Generating " + filterType + " filter in " + space + " space for " + waveFunction)
          N_roi = self.paraDict['N_roi']
          result_dir = result_dir + "/filter_" + space + "_"
          filerKey = 'filter_' + space + '_' + waveFunction
          if filterType == "cbed":
               self.paraDict[filerKey]  = filters.cbed(N_roi, self.paraDict['N_dp'] )
               result_dir = result_dir + "cbed" + str(self.paraDict['N_dp'])

          elif filterType == "square":
               cutoff = self.paraDict['filter_' + space + '_inner_cutoff_' + waveFunction]
               self.paraDict[filerKey]  = filters.square(N_roi, cutoff)
               result_dir = result_dir + "square_cutoff"+str(cutoff)

          elif filterType == "gaussian":
               inner_cutoff = self.paraDict['filter_' + space + '_inner_cutoff_' + waveFunction]
               outer_cutoff = self.paraDict['filter_' + space + '_outer_cutoff_' + waveFunction]
               sigma = self.paraDict['filter_' + space +'_gaussian_sigma_' + waveFunction]
               self.paraDict[filerKey] = filters.gaussian(N_roi, inner_cutoff, outer_cutoff, sigma)
               result_dir = result_dir + "gaussian_cutoff"+str(inner_cutoff)+"_outer_cutoff"+str(outer_cutoff)+"_sigma"+str(sigma)

          elif filterType  == "cosine":
               inner_cutoff = self.paraDict['filter_' + space + '_inner_cutoff_' + waveFunction]
               outer_cutoff = self.paraDict['filter_' + space + '_outer_cutoff_' + waveFunction]
               self.paraDict[filerKey] = filters.cosine(N_roi, inner_cutoff, outer_cutoff)
               result_dir = result_dir + "cosine_inner_cutoff"+str(inner_cutoff)+"_outer_cutoff"+str(ff_outer_cutoff)

          elif filterType == "disk":
               cutoff = self.paraDict['filter_' + space + '_disk_cutoff_' + waveFunction]
               self.paraDict[filerKey] = filters.disk(N_roi, cutoff)
               result_dir = result_dir + "disk_cutoff"+str(cutoff)
          else:
              raise RuntimeError('Unknown filter type!')
          if space == 'f':
              self.paraDict[filerKey] = ifftshift( self.paraDict[filerKey] )

          result_dir = result_dir + "_" + waveFunction
          return result_dir

