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
          self.dp = dp
          self.initialProbe = initialProbe
          self.paraDict = {'dk':dk}
          self.paraDict['dk_y'] = dk
          self.paraDict['dk_x'] = dk
          self.paraDict['ppX'] = ppX
          self.paraDict['ppY'] = ppY
          self.paraDict['N_dp'] = dp.shape[1]
          self.paraDict['N_scan'] = dp.shape[0]
          self.paraDict['badPixels'] = zeros((dp.shape[1],dp.shape[1]))
          self.paraDict['Niter'] = 200
          self.paraDict['Niter_save'] = 50
          self.paraDict['Niter_print'] = 1
          self.paraDict['beta'] = 1.0
          self.paraDict['alpha'] = 1.0
          
          self.paraDict['Niter_update_probe'] = 10
          self.paraDict['uniformInitialObject'] = False
          self.paraDict['normalizeInitialProbe'] = False

          self.paraDict['filter_r_type_psi'] = 'none'
          self.paraDict['filter_r_type_probe'] = 'none'
          self.paraDict['filter_f_type_psi'] = 'cbed'
          self.paraDict['filter_f_type_probe'] = 'none'
         
          self.paraDict['saveData'] = False
          self.paraDict['loadData'] = False
          
          self.paraDict['reconID'] = 0
          self.paraDict['printID'] = ''
          
          #position correction
          self.paraDict['offsetMode'] = ''
          self.paraDict['errorMetric'] = ''
      
      def recon(self):
          print("begin ptychographic reconstruction")
          pie_mixed_states.reconPIE_mixed_state(self.dp, self.paraDict)

      def initialize(self, result_dir):
           ################################## probe #################################
           #create initial probe
           self.initialProbe.dx = 1.0/(self.paraDict['dk']*self.paraDict['N_roi'])
           #print self.initialProbe.dx
           self.initialProbe.Nside = self.paraDict['N_roi']

           if not 'probe0' in self.paraDict: self.paraDict['probe0'] = self.initialProbe.generateProbe()
           #self.initialProbe.printParameters()
           result_dir = result_dir + "_dk"+str(np.round(self.paraDict['dk'],4))
           if self.initialProbe.df != 0: result_dir = result_dir + "_df"+str(self.initialProbe.df) + "A"
           if self.initialProbe.cs != 0: result_dir = result_dir + "_cs"+str(self.initialProbe.cs) + "mm"

           ################################## save data #################################
           if self.paraDict['saveData']:
               print('saving dp_recon...')
               if 'dataDir' in self.paraDict:
                   if not os.path.exists(self.paraDict['dataDir']): os.makedirs(self.paraDict['dataDir'])
                   os.chdir(self.paraDict['dataDir'])
               else:
                   if not os.path.exists(result_dir): os.makedirs(result_dir)
                   os.chdir(result_dir)
               sio.savemat('dp_recon',{'dp_recon':self.dp,'dk':self.paraDict['dk']})

           ################################## filters #################################
           result_dir = self.generateFilters(result_dir, self.paraDict['filter_r_type_psi'], 'r', 'psi')
           result_dir = self.generateFilters(result_dir, self.paraDict['filter_r_type_probe'], 'r', 'probe')
           result_dir = self.generateFilters(result_dir, self.paraDict['filter_f_type_psi'], 'f', 'psi')
           result_dir = self.generateFilters(result_dir, self.paraDict['filter_f_type_probe'], 'f', 'probe')

           result_dir  = result_dir + "/NiterUpdateProbe" + str(self.paraDict['Niter_update_probe'])

           ################################# miscellaneous ############################### 
           if self.paraDict['normalizeInitialProbe']:
               result_dir = result_dir + "_normIniProbe"
           if np.any(self.paraDict['badPixels']): result_dir = result_dir + "_excludeBadPixels"
           if 'probeMask' in self.paraDict: result_dir = result_dir + "_probeMask"

           if 'probe_profile' in self.paraDict: result_dir = result_dir + "_imposeProbeProfile"

           if self.paraDict['Niter_update_position'] < self.paraDict['Niter']: 
               result_dir = result_dir + "_position_correction"+str(self.paraDict['Niter_update_position'])

           self.paraDict['saveName'] =  "recon_N_roi"+str(self.paraDict['N_roi'])

           ################################# scan position correction ###############################
           if self.paraDict['offsetMode'] == 'random':
               dirPositionCorrection = "/scanPositionCorrection_Niter"+str(self.paraDict['N_correct_sp'])+"_Noffset"+str(self.paraDict['N_offset'])+ "_maxOffset"+str(self.paraDict['maxOffset'])+"_totalShiftUB"+str(self.paraDict['totalShiftUpperBound'])+"_maxOffsetLB"+str(self.paraDict['maxOffsetLowerBound'])
           elif self.paraDict['offsetMode'] == 'uniformDirection':
               dirPositionCorrection = "/scanPositionCorrection_Niter"+str(self.paraDict['N_correct_sp'])+"_Noffset"+str(self.paraDict['N_offset'])+ "_maxOffset"+str(self.paraDict['maxOffset'])+"_uniformDirection"
           else:
               dirPositionCorrection = ''


           ############################### mixed states ptychography ######################################
           sult_dir = result_dir + "/N_mixed_states" + str(self.paraDict['Niter_update_states'])
           if self.paraDict['N_probe'] > 1:
               result_dir = result_dir + "_N_prob" + str(self.paraDict['N_probe'])
           
           ####################################### create result dir ######################################
           if not os.path.exists(result_dir): os.makedirs(result_dir)
           os.chdir(result_dir)
           return result_dir

      def generateFilters(self, result_dir, filterType, space, waveFunction):
          if self.paraDict['filter_' + space + '_type_' + waveFunction] == 'none':
              return result_dir
          print("Generating " + filterType + " filter in " + space + " space for " + waveFunction)
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

