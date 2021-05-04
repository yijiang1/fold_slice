import numpy as np
import scipy.io as sio #for read/write matlab file
from scipy import ndimage

class STEMprobe:
      """Probe function for STEM"""
      def __init__(self):
          self.dx = 1.0
          self.Nside = 256
          self.px = 0
          self.py = 0
          self.voltage = 300 #keV
          self.alpha_max = 30 #mrad
          self.df = 0 #angstrom
          self.cs = 0 #mm
          self.f_a2 = 0 #angstrom
          self.theta_a2 = 0
          self.f_a3 = 0 #angstrom
          self.theta_a3 = 0
          self.f_c3 = 0 #angstrom
          self.theta_c3 = 0
          self.Fourier_mag = 1

      def printParameters(self):
          self.wavelength = 12.398/np.sqrt((2*511.0+self.voltage)*self.voltage)  #angstrom
          #print out all the parameters in the dm reconstruction
          print("probe size:",self.Nside,"x",self.Nside)
          print("distance between adjacent pixels: dx =", self.dx)
          print("distance between adjacent pixels in Fourier space: dk =", 1.0/(self.dx*self.Nside))

          print("probe position: (px,py) =(",self.px,",",self.py,")")
          print("beam voltage =", self.voltage, "keV")
          print("beam wavelength =",self.wavelength, "angstrom")
          print("semi-convergence angle =", self.alpha_max, "mrad")
          print("defocus=", self.df, "angstrom")
          print("spherical aberration =", self.cs, "angstrom")
          print("two-fold astigmatism =", self.f_a2, "angstrom.", "azimuthal orientation=",self.theta_a2, "rad")
          print("three-fold astigmatism =", self.f_a3, "angstrom.", "azimuthal orientation=",self.theta_a3, "rad")
          print("coma =", self.f_c3, "angstrom.", "azimuthal orientation=",self.theta_c3, "rad")

      def generateProbe(self):
          print("generating probe function...")
          self.wavelength = 12.398/np.sqrt((2*511.0+self.voltage)*self.voltage)  #angstrom
          amax = self.alpha_max*1e-3 # in rad
          amin = 0.0

          k_max = amax/self.wavelength
          k_min = amin/self.wavelength

          dk= 1.0/(self.dx*self.Nside)
          kx = np.linspace(-np.floor(self.Nside/2.0),np.ceil(self.Nside/2.0)-1,self.Nside)
          [kY,kX] = np.meshgrid(kx,kx)
          kX = kX*dk; kY = kY*dk;
          kR = np.sqrt(kX**2+kY**2)
          theta = np.arctan2(kY,kX)
          
          chi = -np.pi*self.wavelength*kR**2*self.df + np.pi/2*self.cs*1e7*self.wavelength**3*kR**4+np.pi*self.f_a2*self.wavelength*kR**2*np.sin(2*(theta-self.theta_a2))+2*np.pi/3*self.f_a3*self.wavelength**2*kR**3*np.sin(3*(theta-self.theta_a3))+2*np.pi/3*self.f_c3*self.wavelength**2*kR**3*np.sin(theta-self.theta_c3)

          probe = np.exp(-1j*chi)*np.exp(-2*np.pi*1j*self.px*kX)*np.exp(-2*np.pi*1j*self.py*kY)
          probe[kR>k_max] = 0
          probe[kR<k_min] = 0
          
          if self.Fourier_mag != 1:
              probe = probe/abs(probe) * self.Fourier_mag

          #probe = probe/np.sum((np.abs(probe)**2)) #normalize probe
          
          probe = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(probe)))
          probe = probe/np.sqrt(np.sum((np.abs(probe)**2)*self.dx*self.dx)) #normalize probe
          #probe = probe/np.sqrt(np.sum((np.abs(probe)**2))) #normalize probe
          
          mask = np.ones(kR.shape)
          mask[kR>k_max] = 0
          mask[kR<k_min] = 0

          return probe
