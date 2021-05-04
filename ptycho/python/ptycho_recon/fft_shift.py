from numpy import *
from numpy.fft import *
import pyfftw

class shift_fftw:
      """Shift function via FFT"""
      def __init__(self, dk_x, dk_y, N):
          kx = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)
          kx = ifftshift(kx)
          [self.kX, self.kY] = meshgrid(kx,kx)
          self.kX = self.kX*dk_x
          self.kY = self.kY*dk_y
          self.f = pyfftw.empty_aligned((N,N),dtype='complex128',n=16) 
          self.r = pyfftw.empty_aligned((N,N),dtype='complex128',n=16)
          self.N_tot = N*N
          self.fft_forward = pyfftw.FFTW(self.r, self.f, axes=(0,1))
          self.fft_inverse = pyfftw.FFTW(self.f, self.r, direction='FFTW_BACKWARD', axes=(0,1))

      def shift(self, func, px, py):
          self.r[:,:] = ifftshift(func)
          self.fft_forward.update_arrays(self.r, self.f)
          self.fft_forward.execute()
          self.f = self.f*exp(-2*pi*1j*px*self.kX)*exp(-2*pi*1j*py*self.kY)
          self.fft_inverse.update_arrays(self.f, self.r)
          self.fft_inverse.execute();
          return fftshift(self.r) / self.N_tot #fix normalization
