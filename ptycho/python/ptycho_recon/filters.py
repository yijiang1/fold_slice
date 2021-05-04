from numpy import *
from numpy.fft import *

def disk(N,radius):
    output = ones((N,N))
    x = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)
    [X,Y] = meshgrid(x,x)
    S = sqrt(X**2+Y**2)
    output[S>radius] = 0
    return output

def gaussian(N,innerRadius,outerRadius,sigma):
    x = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)
    [X,Y] = meshgrid(x,x)
    S = sqrt(X**2+Y**2)

    output = exp(-(S-innerRadius)**2/(2*sigma**2))
    output[S<=innerRadius] = 1
    if outerRadius>0:
        output[S>outerRadius] = 0
    return output

def cosine(N,innerRadius,outerRadius):
    x = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)
    [X,Y] = meshgrid(x,x)
    S = sqrt(X**2+Y**2)

    T = (outerRadius - innerRadius)*2.0
    output = (cos(2.0*pi/T*(S-innerRadius)) + 1 ) / 2
    output[S<=innerRadius] = 1
    output[S>outerRadius] = 0

    return output

def cbed(N_roi,N_dp):
    center_index_roi = floor(N_roi/2.0)
    index_dp_lb = int(-floor(N_dp/2.0) + center_index_roi)
    index_dp_ub = int(ceil(N_dp/2.0) + center_index_roi)
    output = zeros((N_roi,N_roi))
    output[index_dp_lb:index_dp_ub,index_dp_lb:index_dp_ub] = 1
    return output

def square(N, halfSideLength):
    x = linspace(-floor(N/2.0),ceil(N/2.0)-1,N)    
    [X,Y] = meshgrid(x,x)
    print(halfSideLength)
    output = ones((N,N))
    output[X>halfSideLength] = 0
    output[X<-halfSideLength] = 0
    output[Y>halfSideLength] = 0
    output[Y<-halfSideLength] = 0
    return output
