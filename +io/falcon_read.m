 
% Read data from Falcon readout electronics
% Input is the filename with path
% Output is a structure containing fields:
%    Data contains one spectrum per measurement point
%    Metadata contains some other information like the number of points per
%    data transfer and the total number of points
% 12 March 2019
%
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function dataout = falcon_read(filename)

data1 = io.HDF.hdf5_load(filename);

spectra_per_transfer    = data1.entry.instrument.FalconX1.PixelsPerBuffer(end);%meta.spectra_per_transfer; % aka pixels per buffer in the MEDM
numberofpositions       = data1.entry.instrument.FalconX1.CurrentPixel(end);%data1.entry.instrument.NDAttributes.CurrentPixel(end);%meta.numberofpositions;
data = data1.entry.data.data;

data = data(1+256:end,1,:);

N = size(data);
if size(data,3) > 1
    data = reshape(data,N(1)/spectra_per_transfer,spectra_per_transfer*N(3));
else
    data = reshape(data,N(1)/spectra_per_transfer,spectra_per_transfer);
end

dataout.data = data(1:2:end,1:numberofpositions);

end




