%im_output = remove_linear_phase_smart(im_input, [,<name>,<value>] ...])
%
% Receives a complex valued array, returns the complex array after removing
% the linear and the constant phase offset. Different options are given for
% defining a reference area on the sample
%
% im_input  Input complex valued image
% 
% The optional <name>,<value> pairs are:
%
% 'mask', maskarray     Binary array with ones where the linear phase should be
%                       removed. Alternatively the mask can be an array of
%                       non unitary weights for the computation of average
%                       phase ramp and average phase. For example, the
%                       magnitude can be use to weight 
%
% June 21,2016


% Upcoming features
% 'roi', []             
% method
% upsamp    Linear phase will be removed within 2*pi/upsamp peak to valley
%           in radians
% errorm    Optional ouput with all outputs from dftregistration
%
% Many lines taken from rmsphaseramp.m
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

function [im_output, ph_err] = remove_linear_phase_smart(im, varargin)
import utils.auto_mask_find
import utils.dftregistration

% Defaults
mask = ones(size(im));

% parse the variable input arguments not handled by auto_mask_find
vararg = cell(0,0);                    
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch lower(name)
        case 'mask'
            mask = value;
        otherwise
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end

if any(mask(:) == 1)
    ph = exp(1i*angle(im));
    [gx, gy] = gradient(ph);
    gx = -real(1i*gx./ph);
    gy = -real(1i*gy./ph);
    
    nrm = sum(sum(mask));
    agx = sum(sum(gx.*mask)) / nrm;
    agy = sum(sum(gy.*mask)) / nrm;
    
    sz = size(ph);
    [xx,yy] = meshgrid(1:sz(2),1:sz(1));
    ph_corr = ph.*exp(-1i*(agx*xx + agy*yy));
    ph_corr = ph_corr.*conj(sum(sum(ph_corr.*mask)) / nrm);
    
    im_output = abs(im).*ph_corr;
    
    if nargout > 1
        ph_err = sum(sum(angle(ph_corr).^2.*mask)) / nrm;
    end
else
    warning('The mask provided has only zeros. Cannot remove ramp')
    im_output = im;
    ph_err = 0;
end

end