% FUNCTION IM = C2IMAGE(A)
% 
% Returns a RGB image of complex array A where
% the phase is mapped to hue, and the amplitude
% is mapped to brightness.

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

function im = c2image(a, varargin)


if ismatrix(a)
    absa = abs(a);
    phasea = angle(a);
    
    % (optional second argument can switch between various plotting modes)
    abs_range = [];
    if nargin==2
        m = varargin{1};
    elseif nargin==3
        m = varargin{1};
        abs_range = varargin{2};
    else
        m = 1;
    end
    
    if isempty(abs_range)
        nabsa = absa/max(max(absa));
    else
        nabsa = (absa - abs_range(1))/(abs_range(2) - abs_range(1));
        nabsa(nabsa < 0) = 0;
        nabsa(nabsa > 1) = 1;
    end
    
    switch m
        case 1
            im_hsv = zeros([size(a) 3]);
            im_hsv(:,:,1) = mod(phasea,2*pi)/(2*pi);
            im_hsv(:,:,2) = 1;
            im_hsv(:,:,3) = nabsa;
            im = hsv2rgb(im_hsv);
        case 2
            im_hsv = ones([size(a) 3]);
            im_hsv(:,:,1) = mod(phasea,2*pi)/(2*pi);
            im_hsv(:,:,2) = nabsa;
            im = hsv2rgb(im_hsv);
    end
elseif ndims(a)==3
    sz = size(a);
    
    im_hsv = zeros([sz 3]);
    im = zeros([sz 3]);
    
    for ii=1:sz(3)
        absa = abs(a(:,:,ii));
        phasea = angle(a(:,:,ii));
        
        % (optional second argument can switch between various plotting modes)
        abs_range = [];
        if nargin==2
            m = varargin{1};
        elseif nargin==3
            m = varargin{1};
            abs_range = varargin{2};
        else
            m = 1;
        end
        
        if isempty(abs_range)
            nabsa = absa/max(max(absa));
        else
            nabsa = (absa - abs_range(1))/(abs_range(2) - abs_range(1));
            nabsa(nabsa < 0) = 0;
            nabsa(nabsa > 1) = 1;
        end
       
        
        switch m
            case 1
                im_hsv(:,:,ii,1) = mod(phasea,2*pi)/(2*pi);
                im_hsv(:,:,ii,2) = 1;
                im_hsv(:,:,ii,3) = nabsa;
            case 2
                im_hsv(:,:,ii,1) = mod(phasea,2*pi)/(2*pi);
                im_hsv(:,:,ii,2) = nabsa;
        end
        im(:,:,ii,:) = hsv2rgb(squeeze(im_hsv(:,:,ii,:)));

    end
    
end
end

