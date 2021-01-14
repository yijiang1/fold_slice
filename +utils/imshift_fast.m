%  IMSHIFT_FAST shift of stack of images by given number of
%  pixels and if needed crop / pad image to fit into Npix_new
%
%  img_new=imshift_fast(img_0, x,y, Npix_new=[], type='linear', default_val=0)
%
% Inputs:
%   **img_0       - stack of images 
%   **x,y         - horizontal / vertical shift in pixels  (scalars)
%   **Npix_new    - empty/missing => keep original size, 2x1 vector => embed new image into given frame size 
%   **type        - linear / nearest neighbor interpolation , (missing/empty => linear)
%   **default_val - default value to fill empty regions created after the image shift 
% returns: 
%   ++ img_new    - shifted image padded to Npix_new


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


function img_new=imshift_fast(img_0, x,y, Npix_new, type, default_val)

    
    if nargin < 4 || isempty(Npix_new)
        Npix_new = size(img_0);
    end
    Npix_new = Npix_new(1:2);
    if nargin < 5
        type = 'linear';
    end
    if nargin < 6
        default_val = 0;
    end
    if length(x) > 1 || length(y) > 1
        error('Only scalar position shifts are accepted')
    end
    
    [Nx, Ny, Nimgs] = size(img_0);


    if x==0 && y == 0 && all([Nx,Ny] == Npix_new)
        %% no change is needed, return original image
        img_new = img_0;
        return
    end
    

    pos = -[x,y];
    if strcmp(type, 'linear') && any(round([x,y]) ~= [x,y]) &&  ~isa(img_0, 'logical')
        shift = make_shift(pos - round(pos));
        
        Npix_tmp = size(img_0);
        Npix_tmp(1:2) = Npix_tmp(1:2) + 2;
        img_tmp = zeros(Npix_tmp,  'like', img_0 );
        for i = 1:Nimgs
            img_tmp(:,:,i) = conv2(img_0(:,:,i), shift, 'full');
        end
        img_0 = img_tmp;
    end
    Npix = [Nx,Ny];

    [oROI, pROI] = find_ROI( pos ,Npix_new, Npix );
    

    if all(x==0) && all(y == 0) && all(Npix_new < Npix) 
        img_new = img_0(pROI{:},:);
    else
        if all(abs(pos) <= 1) && all( Npix_new == Npix)
            img_new = img_0;  % for tiny shift reuse the original array 
        else
            img_new = ones([Npix_new,Nimgs], 'like', img_0 )*default_val;
        end
        img_new(oROI{:}, :) = img_0(pROI{:}, :);
    end
end

function [oROI, pROI, oROI_, pROI_] = find_ROI( position, Nobj_new, Nobj_0 )

    oROI = cell(2,1);
    pROI = cell(2,1);
    
    pos = round(position([2,1]));
    %% correction for odd size of the Nobj_new 
    pos = pos - mod(Nobj_0-Nobj_new,2) .* (Nobj_new > Nobj_0); 
    oROI_ = zeros(2);
    pROI_ = zeros(2);
    for dim = 1:2
        range_0 = round(pos(dim) + [1,Nobj_0(dim)] - Nobj_0(dim)/2 + Nobj_new(dim)/2);
        oROI_(dim,:) = min(max(1,range_0), Nobj_new(dim));
        l = oROI_(dim,2) - oROI_(dim,1) +1;
        p1 = min(Nobj_0(dim), Nobj_0(dim) - (range_0(2) - Nobj_new(dim)));
        pROI_(dim,1) = p1 - l+1;
        pROI_(dim,2) = p1;
        if pROI_(dim,1) < pROI_(dim,2)
            pROI{dim} = pROI_(dim,1):pROI_(dim,2);
        else
            pROI{dim} = [];
        end
        if oROI_(dim,1) < oROI_(dim,2)
            oROI{dim} = oROI_(dim,1):oROI_(dim,2);
        else
            oROI{dim} = [];
        end
    end
end


function shift_mat  = make_shift(shift)
    
    x = shift(2);  % correction on pixel position 
    y = shift(1);
    N =  max(1, ceil(abs([x,y])));
    x = x+N(1);
    y = y+N(2);
    dx = x-floor(x);
    dy = y-floor(y);
    
        
    w(1) = dx * dy;
    w(2) = (1-dx) * dy;
    w(3) = dx * (1-dy);
    w(4) = (1-dx) * (1-dy);
 
    ix = 1+floor(x);
    iy = 1+floor(y);

    shift_mat = zeros(2*N+1);
    shift_mat(ix, iy) = w(4);
    shift_mat(ix+1, iy) = w(3);
    shift_mat(ix, iy+1) = w(2);
    shift_mat(ix+1, iy+1) = w(1);    
end
