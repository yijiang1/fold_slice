% FUNCTION img_stack = imrotate_ax(img_stack, angle, ax, val, method)
% bilinear rotate stack of images along given axis 

% IMROTATE_AX  bilinear rotate stack of images along given axis 
% img_stack = imrotate_fft(img, theta, axis,ax=3 val=0, method='bilinear')
% 
% Inputs: 
%   **img   - stacked array of images to be rotated 
%   **theta - rotation angle 
%   **axis  - rotation axis 
% *optional*
%   **val   - fill missing values by this number 
%   **method- interpolation method
% *returns*:
%   ++img   - rotated image 

    

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


function img_stack = imrotate_ax(img_stack, angle, ax, val, method)       
    if nargin < 4 
        val = 0; 
    end
    if nargin  < 3
        ax = 3;
    end
    if nargin < 5
        method = 'bilinear'; 
    end

    %% bilinear rotation along given axis 
    N = size(img_stack,ax); 
    
    if ~isa(img_stack, 'gpuArray')
        %% for CPU based rotation process the inputs slice by slice 
        ind = {':',':',':'};
        for ii = 1:N
            if utils.verbose > 0; utils.progressbar(ii,N); end
            ind{ax} = ii; 

                    auxslice = squeeze(img_stack(ind{:}));
            if any(size(auxslice) ~= N)
                % pad in case of asymmetric input 
                auxslice_pad =  padarray( auxslice , [N N], val); % val is the background
                auxslice_pad = imrotate(auxslice_pad,angle,method,'crop');
                auxslice = auxslice_pad(N+1:end-N, N+1:end-N);
            else
                auxslice = imrotate(auxslice,angle,method,'crop');
            end
            img_stack(ind{:}) = auxslice;
        end
    else
        %% for GPU call directly the internal code for 2D interpolation and process the image block in one step 
        
        if ax == 1
           img_stack = permute(img_stack, [3,2,1]); angle = - angle; 
        elseif ax == 2
           img_stack = permute(img_stack, [1,3,2]);
        end
    
        outputSize = size(img_stack);
        if isreal(img_stack)
            img_stack = images.internal.gpu.imrotate(img_stack, angle, method, outputSize);
        else
            img_stack = complex(images.internal.gpu.imrotate(real(img_stack), angle, method, outputSize),...
                    images.internal.gpu.imrotate(imag(img_stack), angle, method, outputSize));
        end
        
        if ax == 1
           img_stack = permute(img_stack, [3,2,1]);
        elseif ax == 2
           img_stack = permute(img_stack, [1,3,2]);
        end
    

    
end