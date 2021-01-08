% get_deformation_fields - calculate the DVF from the observed deformation
% field arrays bu deconvolution 
%
% [deform_tensors_linear, inv_deform_tensors_linear] = ...
%                    get_deformation_fields(shift_tensors, regularize_lambda, Npix_vol)
%
% Inputs:
%    **shift_tensors       observed deformation 
%    **regularize_lambda   regularization constant for the deconvolution 
%    **Npix_vol            size of the reconstructed volume  
% Outputs: 
%    ++deform_tensors_linear        calculate forward DVF 
%    ++inv_deform_tensors_linear    calculate inverse DVF 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function [deform_tensors_linear, inv_deform_tensors_linear] = ...
                    get_deformation_fields(shift_tensors, regularize_lambda, Npix_vol)
    % simple deconvolution of the recovered shift arrays to the object
    % defomration arrays for linear deformation model 
    
    
    Nblocks = length(shift_tensors); 
    
    if Nblocks == 1
       error('Number of blocks (subtomos) has to be > 1') 
    end
    
    if Nblocks > 1
        conv_mat = spdiags(ones(Nblocks,1),0,Nblocks, Nblocks+1) + spdiags(ones(Nblocks,1),1,Nblocks, Nblocks+1); 
        conv_mat = conv_mat ./ sum(conv_mat,2);
        regul_mat = spdiags(2*ones(Nblocks,1), 0, Nblocks, Nblocks+1) - spdiags(ones(Nblocks,1), 1, Nblocks, Nblocks+1)-spdiags(ones(Nblocks,1), -1, Nblocks, Nblocks+1);
        regul_mat(1,1:2) = 0;

        deform_mat = [];
        for block = 1:Nblocks
            for ax = 1:3
                deform_mat(:,:,:,ax,block) = gather(shift_tensors{block}{ax}); 
            end
        end
        size_deform_mat = size(deform_mat); 
        deform_mat = reshape(deform_mat, [], Nblocks); 

        %% perform Tikhonov based deconvolution 

%         N = 50;
%         lams = logspace(-5,1,N);
%         for i = 1:N
            deconv_def_mat = ((conv_mat'*conv_mat + regularize_lambda*regul_mat'*regul_mat)\(conv_mat'*deform_mat'))';
            % enforce zero for the first deformation 
            deconv_def_mat = deconv_def_mat - deconv_def_mat(:,1); 

%             err(i) = math.mean2((conv_mat*deconv_def_mat' - deform_mat').^2);
%         end
        

        deconv_def_mat = reshape(deconv_def_mat,[size_deform_mat(1:4), Nblocks+1] ); 
        for block = 1:Nblocks+1
            for ax = 1:3
                deform_tensors{block}{ax} = single(deconv_def_mat(:,:,:,ax,block)); 
            end
        end
    else
       deform_tensors =  shift_tensors; 
    end
    
    [deform_tensors,inv_deform_tensors] =  nonrigid.invert_DVF(deform_tensors, Npix_vol)   ;
    
    % join blocks to keep initial and final deform for each block together
    % -> linear deformation evolution is assumed in between 
    if Nblocks > 1
        for block = 1:Nblocks
            deform_tensors_linear{block} = [deform_tensors{block}; deform_tensors{block+1}]; 
            inv_deform_tensors_linear{block} = [inv_deform_tensors{block}; inv_deform_tensors{block+1}]; 
        end
    else  % just assume one single deformation 
        deform_tensors_linear = deform_tensors; 
        inv_deform_tensors_linear = inv_deform_tensors; 
    end
    

end
