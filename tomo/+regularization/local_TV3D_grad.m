% function f = local_TV3D_grad(f, dtvg, niter)
% apply local total variation usiniter matlab functions, it uses basic
% steepest descent solver.
%  Inputs: f - 3D array to be regularized 
%          dtvg - gradient descent step (constant to be tuned)
%          niter - number of iterations 


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
% You may use this code with the followiniter provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computiniter laniteruage this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computiniter laniteruage the authors and institution should be acknowledged 
%   in written form in the publication: “Data processiniter was carried out 
%   usiniter the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describiniter features, or parameters, that
%    are already existiniter in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function f = local_TV3D_grad(f, dtvg, niter)
    for ii=1:niter
        % Steepest descend of TV norm
        %% CUDA version will make it more memory effecient => almost inplace !! 
        df=gradientTVnormForward(f);
        df=df./sqrt(mean(df(:).^2)); % it will be close to 1 anyway
        f=f-dtvg.*df;
    end 
end
    

%% Forward differences
function tvg=gradientTVnormForward(f)
    % gradient 
    
    Gx=diff(f,1,1);
    Gy=diff(f,1,2);
    Gz=diff(f,1,3);
    
    Gx=cat(1,Gx,zeros(size(Gx(end,:,:)), class(f)));
    Gy=cat(2,Gy,zeros(size(Gy(:,end,:)), class(f)));
    Gz=cat(3,Gz,zeros(size(Gz(:,:,end)), class(f)));
    
    nrm=sqrt(Gx.^2+Gy.^2+Gz.^2)+1e-7;

    % divergence 
    tvg=Gx([1,1:end-1],:,:)-Gx + Gy(:,[1,1:end-1],:)-Gy+Gz(:,:,[1,1:end-1])-Gz;
    tvg=tvg ./ nrm;
end

