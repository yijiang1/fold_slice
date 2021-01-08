% show_sinograms - compare reconstructed and mesaured projections
%
%   show_sinograms(rec_avg,dphase, vectors,cfg, angles,  shift_3D_total, par, show_derivative = false )
%
% Inputs:
%    **rec_avg - optimal NCT reconstruction       
%    **dphase    phase difference calculated from the measured data 
%    **vectors   ASTRA configuration vectors
%    **cfg       ASTRA configuration structure
%    **angles    angles of the projections 
%    **shift_3D_total  recovered deformation vector field 
%    **par       parameter structure of the nonrigid tomo 
%  *optional*
%    ++show_derivative  if false, show directly the recovered signal otherwise the phase derivative  

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function show_sinograms(rec_avg,dphase, vectors,cfg, angles,  shift_3D_total, par, show_derivative )
    if nargin < 8
        show_derivative = false;  % show directly the recovered signal 
    end

    Nangles = cfg.iProjAngles;
    Nblocks = length(shift_3D_total); 
    Bsize = ceil(Nangles/Nblocks); 

    % generate deformation tensors from the shift tensors 
    deform_tensors = nonrigid.get_deformation_fields(shift_3D_total, par.regularize_deform_evol, size(rec_avg));

    

    % SHOW ANIMATION OF PROJECTIONS 
    rec_avg = gather(rec_avg);

    split = astra.ASTRA_find_optimal_split(cfg,1,Nblocks); 
        
    for ll = 1:Nblocks
        ids = 1+(ll-1)*Bsize:min(Nangles, ll*Bsize); 
        sinogram_corr(:,:,ids) = gather(tomo.Ax_sup_partial(rec_avg,cfg, vectors(ids,:),split, 'deformation_fields', deform_tensors{ll}));
    end

    split = astra.ASTRA_find_optimal_split(cfg); 
    sinogram_ideal = tomo.Ax_sup_partial(rec_avg,cfg, vectors,split); 
    
    if show_derivative
        dphase_avg = math.get_phase_gradient_1D(exp(-1i*sinogram_ideal),2);
        dphase_corr = math.get_phase_gradient_1D(exp(-1i*sinogram_corr),2);
        sino_diff =  dphase_corr - math.sum2(dphase_corr .* dphase) ./ math.sum2(dphase.^2) .* dphase ; 
%         sino_diff =  dphase_corr - dphase ; 
        dsino_range = math.sp_quantile(sino_diff,[0.001, 0.999],10);
        sino_diff = max(min(sino_diff, dsino_range(2)), dsino_range(1));
        sino_range = math.sp_quantile(dphase_avg,[0.001,0.999],10);
        sino_diff = (sino_diff-dsino_range(1))/diff(dsino_range)*sino_range(2);
        sino_all=  cat(2,dphase_corr, dphase, sino_diff);  
    else
        phase = -math.unwrap2D_fft(dphase, 2, par.air_gap); 
        sino_diff =  sinogram_corr - math.sum2(sinogram_corr .* phase) ./ math.sum2(phase.^2) .* phase ; 
        sino_range = math.sp_quantile(sinogram_corr,[0.001,0.999],10);
        sino_all=  cat(2,sinogram_corr,phase, sino_diff+sino_range(2)/2);  
    end

    [~,order] = sort(angles);
    figure(45864)
    plotting.imagesc3D(sino_all, 'order', order)
    axis off image xy
    colormap bone 
    set(gca, 'clim', sino_range)
    str = 'Model / Data / Difference'; 
    if show_derivative
        str = [str, ' - showing phase-gradient'];
    else
        str = [str, ' - showing unwrapped phase'];
    end
    title(str)
    drawnow 
    
    
    
end