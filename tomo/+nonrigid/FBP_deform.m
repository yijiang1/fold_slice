% [rec, rec_blocks] = FBP_deform(sinogram, cfg, vectors, varargin)
% FUNCTION filtered back projection 
% Inputs:
%     sino - sinogram (Nlayers x width x Nangles)
%     cfg - config struct from ASTRA_initialize
%     vectors - vectors of projection rotation generated by ASTRA_initialize
%     varargin - see the code 

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


function [rec, rec_blocks] = FBP_deform(sinogram, cfg, vectors,  block_size, inv_deform_tensors, varargin )

    par = inputParser;
    par.KeepUnmatched = true; 
    par.addOptional('valid_angles', [], @isnumeric)
    par.addOptional('verbose', 1)   % verbose = 0 : quiet, verbose : standard info , verbose = 2: debug 

    par.parse(varargin{:})
    r = par.Results;
        
    Nblocks = length(inv_deform_tensors); 
    Nangles = cfg.iProjAngles; 

    if isempty(r.valid_angles)
        r.valid_angles = 1:Nangles;
    end

    rec = 0; 
    for ll = 1:Nblocks
        if r.verbose ; utils.progressbar(ll, Nblocks); end
        
        ids = 1+(ll-1)*block_size:min(Nangles, ll*block_size); 
        if ~isempty(r.valid_angles)
            ids = intersect(ids, r.valid_angles); 
        end
        if isempty(ids)
            continue
        end
        if isempty(inv_deform_tensors{ll})
            warning('Empty inv_deform_tensors, skipping %i projections', length(ids))
        end
        rec_blocks{ll} = tomo.FBP(sinogram, cfg, vectors,'deformation_fields',inv_deform_tensors{ll}, varargin{:}, 'valid_angles', ids, 'verbose', 0);
        rec  =  rec + rec_blocks{ll} * length(ids) / length(r.valid_angles);
        if nargout == 1
             rec_blocks{ll} = []; % save memory
        end
    end

end



