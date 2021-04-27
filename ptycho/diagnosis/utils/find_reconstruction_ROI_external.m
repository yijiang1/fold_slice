% FIND_RECONSTRUCTION_ROI_EXTERNAL precalculate the reconstruction regions 
% Modified by YJ for external use outside the GPU engines
% ROI is consistent with "obj_proj" in LSQML.m
%
% [oROI, oROI_vec, sub_px_shift] = find_reconstruction_ROI2( positions,Np_o, Np_p )
%
% ** positions    Npox*2 vector of scanning positions 
% ** Np_o       object size 
% ** Np_p       probe size 
%
% returns:
% ++ oROI               cell array contaning range for each view 
% ++ oROI_vec           cell array contaning range for each view in vector shape 
% ++ sub_px_shift       subpixel rounding errors, used for subpixel shift 
%

function [oROI, oROI_vec, sub_px_shift] = find_reconstruction_ROI_external( positions,Np_o, Np_p )

    positions = positions(:,[2,1]); 
    positions = positions + ceil(Np_o/2-Np_p/2); 
    sub_px_shift = positions - round(positions); 
    
    sub_px_shift = sub_px_shift(:,[2,1]);  % return to the original XY coordinates 
    
    positions = round(positions);

    range = [min(positions), max(positions)+ Np_p];
            
    if any(range(1:2) < 0) || any(range(3:4) > Np_o)
        error('Object size is too small, not enough space for probes !! \nposition range: %i %i %i %i, \nobject size: %i %i ', range(1), range(2), range(3), range(4), Np_o(1), Np_o(2)) 
    end
    
    oROI = cell(2,1);
    for dim = 1:2
        oROI{dim} = [positions(:,dim),positions(:,dim)+ Np_p(dim)-1];
        oROI{dim} = uint32(oROI{dim}); 
    end
    
    if nargout > 1
        Npos = length(positions);
        oROI_vec = cell(Npos,2);
        for ii = 1:Npos
            for i = 1:2
                oROI_vec{ii,i} = (oROI{i}(ii,1)):(oROI{i}(ii,2));
            end
        end
    end
     
end