%   SAVE_MERGED_PROJECTIONS save projections generated from tomography to be loaded as initial guess in
%   ptychography , projections are saved to the same path from where they were loaded 
%   just with a differenent suffix 
% 
%   save_merged_projections(par,stack_object, volData_c, theta, total_shift, name_sufix)
%
% Inputs: 
%   **par           - tomography parameters structure 
%   **stack_object  - array of complex valued projections 
%   **volData_c     - complex valued sample reconstruction 
%   **theta         - measured angles 
%   **total_shift   - reconstructed shift of the projections 
%   **name_sufix    - extra sufix added to the saved projections 
% *returns* 
%   ++prepared_objects - prepared lsit of structures for 3D ptychotomo 

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



function prepared_objects = save_merged_projections(par,stack_object, theta, total_shift, name_sufix)

import ptycho.* 
import utils.* 
import io.*
import plotting.imagesc3D


scanstomo = par.scanstomo; 
num_proj = length(scanstomo); 

verbose(0,'Checking available files')
verbose(0); % make it quiet 
missing_projections = []; 
proj_file_names = cell(length(scanstomo),1);
for num = 1:length(scanstomo)
    progressbar(num, length(scanstomo))
    path = find_projection_files_names(par, scanstomo(num)); 
    if isempty(path)
        missing_projections(end+1) = num; 
        continue
    end
    % take the last file fitting the constraints 
    proj_file_names{num} = path;
end

verbose(par.verbose_level); % return to original settings 
if ~isempty(missing_projections)
   verbose(1,['Did not find following scans:', num2str(missing_projections)]) 
end

[Nx,Ny,~] = size(stack_object);
verbose(1,'Preparing projections')
for num = 1:num_proj
    progressbar(num, num_proj)

    [filepath,name, ext] = fileparts(proj_file_names{num});
    
    out_name = [name,'_',name_sufix]; 
    fullpath = [filepath, '/', out_name,'.mat']; 
    
    % load original data     
    d = io.load_ptycho_recons(proj_file_names{num}, 'recon'); 
    
    
    positions = h5read(proj_file_names{num}, '/reconstruction/p/positions')'; 
    
    obj_size = size(d.object); 
    object = zeros(obj_size,'like',stack_object);
    object(1:min(end,Nx),1:min(end,Ny)) = ...
        stack_object(1:min(end,obj_size(1)),1:min(end,obj_size(2)),num); 
    
    % load the complex projections from fp16 precision if used 
    object = fp16.get(object); 

    
    %% prepare structure for ptychtomo solver 
    
    prepared_objects{num}.object = object; 
    prepared_objects{num}.positions = positions; 
    prepared_objects{num}.illum_sum = []; 
    prepared_objects{num}.weight = []; 
    prepared_objects{num}.angle = theta(num);
    prepared_objects{num}.position_offset = round(total_shift(num,:));
    % testme 
    prepared_objects{num}.probe =  utils.imshift_fft(d.probe, total_shift(num,:) -round(total_shift(num,:))) ; 
    prepared_objects{num}.scan_id = scanstomo(num); 
    prepared_objects{num}.proj_id = num; 
    
    % phase removal in probe 
    [~,~, gamma_x, gamma_y] = utils.stabilize_phase(d.object,object, 'binning', 8);
    
    % remove ramp from probe as well 
    xramp = pi*(linspace(-1,1,par.asize(1)))';
    yramp = pi*(linspace(-1,1,par.asize(2)));
    c_offset = xramp.*gamma_x*par.asize(1) + yramp.*gamma_y*par.asize(2);
    prepared_objects{num}.probe = prepared_objects{num}.probe .* exp(1i*c_offset);
    
   % plotting.imagesc3D( angle(d.object .* conj(object(1:4:end,1:4:end))  ))
    %drawnow 
    
end


end
