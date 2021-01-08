%  REMOVE_PROJECTIONS remove projections and adjust the other relevant paramaters 
% 
%  [stack_object,theta,total_shift,par] = remove_projections(stack_object,theta,total_shift,par, which_remove, plot_fnct = @(x)x)
%
% Inputs: 
%  **stack_object        measured projections    
%  **theta               measured angles 
%  **total_shift         Nx2 vector or projection shifts 
%  **par                 parameter structure 
%  **which_remove        indices or logical array denoting the projection to be removed 
%  **plot_fnct = @(x)x   function to be used for plotting, default == @(x)x
% *returns*
%   ++stack_object,theta,total_shift,par - inputs after projection removal 


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


function [stack_object,theta,total_shift,par] = remove_projections(stack_object,theta,total_shift,par, which_remove, plot_fnct, object_ROI)

    if nargin < 6 
        plot_fnct = @(x)x; 
    end

    if islogical(which_remove)
        which_remove = find(which_remove); 
    end
    
    if isempty(which_remove)
       return 
    end
    
    utils.verbose(1, 'Removing %i/%i projections', length(which_remove), par.num_proj)
    [Nx,Ny,~] = size(stack_object);
    title_extra = {}; 
    for ii = 1:length(which_remove)
        title_extra{end+1} = sprintf(' N residua: %i',par.nresidua_per_frame(which_remove(ii))); 
    end


    %%% Getting rid of unwanted projections %%%
    if ~isempty(which_remove)
        
        tomo.show_projections(stack_object(:,:,which_remove), theta(which_remove), par, ...
            'title', 'Projection to be removed','plot_residua', true, 'title_extra', title_extra, 'fnct', plot_fnct, ...
            'rectangle_pos', [object_ROI{2}(1), object_ROI{2}(end), object_ROI{1}(1), object_ROI{1}(end)]) 

        if strcmpi(input(sprintf('Do you want remove %i missing/wrong projections and keep going (y/N)?',length(which_remove)),'s'),'y') 
            disp('Removing missing/wrong projections. stack_object, scanstomo, theta and num_proj are modified')

            stack_object(:,:,which_remove) = [];
            theta(which_remove)=[];
            total_shift(which_remove,:) = [];
            
            par.scanstomo(which_remove)=[];
            try par.energy(which_remove,:) = []; end
            try par.nresidua_per_frame(which_remove) = []; end
            try par.subtomos(which_remove) = []; end
            par.num_proj = length(par.scanstomo); 
            
            disp('Done')
        end
    end



end