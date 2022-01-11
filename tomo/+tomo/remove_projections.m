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
    %[Nx,Ny,~] = size(stack_object);
    title_extra = {}; 
    for ii = 1:length(which_remove)
        if isfield(par, 'nresidua_per_frame')
            title_extra{end+1} = sprintf(' N residua: %i',par.nresidua_per_frame(which_remove(ii))); 
        end
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
            par.energy(which_remove,:) = [];
            
            try par.nresidua_per_frame(which_remove) = []; end
            try par.subtomos(which_remove) = []; end
            try par.proj_file_names(which_remove) = []; end
            try par.proj_recon_method(which_remove) = []; end
            try par.proj_roi(which_remove) = []; end
            try par.proj_scanNo(which_remove) = []; end
            try par.object_size_orig(:,which_remove) = []; end

            par.num_proj = length(par.scanstomo); 
            
            disp('Done')
        end
    end



end