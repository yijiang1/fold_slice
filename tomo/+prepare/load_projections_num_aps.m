%   LOAD_PROJECTIONS_NUM_APS load reconstructed projections numbers
%  created by YJ based on PSI's function
%   Only read scan numbers, usedful for debugging
%   Inputs: 
%       **par - parameter structure 
%       **exclude_scans - list of scans to be excluded from loading, [] = none
%       **theta - angles of the scans 

function [scanstomo]  = load_projections_num_aps(par, exclude_scans, theta)

import ptycho.* 
import utils.* 
import io.*
import plotting.*

scanstomo = par.scanstomo; 

% avoid loading scans listed in 'exclude_scans'
if  ~isempty(exclude_scans)
    ind = ismember(scanstomo, exclude_scans); 
    scanstomo(ind) = []; 
    theta(ind) = []; 
end

verbose(1,'Checking available files')
missing_scans = []; 
proj_file_names = {};
proj_recon_method = {};
proj_roi = {};
proj_scanNo = {};
for num = 1:length(scanstomo)
    progressbar(num, length(scanstomo))
    %proj_file_names{num} = find_ptycho_filename(par.analysis_path,scanstomo(num),par.fileprefix,par.filesuffix, par.file_extension);
    %proj_file_names{num} = find_projection_files_names_aps(par, scanstomo(num));
    [proj_file_names{num},proj_recon_method{num},proj_roi{num},proj_scanNo{num}] = find_ML_recon_files_names(par, scanstomo(num));
    %disp(proj_file_names{num})
    if isempty(proj_file_names{num})
        missing_scans(end+1) = scanstomo(num); 
    end
end

verbose(par.verbose_level); % return to original settings

if ~isempty(missing_scans)
    ind = ismember(scanstomo, missing_scans); 
    verbose(1,['Scans  not found are ' num2str(missing_scans)])
    verbose(1,['Projections not found are ' num2str(find(ind))])
    scanstomo(ind) = []; 
    theta(ind) = []; 
    proj_file_names(ind) = []; 
    proj_recon_method(ind) = [];
    proj_roi(ind) = [];
    proj_scanNo(ind) = [];
else
    verbose(1,'All projections found')
end

num_proj = length(scanstomo); 


tic

if num_proj == 0
    verbose(0, 'No new projections loaded')
    return
end
    
end
