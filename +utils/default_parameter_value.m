%[default_value] = default_parameter_value(mfile_name,parameter_name,vararg)
%  identify the current system to set useful default parameters
    
% Filename: $RCSfile: default_parameter_value.m,v $
%
% $Revision: 1.10 $  $Date: 2011/08/13 14:10:58 $
% $Author:  $
% $Tag: $
%
% Description:
%  identify the current system to set useful default parameters
%
% Note:
% none
%
% Dependencies:
% none
%
%
% history:
%
% April 28th 2010: 
% add plot_radial_integ, find_files, radial_integ
%
% April 2009: 1st version

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

function [default_value] = ...
    default_parameter_value(mfile_name,parameter_name,vararg)
import utils.identify_system

sys_id = identify_system();

% disable the use of the find command for non Unix/Linux based systems
if (strcmp(parameter_name,'UseFind'))
    if (strcmp(sys_id,'Windows'))
        default_value = 0;
    else
        default_value = 1;
    end
end

switch mfile_name
    case 'image_show'
        switch parameter_name
            case 'FigNo' 
                default_value = 1;
            case 'FigClear'
                default_value = 1;
            case 'ImageHandle'
                default_value = 0;
            case 'AutoScale'
                switch sys_id
                    case {'DPC lab', 'mDPC lab', 'cSAXS-mobile'}
                        default_value = [1 1];
                    otherwise
                        default_value = [0 0];
                end
            case 'AxisMin' 
                default_value = 1;
            case 'AxisMax' 
                default_value = 1e5;
            case 'HistScale' 
                default_value = [ 0.15 0.85 ];
            case 'LogScale' 
                switch sys_id
                    case {'DPC lab', 'mDPC lab', 'cSAXS-mobile'}
                        default_value = 0;
                    otherwise
                        default_value = 1;
                end
            case 'XScale'
                default_value = 1.0;
            case 'YScale'
                default_value = 1.0;
            case 'XOffs'
                default_value = 0.0;
            case 'ColorBar'
                default_value = 1;
            case 'ColorMap'
                default_value = [];
            case 'Axes' 
                default_value = 1;
            case 'DisplayTime'
                default_value = 1;
            case 'DisplayFtime'
                default_value = 1;
            case 'DisplayExptime'
                default_value = 1;
            case 'BgrData'
                default_value = [];
            case 'FrameNumber'
                default_value = 0;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
        
        
        
    case 'image_read'
        switch parameter_name
            case 'OrientByExtension'
                switch sys_id
                    case 'mDPC lab'
                        default_value = 0;
                    otherwise
                        default_value = 1;
                end
            case 'DataType'
                default_value = 'double';
            case 'ForceFileType'
                default_value = [];
            case 'MatlabVar' 
                default_value = 'data';
            case 'RowFrom' 
                default_value = 0;
            case 'RowTo' 
                default_value = 0;
            case 'ColumnFrom' 
                default_value = 0;
            case 'ColumnTo' 
                default_value = 0;
            case 'UnhandledParError'
                default_value = 1;
            case 'IsFmask' 
                default_value = true;
            case 'DisplayFilename' 
                default_value = 1;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
           
            
                
    case 'image_orient'
        switch parameter_name
            case 'Transpose'
                default_value = 0;
            case 'FlipLR'
                switch sys_id
                    case 'mDPC lab'
                        default_value = 1;
                    otherwise
                        default_value = 0;
                end
            case 'FlipUD'
                default_value = 0;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
        

    case 'plot_radial_integ' 
        switch parameter_name
            case 'FigNo'
                % figure number for plotting the integrated data
                default_value = 100;
            case 'NewFig'
                % no new figure for each plot
                default_value = 0;
            case 'ClearFig'
                % clear figure before plotting
                default_value = 1;
            case 'Axis'
                % auto scaling
                default_value = [];
            case 'SleepTime'
                % no sleep after each plot
                default_value = 0.0;
            case 'XLog'
                % linear scaling of the x-axis
                default_value = 0;
            case 'YLog'
                % logarithmic scaling of the y-axis
                default_value = 1;
            case 'PlotQ'
                % plot as a function of q rather than pixel number
                default_value = 0;
            case 'PlotAngle'
                % plot as a function of the azimuthal angle rather than q or radius
                default_value = 0;
            case 'RadiusRange'
                % average over this range in radius for the azimuthal plot
                default_value = [];
            case 'FilenameIntegMasks'
                % location of the integration masks, needed for normalization in case of
                % averaging over radii
                default_value = '~/Data10/analysis/data/pilatus_integration_masks.mat';
            case 'PixelSize_mm'
                % pixel size for q calculation
                default_value = [];
            case 'DetDist_mm'
                % detector distance for q calculation
                default_value = [];
            case 'E_keV'
                % x-ray energy for q calculation
                default_value = [];
                % plot in inverse nm rather than inverse Angstroem
            case 'Inverse_nm'
                default_value = 0;
            case 'QMulPow'
                % do not multiply by q to the power of this value
                default_value = [];
            case 'SegAvg'
                % average over segments
                default_value = 1;
            case 'SegRange'
                % segment range
                default_value = [];
            case 'LegendMulSeg'
                % legend in case of multi segment plots
                default_value = 1;
            case 'PointAvg'
                % plot the average over one Matlab file which is typically a scan line
                default_value = 1;
            case 'PointRange'
                % point range
                default_value = [];
            case 'BgrFilename'
                % background to subtract
                default_value = '';
            case 'BgrScale'
                % scaling factor for background data
                default_value = 1.0;
            case 'BgrPoint'
                % point within the background file to subtract
                default_value = 1;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
        
    case 'find_files'
        switch parameter_name
            case 'UseFind'
                % the default is set above system dependent 
            case 'UnhandledParError'
                % exit with an error message if unhandled named parameters are left at the
                % end of this macro
                default_value = 1;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
        
    case 'radial_integ'
        switch parameter_name
            case 'OutdirData'
                % output directory for integrated intensities
                default_value = '~/Data10/analysis/radial_integration/';
            case 'FilenameIntegMasks'
                % location of the integration masks
                default_value = '~/Data10/analysis/data/pilatus_integration_masks.mat';
            case 'rMaxForced'
                % use the full range of integration masks
                default_value = 0;
            case 'FigNo'
                % do not plot integrated data
                default_value = 0;
            case 'SaveCombinedI'
                % combine integrated intensities from all files within one directory
                default_value = 1;
            case 'Recursive'
                % recursively integrate data from all sub directories
                default_value = 1;
            case 'ParTasksMax' 
                % use parallel processing by default if the toolbox is
                % available
                [dummy, other_system_flags] = identify_system();
                default_value = 1;
                if (other_system_flags.parallel_computing_toolbox_available)
                    default_value = 256;
                end
            case 'UseFind'
                % the default value is set above
            case 'UnhandledParError'
                % exit with an error message if unhandled named parameters are left at the
                % end of this macro
                default_value = 1;
            otherwise
                error('Unknown parameter name %s for m-file %s',...
                    parameter_name,mfile_name);
        end
        

    otherwise
        error('No default parameters set for m-file %s (parameter name %s)',...
            mfile_name,parameter_name);
end
