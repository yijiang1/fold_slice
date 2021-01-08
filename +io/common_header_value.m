% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: common_header_value.m,v $
%
% $Revision: 1.7 $  $Date: 2013/01/25 10:22:26 $
% $Author:  $
% $Tag: $
%
% Description:
% return header information which are common to most file formats used at
% the cSAXS beamline
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% - get_hdr_val
%
%
% history:
%
% September 30th 2010:
% add HDF5
%
% October 3rd 2008:
% update FLI date-field since version 1.20 provides a time stamp string
%
% August 28th 2008:
% correct error display for unknown extensions,
% add extension .dat
%
% July 17th 2008: add mar extension
%
% May 7th 2008: 1st version

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

function [value] = common_header_value(header,extension,signature)
import io.image_read
import utils.get_hdr_val

% initialize return value
value = [];

% check number of input arguments
if (nargin ~= 3)
    common_header_value_help();
    error('invalid number of input arguments');
end
    
switch extension
    case 'cbf'
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'Exposure_time',' %f',1);
            case 'Date'
                % The date is available in the Pilatus comments, without
                % any signature. Searching for the 20 string will fail
                % beyond the year 2099
                value = [ '20' get_hdr_val(header,'# 20',' %[^\r]',1) ];
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case {'h5', 'hdf5'}
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'Exposure_time',' %f',1);
            case 'Date'
                value = get_hdr_val(header,'DateTime',' %[^\n]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case 'dat'
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'Exposure_time',' %f',1);
            case 'Date'
                value = get_hdr_val(header,'DateTime',' %[^\n]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end        
    case 'edf'
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'count_time',' = %f',1);
            case 'Date'
                value = get_hdr_val(header,'Date',' = %[^;]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case {'mar','mccd'}
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'ExposureTime_ms',' %f',1)/1000;
            case 'Date'
                value = get_hdr_val(header,'DateTime',' %[^\n]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case {'mat'}
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'Exposure_time',' %f',1);
            case 'Date'
                value = get_hdr_val(header,'DateTime',' %[^\n]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case 'raw'
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'exptimesec',' %f',1);
            case 'Date'
                value = get_hdr_val(header,'version',' %[^\n]',1);
                if (version < 1.20)
                    value = get_hdr_val(header,'FileTimestamp',' %[^\n]',1);
                else
                    value = get_hdr_val(header,'timestamp_string',' %[^\n]',1);
                end
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case 'spe'
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'exposure',' %f',1);
            case 'Date'
                value = [ get_hdr_val(header,'date',' %[^\n]',1)'; ' '; 
                    get_hdr_val(header,'ExperimentTimeLocal',' %[^\n]',1)' ]';
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    case {'tif', 'tiff'}
        switch signature
            case 'ExposureTime'
                value = get_hdr_val(header,'Exposure_time',' %f',1);
            case 'Date'
                value = get_hdr_val(header,'DateTime',' %[^\n]',1);
            otherwise
                error('unknown signature %s for extension %s',...
                    signature,extension);
        end
    otherwise
        common_header_value_help();
        error('unknown extension ''%s''',extension);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = common_header_value_help()

fprintf('Usage:\n');
fprintf('[value]=%s(header,extension,signature);',mfilename);
fprintf('header and extension are returned by image_read\n');
fprintf('The following signatures are recognized:\n');
fprintf('ExposureTime      exposure time in seconds\n');
fprintf('Date              date string in detector specific format\n');
fprintf('Example:\n');
fprintf('exp_time_sec=common_header_value(frame.header{1},frame.extension{1},''ExposureTime'');\n');
