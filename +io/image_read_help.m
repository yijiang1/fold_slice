% image_read_help(extension,m_file_name,varargin)
% parameter help for image_read

% Filename: $RCSfile: image_read_help.m,v $
%
% $Revision: 1.3 $  $Date: 2013/01/25 10:23:47 $
% $Author:  $
% $Tag: $
%
% Description:
% parameter help for image_read
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
% July 17th 2008: 
% add ForceFileType parameter and support for MAR CCD TIFF
%
% May 9th 2008: 1st version

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

function [] = image_read_help(extension,m_file_name,varargin)
import io.*

% check minimum number of input arguments
if (nargin < 2)
    error('At least the extension and m-file name have to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 3)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 2 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

% parse the variable input arguments
examples = 1;
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'Examples' 
            examples = value;
        otherwise
            % pass unknown parameters to image_read_sub_help
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end

% do not display examples from image_read_sub_help
vararg{end+1} = 'Examples';
vararg{end+1} = 0;

image_read_sub_help(m_file_name,extension,vararg)
fprintf('''DataType'',<Matlab class>            default is ''double'', other possibilities are ''single'', ''uint16'', ''int16'', ''uint32'', etc.\n');
fprintf('                                     The conversion is done using ''cast'', i.e, out-of-range values are mapped to the minimum or maximum value\n');
fprintf('''ForceFileType'',<''extension''>        force the file types to be recognized by the here specified extension,\n');
fprintf('                                     useful in case of no or other types of extensions, used by default as OrientExtension as well.\n');
fprintf('                                     The extension ''mar'' and ''mccd'' can be used to read MAR CCD TIFF data.\n');
fprintf('''RowFrom'',<0-max>                    region of interest definition, 0 or 1 for full frame\n');
fprintf('''RowTo'',<0-max>                      region of interest definition, 0 for full frame\n');
fprintf('''ColumnFrom'',<0-max>                 region of interest definition, 0 or 1 for full frame\n');
fprintf('''ColumnTo'',<0-max>                   region of interest definition, 0 for full frame\n');
image_orient_help(m_file_name,'ParametersOnly',1);
fprintf('''IsFmask'',<0-no,1-yes>               interprete the filename(s) as search mask that may include wildcards, default true\n');
fprintf('''DisplayFilename'',<0-no,1-yes>       display filename of a file before loading it, default yes\n');
fprintf('''UnhandledParError'',<0-no,1-yes>     exit in case not all named parameters are used/known, default is yes\n');
fprintf('\n');
fprintf('HDF5, H5 or NeXus specifics           These files contain data and metadata hierarchically organized in groups and datasets,\n');
fprintf('                                      each group or dataset can also have attributes. Such files are thus here treated in a special way.\n');
fprintf('                                      If you provide only filename then the file contents, including links but excluding attributes,\n');
fprintf('                                      will be recursively read and returned as a Matlab structure. See also hdf5_load.m\n');
fprintf('''H5Location'',<location>               If <location> is a group then it will be read recursively and returned as a Matlab structure.\n');
fprintf('                                      If <location> is a dataset, the dataset will be read and returned within the field ''data'',\n');
fprintf('                                      this is done in an effort to be compatible with the output of image_read for other file extensions. \n');
fprintf('                                      Only in this case the data region options will be used, e.g ''RowFrom'', ''RowTo'', etc. \n');
fprintf('''FrameRange'',<[first_fr last_fr]>     Read only a subset of the frames available in the HDF5 file dataset specifed with ''H5Location''\n');
fprintf('                                      This will only have an effect if ''H5Location'' is a dataset and not a group \n');
fprintf('''ReadAttr'',<0-no,1-yes>               Read the attributes of a dataset or group (default 0). The Name and Value of the attributes are \n');
fprintf('                                      returned in a structure. Note with this option only the attributes (and not the dataset) are read\n');

if (examples)
    fprintf('\n');
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[frame]=%s(''~/Data10/pilatus/image_1_ct.cbf'');\n',...
        m_file_name);
    fprintf('[frame]=%s({''~/Data10/pilatus/image_1_ct1.cbf'',''~/Data10/pilatus/image_1_ct2.cbf''});\n',...
        m_file_name);
    fprintf('[frame]=%s(''~/Data10/pilatus/S00010/*.cbf'',''IsFmask'',1);\n',...
        m_file_name);
    fprintf('[frame]=%s(''~/Data10/pilatus/image_1_ct.cbf'',''RowFrom'',500,''RowTo'',600);\n',...
        m_file_name);
    fprintf('\n');
    fprintf('The returned structure has the fields data, header, filename and extension.\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('Examples for HDF5:\n');
    fprintf('[data] = image_read(''scan_00300.hdf5'')                                               Read all data in the file.\n')
    fprintf('[data] = image_read(''scan_00300.hdf5'',''H5Location'',''/entry/instrument'')              Reads NeXus instrument group.\n')
    fprintf('[data] = image_read(''scan_00300.hdf5'',''H5Location'',''/entry/collection/data/spec'')    Reads spec data which includes counters and motors that change during a scan.\n')   
    fprintf('[data] = image_read(''scan_00300.hdf5'',''H5Location'',''/entry/collection/data/spec'',''ReadAttr'',1)      Reads spec data that did not change during the scan, e.g. static motors.\n')   
    fprintf('[data] = image_read(''scan_00300.hdf5'',''H5Location'',''/entry/instrument/Pilatus_2M/data'')             Reads all Pilatus frames from the scan.\n')   
    fprintf('[data] = image_read(''scan_00300.hdf5'',''H5Location'',''/entry/instrument/Pilatus_2M/data'',''FrameRange'',[5 10], ''RowFrom'',500,''RowTo'',Inf,''ColumnFrom'',200,''ColumnTo'',800 )\n')   
    fprintf('                                      Reads the specified frame range and region of interest of the pilatus frames.\n')   
end
