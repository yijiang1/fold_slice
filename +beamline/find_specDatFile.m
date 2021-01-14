% function specDatFile = find_specDatFile(specDatFile)
% find location of the spec file in the provided folder / path 
% if the variable specDatFile is not a complete path to a file
%   try to guess where a spec data file can be found, by
%   - look for directories called 'spec' or 'dat-files'
%   - look for files called '*.dat'
%   - take the newest one


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

function specDatFile = find_specDatFile(specDatFile)
  while (exist(specDatFile,'file') ~= 2)
% if the variable specDatFile is not a complete path to a file
%   try to guess where a spec data file can be found, by
%   - look for directories called 'spec' or 'dat-files'
%   - look for files called '*.dat'
%   - take the newest one
    compare_str = specDatFile;
    fname = dir(specDatFile);
    if (exist(specDatFile,'dir'))
      if (specDatFile(end) ~= '/')
        specDatFile = strcat(specDatFile,'/');
      end
      
      for ii=1:numel(fname)
        if (regexp(fname(ii).name,'.dat$'))
          specDatFile = strcat(specDatFile,'*.dat');
          fname = [];
          break;
        end
      end
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,'dat-files'))
          specDatFile = strcat(specDatFile,fname(ii).name);
          fname = [];
          break;
        end
      end
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,'specES1'))
          specDatFile = strcat(specDatFile,fname(ii).name);
          break;
        end
        if (strcmp(fname(ii).name,'spec'))
          specDatFile = strcat(specDatFile,fname(ii).name);
          break;
        end
      end
    else
      if (numel(fname)>0)
        [~,ii] = max(cell2mat({fname.datenum}));
        specDatFile = regexprep(specDatFile,'\*\.dat$',fname(ii).name);
      else
        error('''%s'' cannot be found.', specDatFile);
        break
      end
    end
    if (strcmp(specDatFile,compare_str))
      break
    end
  end
  
  if ~exist(specDatFile, 'file') || exist(specDatFile, 'dir')
     error('Spec dat file not found in the provided path %s', specDatFile) 
  end
  
end
