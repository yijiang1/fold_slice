% function [varargout] = spec_read(specDatFile, varargin)
% call 'spec_read()' for further help
%
% please reports bugs, problems, suggestions for improvements to:
% CXS group
%
% the case of duplicate scannumbers in the spec file, and how to
%   address them remains to be implemented

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

function [varargout] = spec_read(specDatFile, varargin)
import io.*
import utils.find_files

% check minimum number of input arguments
  if(~isunix)
    error(['The function ''spec_read'' relies on various system calls,\n' ...
      '  currently available only on UNIX-like systems\n']);
  end
  if (nargin<1) 
    spec_help(mfilename);
    error('At least the spec data filename has to be specified as input parameter.')
  end
  
% accept cell array with name/value pairs as well
  no_of_in_arg = nargin;
  if (nargin == 2)
    if (isempty(varargin))
      % ignore empty cell array
      no_of_in_arg = no_of_in_arg -1;
    else
      if (iscell(varargin{1}))
        % use a filled one given as first and only variable parameter
        varargin = varargin{1};
        no_of_in_arg = 1 + length(varargin);
      end
    end
  end
% check number of input arguments
  if (rem(no_of_in_arg,2) ~= 1)
      error('The optional parameters have to be specified as ''name'',''value'' pairs');
  end
  
% here the default parameters  
  scannr     = -1;
  output_arg = {'meta','motor','counter'};
  unhandled_par_error = 1;
  force_cell = 0;
% other parameters
  burstn0 = 1;
  mexpn0  = 1;
  creader = 1;
  orchestraPath = [];
  specReadOptions = [];
  
% parse the variable input arguments
  vararg = cell(0,0);
  for jj = 1:2:length(varargin)
    name = varargin{jj};
    value = varargin{jj+1};
    switch name
      case 'ScanNr'
        if (~all(value==round(value)) || any(value==0))
          error('The value of ''%s'' needs to be a (vector of) integer(s) not equal to zero', name)
        end
        scannr = value;
      case 'Burst'
        if (~all(value==round(value)) || any(value<=0))
          error('The value of ''%s'' needs to be a positive integer', name)
        end
        burstn0 = max(1,value);
        mexpn0 = 1;
      case 'MultExp'
        if (~all(value==round(value)) || any(value<=0))
          error('The value of ''%s'' needs to be a positive integer', name)
        end
        burstn0 = 1;
        mexpn0 = max(1,value);
      case 'OutPut'
        output_arg = set_TextFlag(output_arg,value);
      case 'PilatusMask'
        pilatusMask = value;
        output_arg = set_TextFlag(output_arg,'+pilatus');
      case 'PilatusDir'
        pilatusDir0 = value;
        output_arg = set_TextFlag(output_arg,'+pilatus');
      case 'Cell'
        force_cell = value;
      case 'CReader'
        creader = value;
      case 'Orchestra'
        orchestraPath = value;
      case 'UnhandledParError'
        unhandled_par_error = value;
      otherwise
        vararg{end+1} = name;                                   %#ok<AGROW>
        vararg{end+1} = value;                                  %#ok<AGROW>
    end
  end
    
  specDatFile = beamline.find_specDatFile(specDatFile);
  if creader
      try
          if ~isempty(orchestraPath)
              specReadOptions = [specReadOptions ' --orchestra ' orchestraPath];
          end
          sysCall = [ fullfile(fileparts(mfilename('fullpath')), './spec_reader/spec_reader') ' -s ' specDatFile ' --scanNr ' strjoin(arrayfun(@(x) num2str(x),scannr,'UniformOutput',false),',') specReadOptions];
          [stat, out] = system(sysCall);
          if stat
              error(['Fast spec reader failed: ' strtrim(out)])
          end
          res = jsondecode([ out ]);
          if length(res)==1
              varargout{1} = res;
          else
              varargout{1} = num2cell(res)';
          end
      catch ME
          if creader
              warning('Fast spec reader failed.')
              disp(ME)
              creader = 0;
          end
      end
  end
  if ~creader
      % reading the spec file for the scans
      if ~exist(specDatFile, 'file')
          error('Spec file was not found in %s', specDatFile)
      end
      [scans,scanline] = spec_readscan(specDatFile, scannr);
      
      % reading the spec file for configuration information
      [motorc,motorline] = spec_readconfig(specDatFile);
      motorn = cell(size(scannr));
      for jj=1:numel(scans)
          motorn{jj} = motorc{find(motorline<scanline(1),1,'last')};
      end
      
      
      varargout{1} = cell(size(scans));
      % copying vararg in order to keep its content for each iteration
      vararg_1 = vararg;
      vararg_1{end+1} = 'UnhandledParError';
      vararg_1{end+1} = 0;
      for jj=1:numel(scans)
          scanstr = scans{jj}';
          arrout = regexp(scanstr{1},' +','split');
          scannr(jj) = str2double(arrout{2});
          
          % get the motors
          line = ~cellfun(@isempty,regexp(scanstr,'^#P'));
          arrout = reshape(transpose(strvcat((scanstr(line)))),1,[]);
          arrout = regexp(arrout,' +|(#P[0-9]+)','split');
          motorv = str2double(arrout(~cellfun(@isempty,arrout)));
          
          % get the counters
          line = ~cellfun(@isempty,regexp(scanstr,'^#L'));
          counters = regexp(scanstr{line},' +','split');
          counters = counters(2:end);
          
          % get the data
          line  = ~cellfun(@isempty,regexp(scanstr,'^[0-9+-]'));
          data = cellfun(@str2num,scanstr(line),'UniformOutput',false);
          data = cell2mat(data(:));
          
          
          % get Pilatus information if required
          if (exist('pilatusMask','var'))
              if (exist('pilatusDir0','var'))
                  pilatusDir = find_pilatusDir(pilatusDir0,scannr(jj));
              else
                  pilatusDir = '';
              end
              
              [pilatusDir, pilatusName, vararg] = ...
                  find_files(strcat(pilatusDir,pilatusMask), vararg_1);
              
              pilatusInd = zeros(numel(pilatusName),2);
              for ii=1:numel(pilatusName)
                  arrout = regexp(pilatusName(ii).name,'[_\.]','split');
                  pilatusInd(ii,:) = str2double(arrout(end-2:end-1));
              end
          end
          
          % make a structure
          scan_structure = struct;          
          if (any(strcmp(output_arg,'motor')) || ...
                  any(strcmp(output_arg,'motors')))
              for ii=1:numel(motorn{jj})
                  scan_structure.(motorn{jj}{ii}) = motorv(ii);
              end
          end
          if (any(strcmp(output_arg,'meta')))
              scan_structure.S = scanstr{1};
              
              % get the time stamp and convert it to a MATLAB-compatible date string
              line = ~cellfun(@isempty,regexp(scanstr,'^#D'));
              datestr = regexp(scanstr{line},' +','split');
              scan_structure.D = sprintf('%s-%s-%s %s', ...
                  datestr{4}, ...
                  datestr{3}, ...
                  datestr{6}, ...
                  datestr{5});
              
              line = ~cellfun(@isempty,regexp(scanstr,'^#C meta'));
              meta = scanstr(line);
              for ii=1:length(meta)
                  metaEntry = strsplit(strtrim(meta{ii}), ' ');
                  switch metaEntry{4}
                      case {'int', 'float'}
                          converterFunc = @str2double;
                      case 'string'
                          converterFunc = @(x) x;
                  end
                  
                  tmp = [];
                  for jjj=1:str2double(metaEntry{5})
                      tmp = [tmp converterFunc(metaEntry{5+jjj})];
                  end
                  scan_structure.(metaEntry{3}) = tmp;
              end
                          
              
          end        
          if (any(strcmp(output_arg,'counter')) || ...
                  any(strcmp(output_arg,'counters')))
              if (numel(data)>0)
                  for ii=1:numel(counters)
                      scan_structure.(counters{ii}) = data(:,ii);
                  end
              end
          end
          
          if (any(strcmp(output_arg,'pilatus')))
              if (burstn0>1)
                  burstn = burstn0;
              elseif (isfield(scan_structure,'burstn'))
                  burstn = max(ones(size(scan_structure.burstn)),scan_structure.burstn);
              else
                  burstn = burstn0;
              end
              if (mexpn0>1)
                  mexpn = mexpn0;
              elseif (isfield(scan_structure,'mexpn'))
                  mexpn = max(ones(size(scan_structure.burstn)),scan_structure.mexpn);
              else
                  mexpn = mexpn0;
              end
              
              if (numel(burstn)==1 && numel(mexpn)==1)
                  % these pseudo motors shouldn't be scanned and are
                  %   therefore presumably constant
                  Pilatus = cell(size(data,1),1);
                  scanInd = [transpose(floor((0:size(data,1)-1)/mexpn)), ...
                      transpose(mod(0:size(data,1)-1,mexpn))];
                  for ii=1:size(data,1)
                      if (burstn<=1)
                          indices = all(ones(size(pilatusName))'*scanInd(ii,:)==pilatusInd,2);
                      else
                          indices = all(ones(size(pilatusName))'*scanInd(ii,1)==pilatusInd(:,1),2);
                      end
                      Pilatus{ii} = {pilatusName(indices).name};
                      if (numel(indices)<burstn)
                          % accounting for missing frames
                          % untested, as of Dec 16
                          tmp = Pilatus{ii};
                          Pilatus{ii} = cell(1,burstn);
                          for kk=1:numel(tmp)
                              Pilatus{ii}{pilatusInd(2,kk)} = tmp{kk};
                          end
                      end
                  end
              end
              scan_structure.Pilatus = Pilatus;
              scan_structure.PilatusDir = pilatusDir;
          end
          
          % append orchestra data
          if ~isempty(orchestraPath)
              [config, data] = readOrchestraFile(orchestraPath, scannr(jj));
              for oidx=1:length(config)
                  scan_structure.(config{oidx}) = data{oidx};
              end
          end
          varargout{1}{jj} = scan_structure;
      end
      if ((unhandled_par_error) && (~isempty(vararg)))
          vararg                                                      %#ok<NOPRT>
          error('Not all named parameters have been handled.');
      end
      if (numel(varargout{1})==1 && ~force_cell)
          varargout{1} = varargout{1}{1};
      end
      varargout{2} = vararg;
      return
  end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = spec_readscan(specDatFile, scannr)
import io.*
import utils.*
% check minimum number of input arguments
  if (nargin<2) 
    spec_read_help(mfilename);
    error('usage ''spec_readscan(specDatFile, scannr)''');
  end
  
  cmd = sprintf('grep -n ''#S '' %s', specDatFile);
  [~,sysout] = system(cmd);
  arrout = regexp(sysout,'[:\n]','split');
  arrout = arrout(~cellfun(@isempty,arrout));
  arrout = transpose(reshape(arrout,2,[]));
  [s1,~] = size(arrout);
  
  if (any(scannr>0))
    tmp = regexp(arrout(:,2),' +','split');
    scanID = zeros(size(tmp));
    for ii=1:numel(tmp)
      scanID(ii) = str2double(tmp{ii}{2});
    end
  end
  
% creating a look up table where the scans begin and end
  lineI = zeros(size(scannr));
  scans = cell(size(scannr));
  for ii=1:numel(scannr)
    if (scannr(ii)<0)
      if (-scannr(ii)>size(arrout,1))
        error('There does not seem to be sufficiently many scans in ''%s''', ...
          specDatFile);
      end
      tmp = regexp(arrout(end+scannr(ii)+1,:),' +','split');
      lineI(ii) = str2double(tmp{1});
      if (scannr(ii) < -1)
        tmp = regexp(arrout(end+scannr(ii)+2,:),' +','split');
        lineF = str2double(tmp{1})-1;
      else
        cmd = sprintf('wc -l %s', specDatFile);
        [~,sysout] = system(cmd);
        lineF = str2double(regexp(sysout,'[0-9]+ ','match'));
      end
      scannr(ii) = str2double(tmp{2}(2));
    else
      tmp = find((scannr(ii)==scanID),1,'last');
      try 
          verbose(2,'reading ''%s''',arrout{tmp,2});
      catch
          fprintf('reading ''%s''\n',arrout{tmp,2});
      end

      if (isempty(tmp))
        fprintf('Do not find scan #S %d in ''%s''.\n', ...
          scannr(ii),specDatFile);
        continue
      end
      lineI(ii) = str2double(arrout(tmp,1));
      if (tmp<s1)
        lineF = str2double(arrout(tmp+1,1))-1;
      else
        cmd = sprintf('wc -l %s', specDatFile);
        [~,sysout] = system(cmd);
        lineF = str2double(regexp(sysout,'[0-9]+ ','match'));
      end
    end

    cmd = sprintf('head -n +%d %s | tail -n -%d', ...
      lineF, specDatFile, lineF+1-lineI(ii));

  % read the scan
    [~,sysout] = system(cmd);
    scans{ii} = regexp(sysout,'\n','split');
    scans{ii} = scans{ii}(~cellfun(@isempty,scans{ii}));
  end
  varargout{1} = scans;
  varargout{2} = lineI;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [motors,linenr] = spec_readconfig(specDatFile)
import io.*
import utils.find_files
  % finding the end of the last config
  cmd = sprintf('grep -nE ''^#O[0-9]+ '' %s', ...
      specDatFile);
  [~,sysout] = system(cmd);
  arrout = regexp(sysout,'[:\n]','split');
  arrout = arrout(~cellfun(@isempty,arrout));
  arrout = transpose(reshape(arrout,2,[]));
  [s1,~] = size(arrout);
  
  % devide this collection of configuration information into single chuncks
  linenr = '';
  motors = cell(0);
  for ii=1:s1
    if (~isempty(regexp(arrout{ii,2},'^#O0+','match')))
      if (isempty(linenr))
        linenr = str2double(arrout{ii,1});
      else
        linenr = vertcat(linenr, ...
          str2double(arrout{ii,1}));
      end
      motors{numel(linenr)} = regexp(arrout{ii,2},'(^#O[0-9]+ *)| *','split');
    else
      motors{numel(linenr)} = horzcat(motors{numel(linenr)}, ...
        regexp(arrout{ii,2},'(^#O[0-9]+ *)| *','split'));
    end
  end
  for ii=1:numel(linenr)
    motors{ii} = motors{ii}(~cellfun(@isempty,motors{ii}));
  end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pilatusDir = find_pilatusDir(pilatusDir,scannr)
import io.*
import utils.find_files
  while (any(exist(pilatusDir,'dir') == [2,7]))
% if the variable pilatusDir is not a complete path to a file
%   try to guess where a spec data file can be found, by
%   - look for directories called 'spec' or 'dat-files'
%   - look for files called '*.dat'
%   - take the newest one
    compare_str = pilatusDir;
    fname = dir(pilatusDir);
    if (exist(pilatusDir,'dir'))
      if (pilatusDir(end) ~= '/')
        pilatusDir = strcat(pilatusDir,'/');
      end
      
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,sprintf('S%02d000-%02d999',floor(scannr/1e3),floor(scannr/1e3))) || ...
            strcmp(fname(ii).name,sprintf('S%05d',scannr)))
          pilatusDir = strcat(pilatusDir,fname(ii).name);
          fname = [];
          break;
        end
      end
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,'pilatus'))
          pilatusDir = strcat(pilatusDir,fname(ii).name);
          break;
        end
      end
    end
    if (strcmp(pilatusDir,compare_str))
      break
    end
  end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flags] = set_TextFlag(flags,value)
import io.*
import utils.find_files
  if (~iscell(value))
    value = {value};
  end
  for ii=1:numel(value)
    value_i = value{ii};
    if (regexp(value_i,'^[a-zA-Z]+$'))
      flags = {value_i};
      if (ii+1<=numel(value))
        value{ii+1} = regexprep(value{ii+1},'^([a-zA-Z]*)','+$1');
      end
    elseif (regexp(value_i,'^\+[a-zA-Z]+$'))
      flags = unique([flags,value_i(2:end)]);
      if (ii+1<=numel(value))
        value{ii+1} = regexprep(value{ii+1},'^([a-zA-Z]*)','+$1');
      end
    elseif (regexp(value_i,'^-[a-zA-Z]+$'))
      flags = flags(~strcmp(flags,value_i(2:end)));
      if (ii+1<=numel(value))
        value{ii+1} = regexprep(value{ii+1},'^([a-zA-Z]*)','-$1');
      end
    end
  end
end

function [config, data] = readOrchestraFile(orchestraPath, scanNr)
config = [];
data = [];
fname = fullfile(orchestraPath, sprintf('scan_%05d.dat', scanNr));
if ~exist(fname, 'file')
    warning('Could not find any orchestra file for scan %d', scanNr);
    return;
end
fid = fopen(fname);
fgetl(fid);
config = strsplit(fgetl(fid), ' ');

data=textscan(fid,repmat('%f ', 1, length(config)));

end