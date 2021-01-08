% [] = spec_help(m_file_name)

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

function [] = spec_help(m_file_name)
import io.*

  switch m_file_name
    case 'spec_read'
      fprintf(['usage: spec_read(specDatFile)\n' ...
        '       spec_read(specDatFile,varargin)\n' ...
        'where ''varagin'' is in ''name''-''value''pairs\n' ...
        'the output is a structure or a cell of structures\n\n']);
      name = {'ScanNr', ...
        'Burst', 'MultExposure', ...
        'OutPut','Cell', ...
        'PilatusMask','PilatusPath', ...
        'ValidMask', ...
        'UnhandledParError'};
    case 'spec_plot'
      fprintf(['usage: spec_plot(specDatFile)\n' ...
        '       spec_plot(specDatFile,varargin)\n' ...
        'where ''varagin'' is in ''name''-''value''pairs\n' ...
        'the output is a cell containing graphics handles\n\n']);
      name = {'ScanNr', ...
        'Burst', 'MultExposure', ...
        'Counter','FigNo','Sleep', ...
        'PilatusMask','PilatusPath', ...
        'ValidMask', ...
        'UnhandledParError'};
    otherwise
      fprintf('unknown ''M-file''\n')
  end
  outputstr = cell(1,numel(name)+1);
  outputstr{1} = {'name','value','note',''};
  filler = 4*ones(1,2);
  for ii=1:numel(name)
    note2 = cell(0);
    switch name{ii}
      case 'Burst'
        value = '<pos. integer>';
        note  = 'number of bursts per point';
      case 'MultExp'
        value = '<pos. integer>';
        note  = 'number of spec-controlled exposures per point';
      case 'Cell'
        value = '[01]';
        note  = 'to force the output to be a cell even in case of a single scan.';
      case 'Counter'
        value = '<spec counter>';
        note  = 'counter whose scalar output is to be displayed';
      case 'FigNo'
        value = '<handle>';
        note  = 'handle to determine graphical output';
      case 'PilatusMask'
        value = '<mask>';
        note  = 'mask (including wild cards) how a Pilatus data file will be named';
      case 'PilatusPath'
        value = '<path>';
        note  = 'path in which to expect the corresponding Pilatus data files';
      case 'OutPut'
        value = '[''meta''|''counter''|''motor'']';
        note  = 'what part of the spec data file is being output';
        note2 = {'These parameters can be combined in a cell.', ...
          'If preceeded by ''-'', the corresponding output is suppressed.' ...
          'If preceeded by ''+'', the corresponding output is enabled.' ...
          'Without prefix, only the output corresponding to ''value'' is enabled.'};
      case 'ScanNr'
        value = '<integer>';
        note  = 'indicating the scan number to be read\n';          
        note2 = {'In case of duplicate scan numbers, the last one is used only.', ...
          'In case of negative input, scans are counted backwords from the end of the  data file.'};
      case 'Sleep'
        value = '<in seconds>';
        note = 'delay in seconds between updating subsequent plots';
      case 'UnhandledParError'
        value = '[01]';
        note  = 'if an error occurs if not all named parameters have been handled.';
      case 'ValidMask'
        value = '<path>';
        note  = 'valid pixel mask';
      otherwise
        fprintf('unknown argument ''%s''\n', ...
          name{ii});
    end
    outputstr{ii+1} = {name{ii}, value, note, note2};
    if (numel(name{ii}) > filler(1))
      filler(1) = numel(name{ii});
    end
    if (numel(value) > filler(2))
      filler(2) = numel(value);
    end
  end
  outputformat = {sprintf('%%-%ds %%-%ds -> %%s\\n', filler), ...
    sprintf('%%-%ds %%s\\n', sum(filler)+4)};
  for ii=1:numel(outputstr)
    fprintf(outputformat{1},outputstr{ii}{1:3})
    if (~isempty(outputstr{ii}{4}))
      for jj=1:numel(outputstr{ii}{4})
        fprintf(outputformat{2},'',outputstr{ii}{4}{jj})        %#ok<CTPCT>
      end
    end
  end
  fprintf('\n');
end