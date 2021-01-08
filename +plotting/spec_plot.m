% function [varargout] = spec_plot(specDatFile, varargin)
% call 'spec_plot()' for further help
%
% please reports bugs, problems, suggestions for improvements to:
% CXS group

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

function [varargout] = spec_plot(specDatFile, varargin)
import io.spec_help
import io.spec_read

% check minimum number of input arguments
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
  unhandled_par_error =  1;
  wait_for            =  0;
  axno                = -1;
  
% parse the variable input arguments
  vararg = cell(0,0);
  for jj = 1:2:length(varargin)
    name = varargin{jj};
    value = varargin{jj+1};
    switch name
      case 'Counter'
        counter = value;
      case 'FigNo' 
        figno = value;        
      case 'Sleep'
        wait_for = value;
      case 'ContMesh'
        saxis = value;
      case 'PFunction'
        pfunction = value;
        counter = 'Pilatus';
        vararg{end+1} = 'OutPut';                               %#ok<AGROW>
        vararg{end+1} = '+pilatus';                             %#ok<AGROW>
      case 'UnhandledParError'
        unhandled_par_error = value;
      otherwise
        vararg{end+1} = name;                                   %#ok<AGROW>
        vararg{end+1} = value;                                  %#ok<AGROW>
    end
  end
  
  vararg{end+1} = 'Cell';
  vararg{end+1} = 1;
  vararg{end+1} = 'UnhandledParError';
  vararg{end+1} = 0;
  [S,vararg_new] = spec_read(specDatFile,vararg);

  vararg = cell(0);
  pilatusDir = cell(0);
  for jj = 1:2:length(vararg_new)
    name = vararg_new{jj};
    value = vararg_new{jj+1};
    switch name
      case 'PilatusDir'
        pilatusDir{end+1} = value;
    otherwise
        vararg{end+1} = name;                                   %#ok<AGROW>
        vararg{end+1} = value;                                  %#ok<AGROW>
    end
  end

  varargout{1} = cell(size(S));
  for jj=1:numel(S)
    S_j = S{jj};
    titlestr = S_j.S;
    if (isfield(S_j,counter))
      scanstr = regexp(S_j.S,' *','split');
      data = S_j.(counter);
      if (iscell(data))
        data = PilatusCounter(data,S_j.PilatusDir,pfunction);
        numall = zeros(size(data));
        for ii=1:numel(data)
          numall(ii) = numel(data{ii});
        end
        if (all(numall==1))
          for ii=1:numel(data)
            numall(ii) = data{ii}{:};
          end
          data = numall;
        end
      end
      if (numel(data)<=1 && ~iscell(data))
        continue;
      end
      if (~exist('figno','var'))
        figno = figure;
      elseif (~ishandle(figno))
        figure(figno)
      else
        switch get(figno,'type')
          case 'figure'
            figure(figno)
          case 'axes'
            axno = figno;
            figno = get(figno,'parent');
          otherwise
            while (strcmp(get(figno,'type'),'axes'))
              figno = get(figno,'parent');
            end
            axno = figno;
            figno = get(figno,'parent');
        end
        if (ishandle(axno))
          axes(axno); 
          cla;
        end
      end
      if (~isempty(regexp(scanstr{3},'[ad][0-3]?scan','match')))
        faxis = scanstr{4};
        fpos = S_j.(faxis);

        h = plot(fpos,data);
        xlabel(faxis); ylabel(counter);
      elseif (~isempty(regexp(scanstr{3},'loopscan','match')))
        faxis = 'Time';
        fpos = S_j.(faxis);

        h = plot(fpos,data);
        xlabel(faxis); ylabel(counter);
      elseif (~isempty(regexp(scanstr{3},'cont_line','match')))
        faxis = scanstr{4};
        if (~exist('saxis','var'))
          fpos = linspace(str2double(scanstr{5}),str2double(scanstr{6}), ...
            str2double(scanstr{7})+1);
          sdata = cell2mat(data{1});
          h = plot(fpos,sdata);
          xlabel(faxis); ylabel(pfunction,'Interpreter','None');
        else
          scannr = str2double(scanstr{2});
          if (~exist('fpos','var'))
            fpos = linspace(str2double(scanstr{5}),str2double(scanstr{6}), ...
                str2double(scanstr{7})+1);
            spos = S_j.(saxis);
            sdata = cell2mat(data{1});
            minnr = scannr;
            maxnr = scannr;
          else
            fpos = [fpos, ...
                linspace(str2double(scanstr{5}),str2double(scanstr{6}), ...
                str2double(scanstr{7})+1)];
            spos = [spos, S_j.(saxis)];
            fpos = unique(fpos);
            spos = unique(spos);
            sdata = [sdata;cell2mat(data{1})];
            minnr = min(minnr,scannr);
            maxnr = max(minnr,scannr);
            titlestr = sprintf('#S %05d-%05d cont_mesh %s %f %f %d %s %f %f %d %s %s', ...
              minnr,maxnr, ...
              faxis, min(fpos), max(fpos), numel(fpos)-1, ...
              saxis, min(spos), max(spos), numel(spos)-1, ...
              scanstr{8},scanstr{9});
          end

          h = imagesc(fpos,spos,sdata);
          if (numel(spos)>1)
            axis equal tight;
            hold on
            contour(fpos,spos,sdata,'k');
            hold off
          end
          xlabel(faxis); ylabel(saxis);
          colorbar;
        end
      elseif (~isempty(regexp(scanstr{3},'mesh','match')))
        faxis = scanstr{4};
        saxis = scanstr{8};
        fpos = unique(S_j.(faxis));
        spos = unique(S_j.(saxis));

        if (numel(data) == numel(fpos)*numel(spos))
          data = transpose(reshape(data,numel(fpos),numel(spos)));
        elseif (mod(numel(data)/(numel(fpos)*numel(spos)),1) == 0)
          data = reshape(data,[],numel(fpos),numel(spos));
          data = transpose(squeeze(sum(data,1)));
        end
        if ((mod(numel(data)/(numel(fpos)*numel(spos)),1) == 0) && ...
            numel(data) > 1)
          h = imagesc(fpos,spos,data);
        else
          [fpos,spos,data,h] = plot_grid(S_j.(faxis),S_j.(saxis),data);
        end
        if (size(data,1)>1)
          axis equal tight;
          hold on
          contour(fpos,spos,data,'k');
          hold off
        end
        xlabel(faxis); ylabel(saxis);
        colorbar;
      elseif (~isempty(regexp(scanstr{3},'round_scan','match')))
        faxis = scanstr{4};
        saxis = scanstr{5};
        [fpos,spos,data,h] = plot_grid(S_j.(faxis),S_j.(saxis),data);
        axis equal tight;
        hold on
        contour(fpos,spos,data,'k');
        hold off
        xlabel(faxis); ylabel(saxis);
        colorbar;
      elseif (~isempty(regexp(scanstr{3},'round_roi_scan','match')))
        faxis = scanstr{4};
        saxis = scanstr{7};
        [fpos,spos,data,h] = plot_grid(S_j.(faxis),S_j.(saxis),data);
        axis equal tight;
        hold on
        contour(fpos,spos,data,'k');
        hold off
        xlabel(faxis); ylabel(saxis);
        colorbar;
      else
        fprintf('unknown scan ''%s'' encountered\n', ...
          scanstr{3})
      end
      title(titlestr,'Interpreter','None')
      drawnow;
      if (wait_for && jj<numel(S))
        pause(wait_for);
      end
    else
      error('Counter ''%s'' not found in scan ''%s''', ...
        counter,S_j.S);
    end
    if (exist('h','var'))
      varargout{1}{jj} = h;
    end
  end
  
  if ((unhandled_par_error) && (~isempty(vararg)))
    vararg                                                      %#ok<NOPRT>
    error('Not all named parameters have been handled.');
  end
  varargout{2} = vararg;
  return
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,Z,h] = plot_grid(xaxis,yaxis,data)
  pos = [xaxis,yaxis];
  diff = pos(2:end,:) - pos(1:end-1,:);
  [y,x] = hist(sqrt(sum(diff.^2,2)),200);
  [~,n] = max(y(2:end));
  s = x(n)/2;
  [X,Y] = meshgrid(linspace(min(xaxis)-s,max(xaxis)+s,200), ...
      linspace(min(yaxis)-s,max(yaxis)+s,200));
  Z = zeros(size(X));
  N = Z;
  for ii=1:numel(data)
      Z = Z+data(ii)*exp(-((X-xaxis(ii)).^2 + (Y-yaxis(ii)).^2)/(2*s^2)); %/(2*pi*s^2);
      N = N+         exp(-((X-xaxis(ii)).^2 + (Y-yaxis(ii)).^2)/(2*s^2)); %/(2*pi*s^2);
  end
  Z = Z./N;
  x = X(1,:);
  y = Y(:,1);
  h = imagesc(x,y,Z);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = PilatusCounter(fnames,fpath,evalfunction)
  if (iscell(fnames{1}))
    data = cell(size(fnames));
  else
    data = -ones(size(fnames));
  end
  for ii=1:numel(fnames)
    if (iscell(fnames{ii}))
      data{ii} = cell(size(fnames{ii}));
      for jj=1:numel(fnames{ii})
        data{ii}{jj} = eval(sprintf(evalfunction, ...
          strcat(fpath,fnames{ii}{jj})));
      end
    else
      data(ii) = eval(sprintf(strcat(evalfunction,';'), ...
          strcat(fpath,fnames{ii}{jj})));
    end
  end
end