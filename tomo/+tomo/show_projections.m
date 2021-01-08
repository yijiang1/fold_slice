% SHOW_PROJECTIONS show playable animation of projections 
% show_projections(img_stack, theta,  par, varargin)      
% 
% Inputs:
%     **img_stack - 3D array of projections 
%     **theta     - angles corresponding to the frames 
%     **par       - parameter structure 
% *optional*
%     **windowautopos     % auomtatic window positioning 
%     **baraxis           % plotted range 
%     **rectangle_pos     % draw rectangle at given coordinates [xmin, xmax, ymin, ymax]
%     **title             % constant title prefix for all projections 
%     **title_extra       % variable title suffix, one cell containing a string per projections !!!
%     **fps               % maximal frame rate 
%     **figure_id         % id number of the created figure 
%     **showsorted        % show projection sorted by angle 
%     **fnct              % data processing function 
%     **init_frame        % starting frame number 
%     **show_grid         % plot grid ovelaying the plotted image 
%     **plot_residua      % plot residua in the image 
%     **plot_only_180_range % plot projection with angle > 180 mirrored and merged with projections from 0-180deg
    
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



function show_projections(img_stack, theta, param, varargin)

    import math.*
    import plotting.*
    
    if nargin < 3
        param = struct();
    end

    parser = inputParser;
    parser.addParameter('windowautopos',  true , @islogical )  % auomtatic window positioning 
    parser.addParameter('baraxis',  'auto'  )   % plotted range 
    parser.addParameter('rectangle_pos',  [], @isnumeric  )   % draw rectangle with at coordinates 
    parser.addParameter('title',  '', @isstr  )   % constant title prefix for all projections 
    parser.addParameter('title_extra',  {}, @iscell  )   % variable title suffix, one cell containing a string per projections
    parser.addParameter('fps', 25, @isnumeric)  % maximal frame rate 
    parser.addParameter('figure_id', 1, @isnumeric)  % maximal frame rate 
    parser.addParameter('showsorted', true, @islogical)  % maximal frame rate 
    parser.addParameter('fnct', @(x)x)   % data processing function 
    parser.addParameter('init_frame', 1, @isnumeric)  % starting frame number 
    parser.addParameter('show_grid', true, @islogical)  % starting frame number 
    parser.addParameter('plot_residua', false, @islogical) % plot residua in the image 
    parser.addParameter('plot_only_180_range', false, @islogical) % plot residua in the image 

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(param, name) || ~ismember(name, parser.UsingDefaults)
            % prefer values in varargin structure 
            param.(name{1}) = r.(name{1});
        end
    end
    
    screensize = get( groot, 'Screensize' );

    nframes = size(img_stack,3); 
    frames = 1:nframes;
    
    fig = plotting.smart_figure(param.figure_id);
    clf()
    
    if param.windowautopos
        win_size = [1060 767]; 
        set(fig,'Outerposition',[150 min(270,screensize(4)-win_size(2)) 1060 767]);
    end
    
    
    if strcmpi(param.baraxis,'auto') && nframes > 1 && isreal(img_stack) && strcmp( func2str( r.fnct ),  '(x)x')
        % set the same range for all the frames using quantile range
        param.baraxis = gather(sp_quantile(img_stack, [1e-3, 1-1e-3], ceil(sqrt(numel(img_stack)))));
    end
    
   if r.plot_only_180_range
       img_stack = fliplr(img_stack(:,:,theta >= 180 | theta < 0));
       theta = mod(theta, 180); 
   end
    

    if param.showsorted && nframes > 1
        [~, ind] = sort(theta);
        frames = ind(frames);

    end
        
    if isa(img_stack,'uint16')  % half precision data 
       param.fnct = @(x)(param.fnct(fp16.get(x)));  % covnvert to singles first  
    end
    
    for ii = 1:nframes
        num=ii; % frames(ii);
        try
            title_list{ii} = sprintf('%s Scan: %05d  Projection: %03d Angle=%5.2f deg Energy=%.3f kev', param.title, param.scanstomo(num), num, theta(num), param.energy(num)); 
        catch
            title_list{ii} = sprintf('%s Projection: %03d', param.title, num); 
        end
        if ~isempty(param.title_extra) && length(param.title_extra) == nframes
            %title_list{ii} = [title_list{ii}, ' ',  param.title_extra{num}]; 
            title_list{ii} = {title_list{ii};  param.title_extra{num}}; 

        end
    end
    
    slider_default = [.15 0.01 0.7 0.05];
    play_default = [slider_default(1)-0.1 slider_default(2) 0.08 0.05];
    edit_default = [slider_default(1)+slider_default(3)+0.01 slider_default(2) 0.08 0.05];


    imagesc3D(img_stack, 'title_list', title_list, 'fps', param.fps , 'slider_position',slider_default , ...
         'play_position', play_default , 'edit_position', edit_default,'fnct', param.fnct, 'order', frames,...
         'init_frame', param.init_frame, 'loop', true, 'plot_residua', param.plot_residua )
    colormap bone(256); 
    axis xy equal tight;
    colorbar

    if ~isempty(param.rectangle_pos)
        hold all 
        rectangle('Position',[param.rectangle_pos(1), param.rectangle_pos(3), param.rectangle_pos(2)-param.rectangle_pos(1), param.rectangle_pos(4)-param.rectangle_pos(3)],'edgecolor','r') 
        hold off 
    end
  

    if ~strcmpi(param.baraxis,'auto')
        caxis(param.baraxis)
    end
    if param.show_grid
       grid on  
    end
    drawnow 


    
 
end
