%IMAGESC3D 3D wrapper for imagesc
% imagesc3D supports the same parameters as Matlab's imagesc. In addition, the following 
% parameters can be set
%
% init_frame...       starting frame number (default 1)
% slider_axis...      axis along which you want to use imagesc (default 3)
% fps...              frames per second (default 25); will be adjusted by a factor of 1.2 to account for internal overhead
% title_list...       individual title for each frame (default {})
% loop...             run in a loop (default false)
% reset_frame...      stop resets frame to init_frame (default false)
% autoplay...         stark movie automatically (default false)
% slider_position...  slider position [left bottom width height] (default center of axis)
% play_position...    play button position [left bottom width height]
% edit_position...    edit box position [left bottom width height]
% show_play_button... show/hide button; needs to be visible if loop=true; (default true)
% show_edit_box...    show/hide box
% fnct...             data processing function
% order...            change slice order in stack
% save_movie...       specify filename if a movie shall be written
% movie_quality...    image quality of the saved movie
%
% Complex images will be converted to RGB using c2image.
%
% If you are not using 'autplay', you can also set a global title instead
% of a title list (similar to imagesc) and use '%d' to get the slice number
%    title('Random block - slice %d')
% 
%
% EXAMPLES:
%       imagesc3D(rand(256, 256, 100), 'fps', 10, 'loop', true)
%       imagesc3D(rand(256, 256)*1j)
%       imagesc3D(rand(20, 256, 256), 'slider_axis', 1);
%
%
%    Additionally, you can use imagesc/imagesc3D routines and trigger the movie by
%    calling the play method of a specified axis:
%
%       figure(1); 
%       imagesc3D(rand(256, 256, 100), 'fps', 20);
%       title('Random block - slice %d');
%       colorbar();
%       ax = gca;
%       ax.play();
%


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function imagesc3D(varargin)
import math.isint
import plotting.c2image

ax_img = {};
if nargin == 1
    img = varargin{1};
    vararg = {};
elseif nargin == 3 && isnumeric(varargin{1}) && isnumeric(varargin{2})  && (islogical(varargin{3}) || isnumeric(varargin{3}))
    img = varargin{3};
    ax_img = varargin(1:2);
    vararg = {};
elseif (islogical(varargin{1}) || isnumeric(varargin{1})) && ischar(varargin{2})
    % assume that first argument is images, and next are
    % string+arguments
    img = varargin{1};
    vararg = varargin(2:end);
elseif isnumeric(varargin{1}) && isnumeric(varargin{2})  && (islogical(varargin{3}) || isnumeric(varargin{3})) 
    % assume that first two arguments are axis,and third is images, and next are
    % string+arguments
    ax_img = varargin(1:2);
    img = varargin{3};
    vararg = varargin(4:end);
else
   error('Unknown combination of parameters') 
end

ax = gca;
pos = ax.Position;
slider_default = [pos(1)+pos(3)/2-0.06 pos(2)-0.1 0.14 0.05];
play_default = [slider_default(1)-0.1 pos(2)-0.1 0.08 0.05];
edit_default = [slider_default(1)+slider_default(3)+0.01 slider_default(2) 0.08 0.05];

par = inputParser;
par.addParameter('fps', 25, @isnumeric)  % maximal frame rate 
par.addParameter('init_frame', 1, @isnumeric)  % starting frame number 
par.addParameter('title_list', {}, @iscell)  % list of titles for each frame
par.addParameter('slider_axis',3, @isnumeric)  % array axis
par.addParameter('loop', false, @islogical) % loop
par.addParameter('reset_frame', false, @islogical) % stop resets frame to init_frame
par.addParameter('autoplay', false, @islogical) % start loop automatically
par.addParameter('slider_position',slider_default, @isnumeric)  % slider position; [left bottom width height]
par.addParameter('play_position',play_default, @isnumeric)  % slider position; [left bottom width height]
par.addParameter('edit_position', edit_default, @isnumeric) % edit position; [left bottom width height]
par.addParameter('show_play_button',true, @islogical)  % array axis
par.addParameter('show_edit_box', true, @islogical) % edit box
par.addParameter('fnct', @(x)x)   % data processing function 
par.addParameter('order', 1:size(img,3), @isnumeric) % change slice order in stack
par.addParameter('plot_residua', false, @islogical) % plot residua in the image 
par.addParameter('save_movie', '', @ischar) % specify filename if a movie shall be written
par.addParameter('movie_quality', 80, @isnumeric) % movie quality


par.parse(vararg{:})
vars = par.Results;

vars.fps = vars.fps *1.2; % correct for overhead

% permute the array to slide along diferent axis
switch  vars.slider_axis
    case 1
        img = rot90(permute(img,[2,3,1]));
    case 2
        img = rot90(permute(img,[1,3,2]));
end
if any(cellfun(@(x)(strcmpi(x, 'order')),  par.UsingDefaults))
    % redefine the order just in case that the axis were swapped, but only
    % if there is  not use preference 
    vars.order = 1:size(img,3); 
end


if ~isempty(vars.title_list)
   assert(length(vars.title_list) == size(img,3), 'Number of titles has to correspond to number of frames') 
end

sz = size(img,3);

im = imhandles(gcf);
ax = gca;
if ~isprop(ax, 'index')
    ax.addprop('index');
    ax.index = length(im)+1;
else
    if isprop(ax, 'play_handle')
        delete(ax.play_handle);
    end
    if isprop(ax, 'slider_handle')
        delete(ax.slider_handle);
    end
    if isprop(ax, 'edit_handle')
        delete(ax.edit_handle);
    end
    if isprop(ax, 'vars')
        ax.vars = [];
    end
end

if sz>1
    % checks
    vars.init_frame = round(vars.init_frame); 

    if vars.init_frame > sz || vars.init_frame  < 1
        warning('Initial frame exceeds stack size.')
        vars.init_frame = 1;
    end
        
    vars.vargin = ax_img;
    
    if ~ax.isprop('img')
        ax.addprop('img');
    end

    ax.img = img;

    if ~ax.isprop('play')
        ax.addprop('play');
    end
    ax.play = @(x)play(x);

    if ~ax.isprop('stop')
        ax.addprop('stop');
    end
    ax.stop = @(x)stop(x);
    
    if ~ax.isprop('update_fig')
        ax.addprop('update_fig');
    end    
    ax.update_fig = @(x)update_fig(x);

    
    %%% set handles
    
    % slider
    slider_handle=uicontrol(gcf,'Style','slider','Max',sz,'Min',1,...
        'Value',vars.init_frame,'SliderStep',[1/(sz-1) 10/(sz-1)],...
        'Units','normalized','Position',vars.slider_position);
    if ~isprop(slider_handle, 'ax_index')
        slider_handle.addprop('ax_index');
        slider_handle.ax_index = ax.index;
    end
    if ~ax.isprop('slider_handle')
        ax.addprop('slider_handle');
        ax.slider_handle = slider_handle;
    elseif ax.isprop('slider_handle') && ~ax.slider_handle.isvalid
        ax.slider_handle = slider_handle;
    end
    
    % play button
    if vars.show_play_button
        visible_button = 'on';
    else
        visible_button = 'off';
        if vars.loop
            warning('Loop can not be aborted without buttons. Setting ''loop'' back to ''false''.');
            vars.loop = false;
        end
    end
        
    play_handle=uicontrol(gcf,'Style','pushbutton','string','Play',...
        'Units','normalized','Position',vars.play_position, 'Visible', visible_button);
    if ~isprop(play_handle, 'ax_index')
        play_handle.addprop('ax_index');
        play_handle.ax_index = ax.index;
    end
    if ~ax.isprop('play_handle')
        ax.addprop('play_handle');
        ax.play_handle = play_handle;
    elseif ax.isprop('play_handle') && ~ax.play_handle.isvalid
        ax.play_handle = play_handle;
    end
    if ~ax.isprop('vars')
        ax.addprop('vars');
        ax.vars = vars;
    else
        ax.vars = vars;
    end
    set(play_handle,'Callback',{@play_callback,ax});

    
    % text edit
    if vars.show_edit_box
        visible_box = 'on';
    else
        visible_box = 'off';
    end
    edit_handle = uicontrol('style','edit','units','normalized', 'Position', vars.edit_position, 'Visible', visible_box);
    set(edit_handle, 'Callback', {@edit_callback, ax});
    if ~ax.isprop('edit_handle')
        ax.addprop('edit_handle');
        ax.edit_handle = edit_handle;
    elseif ax.isprop('edit_handle') && ~ax.edit_handle.isvalid
        ax.edit_handle = edit_handle;
    end    
    if ~isprop(edit_handle, 'ax_index')
        edit_handle.addprop('ax_index');
        edit_handle.ax_index = ax.index;
    end

    
    % set callback functions
    set(slider_handle,'Callback',{@slider_callback,ax});
    set(edit_handle, 'String', num2str(get(ax.slider_handle,'Value')));
    
    
    if vars.autoplay
        play_callback(ax, ax, ax);
    end
    
    update_fig(ax)
    
else
    % standard imagesc should be enough
    if ax.isprop('update_title') || ax.isprop('play_handle') || ax.isprop('vars')
        if ax.isprop('vars') && isfield(ax.vars, 'slider_handle')
            ax.vars = rmfield(ax.vars, 'slider_handle');
        end
        if ax.isprop('play_handle')
            delete(ax.play_handle);
        end
        cla(ax);
    end
    
    img = gather(vars.fnct(img));

    if ~isreal(img)
        img = c2image(img);
    end
    if ~isempty(ax_img)
        imagesc(ax_img{:}, img);
    else
        imagesc(img);
    end
    if ~isempty(vars.title_list)
        title(ax, vars.title_list{1}, 'Interpreter', 'none')
    end
end

end

% plotting function
function update_fig(ax)
import math.isint
import plotting.c2image

%     im = imhandles(gcf);
    vars = ax.vars;
    slice = round(get(ax.slider_handle,'Value'));
    slice = max(1, min(length(vars.order), slice));

    % FIXME: everything works better without following lines 
    % sl = gcbo();
    % if ~isempty(sl)
    %     ax = findobj('index', sl.ax_index);
    % end

    img = gather(vars.fnct(squeeze(ax.img(:,:,vars.order(slice),:))));
    
    if ~isreal(img)
        img = c2image(img);
    end
    if vars.plot_residua
        [residua{2},residua{1}] = find(abs(utils.findresidues(img))>0.1); 
        
    end
    
    % if the current axis is empty, use imagesc with remaining arguments
    if ~ax.isprop('update_title')
        ax.addprop('update_title');
        ax.addprop('user_title');
        ax.update_title = true;
        if ~isempty(vars.vargin)
            imagesc(vars.vargin{:}, img);
        else
            imagesc(img);
        end
        hold all 
        if vars.plot_residua && ~isempty(residua{1})
            plot(residua{:},'or')
        elseif vars.plot_residua
            plot(0,0,'or')
        end
        hold off 
        addlistener(ax.Title, 'String', 'PostSet', @(gt, event)callback_title_post(ax, ax));
    else
        % if we just need to update the figure, only update the data
        ax_data = ax.findobj('Type', 'Image');
        ax_data.CData = img;
        if vars.plot_residua
            ax_data = ax.findobj('Type', 'Line');
            ax_data(1).XData = residua{1};
            ax_data(1).YData = residua{2};
        end
    end    

        
    % write title
    if isempty(vars.title_list)
        ax.update_title = false;
        if ~isempty(ax.user_title)
            title_text = sprintf(ax.user_title, vars.order(slice));
            title(ax, title_text, 'Interpreter', 'none');
        end
        ax.update_title = true;
    else
        if ~isempty(vars.title_list)
            title(ax, vars.title_list{vars.order(slice)}, 'Interpreter', 'none')
        end
    end
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Callback subfunctions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider_callback(~,~,ax)
%     ob = gco;
%     vars.ax_index = ob.ax_index;
    if ax.isprop('edit_handle')
        set(ax.edit_handle, 'string', num2str(round(get(ax.slider_handle,'Value'))));
    end
    update_fig(ax)
    drawnow()

end


function callback_title_post(ax, ~, ~)
    if ax.update_title
        ax.user_title = ax.Title.String;
        try
        ax.Title.String = sprintf(ax.user_title, get(ax.slider_handle,'Value'));
        end
    end
end


function play_callback(~,~,ax)
%     ax = findobj('index', ax.slider_handle.ax_index);
    vars = ax.vars;
    if ax.isprop('slider_handle') && ax.slider_handle.isvalid
        update_slider = true;
    else
        update_slider = false;
    end
    if ax.isprop('edit_handle') && ax.edit_handle.isvalid
        update_edit = true;
    else
        update_edit = false;
    end
    try
        switch get(ax.play_handle,'string')
            case 'Play'
                if ~isempty(vars.save_movie)
                    disp(['Saving movie to  ' vars.save_movie]);
                    writeobj = VideoWriter(vars.save_movie);
                    writeobj.Quality=vars.movie_quality;
                    writeobj.FrameRate=vars.fps;
                    
                    open(writeobj);
                    vars.writeobj = writeobj;
                end
                set(ax.play_handle,'string','Stop')
                sz = size(ax.img,3);
                pos = round(get(ax.slider_handle,'Value'));
                if pos == sz
                    set(ax.slider_handle,'Value',1);
                    pos = 1;
                end
                while pos <=sz
                    
                    if strcmp(get(ax.play_handle,'string'), 'Play')
                        break
                    end
                    if update_slider
                        set(ax.slider_handle,'Value',pos)
                    end
                    if update_edit
                        set(ax.edit_handle, 'String', num2str(pos));
                    end
                    
                    update_fig(ax)
                    pause(1/vars.fps)
                    if vars.loop && pos == sz
                        pos = 1;
                    else
                        pos = pos+1;
                    end
                    if vars.save_movie
                        currFrame = getframe;
                        writeVideo(vars.writeobj,currFrame);
                    end
                    
                end
                set(ax.play_handle,'string','Play')
                if vars.reset_frame
                    set(ax.slider_handle,'Value',vars.init_frame)
                end
                if vars.save_movie
                    close(vars.writeobj);
                end
            case 'Stop'
                set(ax.play_handle,'string','Play')
                if vars.save_movie    
                    close(vars.writeobj); 
                end
        end
    catch
        if ~ax.isprop('play_handle')
            fprintf('Lost connection to figure instance.\n')
        end
    end
end

function edit_callback(~,~, ax)
    str=get(ax.edit_handle,'String');
    
    if isempty(str2num(str))
        warndlg('Input must be numerical'); 
        set(ax.edit_handle, 'string', num2str(round(get(ax.slider_handle,'Value'))));
    else
        set(ax.slider_handle,'Value',str2num(str))
        update_fig(ax)
        drawnow()
    end

end

function play(ax)
    play_callback(ax.vars, ax.vars, ax);
end

function stop(ax)
    set(ax.play_handle,'string','Play')
    set(ax.slider_handle,'Value',ax.vars.init_frame)
end