% Usage
% This function can be used standalone using "Load stack"
% Can also be run from a script with optional inputs
% [stack_phase_corr mask] = findmask(stack_object)
%   Inputs
% stack_object      - Stack of complex valued projections
%   Outputs
% stack_phase_corr  - Stack of corrected complex valued projections
% mask              - Mask used for correction
% Manuel Guizar, Cameron Kewish 2012-10-25

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



function varargout = findmask(varargin)
% FINDMASK MATLAB code for findmask.fig
%      FINDMASK, by itself, creates a new FINDMASK or raises the existing
%      singleton*.
%
%      H = FINDMASK returns the handle to a new FINDMASK or the handle to
%      the existing singleton*.
%
%      FINDMASK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINDMASK.M with the given input arguments.
%
%      FINDMASK('Property','Value',...) creates a new FINDMASK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before findmask_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to findmask_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help findmask

% Last Modified by GUIDE v2.5 25-Oct-2012 15:54:12

% Begin initialization code - DO NOT EDIT
 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @findmask_OpeningFcn, ...
                   'gui_OutputFcn',  @findmask_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before findmask is made visible.
function findmask_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to findmask (see VARARGIN)

% Choose default command line output for findmask
handles.output = hObject;
%%%% Initialization with variable %%%
if nargin >= 4
    handles.stack_object = varargin{1};
    if numel(varargin)>1
        handles.mask_ramp = varargin{2};
    end
    if ~isfield(handles,'mask_ramp')
        handles.mask_ramp = logical(handles.stack_object*0);
    else
        if ~all(size(handles.stack_object) == size(handles.mask_ramp))
            msgbox('Size of mask and projections does not match','Wrong mask','error')
        end
    end
    handles.corrected = ones(size(handles.stack_object,3),1);
    handles.see_mask_only = get(handles.checkbox_see_mask_only,'Value');
    handles.projnum = 1;
    handles.linecut = round(size(handles.stack_object,1)/2);
    clear stack_object
    set(handles.edit_line,'String',num2str(handles.linecut))
    handles.region_select = handles.stack_object(:,:,handles.projnum)*0;
    show_image(handles)
    %%% Make shape %%%
    set(handles.make_square,'Enable','on')
    set(handles.make_poly,'Enable','on')
    %%% Save %%%
    set(handles.save_stack,'Enable','on')
    set(handles.save_mask,'Enable','on')
    set(handles.return_to_script,'Enable','on')
    %%% Display %%%
    set(handles.edit1,'Enable','on')
    set(handles.prev_proj,'Enable','on')
    set(handles.stop_button,'Enable','on')
    set(handles.play,'Enable','on')
    set(handles.next_proj,'Enable','on')
    set(handles.edit_line,'Enable','on')
    set(handles.checkbox_see_mask_only,'Enable','on')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes findmask wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = findmask_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in load_stack.
function load_stack_Callback(hObject, eventdata, handles)
% hObject    handle to load_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.stack_filename handles.stack_path Filterindex] = uigetfile('*.mat');
load([handles.stack_path handles.stack_filename]);
handles.stack_object = stack_object;
if ~exist('stack_object')
    msgbox('Variable stack_object does not exist in the file','Variable not found','error')
end
if ~isfield(handles,'mask_ramp')
    handles.mask_ramp = logical(stack_object*0);
else
    if ~all(size(handles.stack_object) == size(handles.mask_ramp))
        msgbox('Size of mask and projections does not match','Wrong mask','error')
    end
end
handles.corrected = ones(size(stack_object,3),1);
handles.see_mask_only = get(handles.checkbox_see_mask_only,'Value');
handles.projnum = 1;
handles.linecut = round(size(stack_object,1)/2);
clear stack_object
set(handles.edit_line,'String',num2str(handles.linecut))
handles.region_select = handles.stack_object(:,:,handles.projnum)*0;
show_image(handles)
%%% Make shape %%%
set(handles.make_square,'Enable','on')
set(handles.make_poly,'Enable','on')
%%% Save %%%
set(handles.save_stack,'Enable','on')
set(handles.save_mask,'Enable','on')
%%% Display %%% 
set(handles.edit1,'Enable','on')
set(handles.prev_proj,'Enable','on')
set(handles.stop_button,'Enable','on')
set(handles.play,'Enable','on')
set(handles.next_proj,'Enable','on')
set(handles.edit_line,'Enable','on')
set(handles.checkbox_see_mask_only,'Enable','on')
%%%%%%%%%%%%%%%
% save the changes to the structure
guidata(hObject,handles)





% --- Executes on button press in load_mask.
function load_mask_Callback(hObject, eventdata, handles)
% hObject    handle to load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.mask_filename handles.mask_path Filterindex] = uigetfile('*.mat');
load([handles.mask_path handles.mask_filename]);
if ~exist('mask')
    msgbox('Variable mask does not exist in the file','Variable not found','error')
end
handles.mask_ramp = mask;
if isfield(handles,'stack_object')
    if ~all(size(handles.stack_object) == size(handles.mask_ramp))
        msgbox('Size of mask and projections do not match','Wrong mask','error')
    end
end
handles.corrected = zeros(size(handles.mask_ramp,3),1);
handles.projnum = 1;
clear mask_ramp
handles.region_select = handles.mask_ramp(:,:,handles.projnum)*0;
if isfield(handles,'stack_object')
    show_image(handles)
end
enable_mask_panel(hObject, handles)
enable_ramp_panel(hObject, handles)
% save the changes to the structure
guidata(hObject,handles)



% --- Executes on button press in save_stack.
function save_stack_Callback(hObject, eventdata, handles)
% hObject    handle to save_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if any(handles.corrected == 0)
    msgbox('Mask has not been applied to all projections. Not saving.','Current mask not applied yet','warn')
else
    stack_object = handles.stack_object;
    if isfield(handles,'stack_path')
        uisave('stack_object',[handles.stack_path 'phase_corr_' handles.stack_filename]);
    else
        uisave('stack_object',['phase_corr_.mat']);
    end
end



% --- Executes on button press in save_mask.
function save_mask_Callback(hObject, eventdata, handles)
% hObject    handle to save_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mask = handles.mask_ramp;
    if isfield(handles,'stack_path')
        uisave('mask',[handles.stack_path 'mask_' handles.stack_filename]);
    else
        uisave('mask',['mask_.mat']);
    end



% --- Executes on button press in return_to_script.
function varargout = return_to_script_Callback(hObject, eventdata, handles)
% hObject    handle to return_to_script (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if any(handles.corrected == 0)
    msgbox('Mask has not been applied to all projections. Not ready to return to script.','Current mask not applied yet','warn')
else
     msgbox('Load variable from script and close figure.','Exit to script','help')
end


%%%%%%%%%%%%%%%%%%%
%%%  Make mask  %%%
%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in make_square.
function make_square_Callback(hObject, eventdata, handles)
% hObject    handle to make_square (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure_text,'string','Click on two corners of the desired rectangle',...
    'foregroundcolor',[0 0 1])
set(gcf,'CurrentAxes',handles.axes_phase)
[x y] = ginput(2);
hold on,
plot([x(1) x(1)],[y(1) y(2)])
plot([x(2) x(2)],[y(1) y(2)])
plot([x(1) x(2)],[y(1) y(1)])
plot([x(1) x(2)],[y(2) y(2)])
hold off,
x = round(x);
y = round(y);
handles.region_select = handles.region_select*0;
handles.region_select(min(y):max(y),min(x):max(x)) = 1;
set(handles.figure_text,'string','Choose to add or remove from mask',...
    'foregroundcolor',[0 0 1])
enable_mask_panel(hObject, handles)
enable_ramp_panel(hObject, handles)
% save the changes to the structure
guidata(hObject,handles)


% --- Executes on button press in make_poly.
function make_poly_Callback(hObject, eventdata, handles)
% hObject    handle to make_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h = msgbox(['Please click on the image to choose the corners ' ...
%     'of the polygon ROI. When finished, double-click the last vertex to ' ...
%     'close the polygon.'], 'Make Poly', 'help');
set(handles.figure_text,'string','Click closed polygon corners. Double-click when done',...
    'foregroundcolor',[0 0 1])
set(gcf,'CurrentAxes',handles.axes_phase)
[handles.region_select x y] = roipoly;
% do we want to plot the roi border?
hold on,
plot(x,y)
hold off,
set(handles.figure_text,'string','Choose to add or remove from mask',...
    'foregroundcolor',[0 0 1])
enable_mask_panel(hObject, handles)
enable_ramp_panel(hObject, handles)
% save the changes to the structure
guidata(hObject,handles)



% --- Executes on button press in next_proj.
function next_proj_Callback(hObject, eventdata, handles)
% hObject    handle to next_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.projnum = min(handles.projnum + 1,size(handles.stack_object,3));
show_image(handles);
% save the changes to the structure
guidata(hObject,handles);

    

% --- Executes on button press in prev_proj.
function prev_proj_Callback(hObject, eventdata, handles)
% hObject    handle to prev_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.projnum = max(handles.projnum - 1,1);
show_image(handles);
% save the changes to the structure
guidata(hObject,handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.projnum = str2double(get(hObject,'String'));
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% persistent playmovie
% playmovie = true;
ii = handles.projnum;
handles.playmovie = true;
guidata(hObject,handles);
% % while playmovie
while handles.playmovie
    handles = guidata(hObject);
%     display(handles.playmovie)
    handles.projnum = ii;
    if ii == size(handles.stack_object,3)
        ii = 1;
    else
        ii = ii + 1;
    end
    show_image(handles);
    pause(0.2)
end
% save the changes to the structure
guidata(hObject,handles);


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% persistent playmovie
% playmovie = false;
handles.playmovie = false;
% display('Trying to stop movie')
% save the changes to the structure
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Mask Operations  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in add_to_mask.
function add_to_mask_Callback(hObject, eventdata, handles)
% hObject    handle to add_to_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask_ramp(:,:,handles.projnum) = handles.mask_ramp(:,:,handles.projnum)|handles.region_select;
handles.corrected(handles.projnum) = 0;
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)



% --- Executes on button press in remove_from_mask.
function remove_from_mask_Callback(hObject, eventdata, handles)
% hObject    handle to remove_from_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask_ramp(:,:,handles.projnum) = handles.mask_ramp(:,:,handles.projnum)&(1-handles.region_select);
handles.corrected(handles.projnum) = 0;
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)



% --- Executes on button press in add_toall_masks.
function add_toall_masks_Callback(hObject, eventdata, handles)
% hObject    handle to add_toall_masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for ii = 1:size(handles.stack_object,3)
    handles.mask_ramp(:,:,ii) = handles.mask_ramp(:,:,ii)|handles.region_select;
end
handles.corrected = zeros(size(handles.stack_object,3),1);
    
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)



% --- Executes on button press in remove_fromall_masks.
function remove_fromall_masks_Callback(hObject, eventdata, handles)
% hObject    handle to remove_fromall_masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for ii = 1:size(handles.stack_object,3)
    handles.mask_ramp(:,:,ii) = handles.mask_ramp(:,:,ii)&(1-handles.region_select);
end
handles.corrected = zeros(size(handles.stack_object,3),1);
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%
%%%  Ramp removal  %%%
%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in apply_removal_current.
function apply_removal_current_Callback(hObject, eventdata, handles)
% hObject    handle to apply_removal_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure_text,'String',['Removing offset on projection ' ...
    num2str(handles.projnum) ' ...'],...
    'foregroundcolor',[0 0 1])
upsamp = 100;
% [handles.stack_object(:,:,handles.projnum) errorm] = utils.remove_linearphase(handles.stack_object(:,:,handles.projnum),handles.mask_ramp(:,:,handles.projnum),upsamp);
[handles.stack_object(:,:,handles.projnum) errorm] = utils.remove_linearphase(exp(1i*angle(handles.stack_object(:,:,handles.projnum))),handles.mask_ramp(:,:,handles.projnum),upsamp);
handles.corrected(handles.projnum) = 1;
% save the changes to the structure
guidata(hObject,handles)
show_image(handles)

% --- Executes on button press in apply_removal_all.
function apply_removal_all_Callback(hObject, eventdata, handles)
% hObject    handle to apply_removal_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
upsamp = 100;
hwait = waitbar(0,'Patience is a virtue');

tmp = handles.stack_object;
tmpmask = handles.mask_ramp;

% CMK: Added a parfor loop to save time here, but it broke the "Patience" waitbar...
% if isempty(gcp('nocreate')) == 0
%     parpool
% end

parfor ii = 1:size(tmp,3)
    [tmp(:,:,ii) errorm] = utils.remove_linearphase(tmp(:,:,ii),tmpmask(:,:,ii), upsamp);
    display(ii)
end
handles.stack_object = tmp;

% for ii = 1:size(handles.stack_object,3)
%     set(handles.figure_text,'String',['Removing offset on projection ' ...
%         num2str(ii) ' ...'],...
%         'foregroundcolor',[0 0 1])
%     guidata(hObject,handles)
%     if ~handles.corrected(ii)
%         [handles.stack_object(:,:,ii) errorm] = utils.remove_linearphase(handles.stack_object(:,:,ii),handles.mask_ramp(:,:,ii),upsamp);
%         handles.corrected(ii) == 1;
%     end
%     waitbar(ii/size(handles.stack_object,3),hwait);
% end

close(hwait)
handles.corrected = ones(size(handles.stack_object,3),1);
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)



function edit_line_Callback(hObject, eventdata, handles)
% hObject    handle to edit_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_line as text
%        str2double(get(hObject,'String')) returns contents of edit_line as a double
handles.linecut = str2double(get(hObject,'String'));
show_image(handles)
% save the changes to the structure
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_see_mask_only.
function checkbox_see_mask_only_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_see_mask_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.see_mask_only = get(handles.checkbox_see_mask_only,'Value');
show_image(handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox_see_mask_only
% save the changes to the structure
guidata(hObject,handles)

function enable_mask_panel(hObject, handles)
set(handles.add_to_mask,         'Enable','on')
set(handles.add_toall_masks,     'Enable','on')
set(handles.remove_from_mask,    'Enable','on')
set(handles.remove_fromall_masks,'Enable','on')
% save the changes to the structure
guidata(hObject,handles)

function enable_ramp_panel(hObject, handles)
set(handles.apply_removal_current,'Enable','on')
set(handles.apply_removal_all,    'Enable','on')

% save the changes to the structure
guidata(hObject,handles)

function show_image(handles)
phase_plot = angle(handles.stack_object(:,:,handles.projnum));
min_phaseplot = min(phase_plot(:));
max_phaseplot = max(phase_plot(:));
phase_plot = phase_plot - min_phaseplot;
phase_plot = phase_plot.*(handles.mask_ramp(:,:,handles.projnum) + ...
    0.5*(~handles.mask_ramp(:,:,handles.projnum)))+min_phaseplot;
% handles.axes_phase = imagesc(phase_plot);
%%%
set(gcf,'CurrentAxes',handles.axes_phase);
if ~handles.see_mask_only
    imagesc(phase_plot);
    caxis([min_phaseplot max_phaseplot]);
else
    imagesc(phase_plot.*handles.mask_ramp(:,:,handles.projnum));
end
axis xy equal tight
colormap bone
%%%
set(gcf,'CurrentAxes',handles.axes_linecut);
plot(angle(handles.stack_object(handles.linecut,:,handles.projnum)));
xlim([1 size(handles.stack_object,2)]);
grid on
%%%
set(handles.figure_text,'String', ...
    ['Projection ' num2str(handles.projnum)],'foregroundcolor',[0 0 0])







%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auxiliary functions  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    output=[error,diffphase];
        
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc); 
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    md2 = fix(m/2); 
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
% Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
  
    % Compute crosscorrelation and locate the peak 
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid 
        [max1,loc1] = max(CC);   
        [max2,loc2] = max(max1); 
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end  

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return



