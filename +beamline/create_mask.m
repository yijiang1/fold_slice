%CREATE_MASK
% create a binary mask for the current figure
% The following arguments have to be given as name/value pairs. However,
% they can also be set within the GUI.
% 
% *optional*
% ** mask               initial mask; either a file, an array or a structure (indicies + asize)
% ** fig                pass figure handle; default: current figure
% ** ind                convert mask to indicies
% ** file               save mask to disk; specify path + filename
% 
% returns:
% ++ out                2D binary mask or structure containing the asize and the indicies
%
% EXAMPLE:
%   img = io.image_read('~/Data10/pilatus_1/S00000-00999/S00170/*.cbf'); % load image stack
%   plotting.imagesc3D(log10(img.data)); axis equal tight xy; % plot image stack
%   beamline.create_mask(); % open the GUI and create the mask
%
% see also: beamline.mask2ind

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [out] = create_mask(varargin)

check_input_mask = @(x) ischar(x) || (isnumeric(x)|| islogical(x)) || isstruct(x);

par = inputParser;
par.addParameter('mask', [], check_input_mask)
par.addParameter('fig', [], @ishandle)
par.addParameter('ind', false, @islogical)
par.addParameter('file', [], @ischar)
par.parse(varargin{:})
vars = par.Results;

% Check screen size
try
    scrsz = get(0,'ScreenSize');
catch
    scrsz = [1 1 2560 1024];
end

% get fig
if isempty(vars.fig)
    fig = gcf;
end

current_pos = fig.Position;
new_fig_pos(2:4) = current_pos(2:4);

if current_pos(1)+current_pos(3)/2 - scrsz(3)/2 > 0
    % figure to the left
    new_fig_pos(1) = current_pos(1)-current_pos(3);
else
    % figure to the right
    new_fig_pos(1) = current_pos(1)+current_pos(3);
end


% get axis
ax = gca; 
% get current data size
if ~isempty(ax.Children)
    asize = size(ax.Children.CData);
else
    fig = gcf;
    close(fig)
    error('Failed to connect to figure instance.')
end

% prepare mask
if isempty(vars.mask)
    mask = ones(asize);
else
    if ischar(vars.mask)
        % load a mask from disk
        try
            f = load(vars.mask);
            mask = f.mask;
            clear f
        catch
            fprintf('Failed to load mask. Using empty mask instead.\n')
            mask = ones(asize);
        end
    elseif isnumeric(vars.mask) || islogical(vars.mask)
        mask = vars.mask;
    elseif isstruct(vars.mask)
        mask = beamline.ind2mask(vars.mask);
    else
        error('Unknown mask data format.')
    end
    assert(all(size(mask)==asize), 'Mask size and data size does not match')

end
pause(0.1)


% apply mask
CData_orig = ax.Children.CData;
if ax.isprop('img')
    img_orig = ax.img;
    ax.img = ax.img .* mask;
    mask_dims = ndims(img_orig);
    if mask_dims==3
        mask3D = true;
    else
        mask3D = false;
    end
    mask_dims = size(img_orig,3);
else
    mask3D = false;
    mask_dims = 1;
end
ax.Children.CData = ax.Children.CData .* mask;

s = create_mask_GUI_export('mask', mask, 'mask3D', mask3D, 'mask_dims', mask_dims);
s.figure1.UserData.ax = ax;
s.figure1.UserData.fig = fig;
s.figure1.UserData.asize = asize;
s.figure1.UserData.CData = CData_orig;
if ax.isprop('img')
    s.figure1.UserData.img_orig = img_orig;
    if s.figure1.UserData.mask3D
        orig_fig_listener = s.figure1.UserData.ax.slider_handle.listener('Value','PostSet',@(src, evnt)orig_fig_slice_update(s));
    end

end

try
    while ~s.figure1.UserData.done
        mask = s.figure1.UserData.mask;
        pause(0.1)
    end
    
    set(groot,'CurrentFigure',fig);
    ax.Children.CData = CData_orig;
    if ax.isprop('img')
        ax.img = img_orig;
    end
catch
    if ~isprop(s, 'figure1')
        fprintf('Lost connection to GUI.\n')
    end
end

try
    if s.figure1.UserData.mask3D
        delete(orig_fig_listener)
    end
    delete(s.figure1)
catch
end


% if needed, convert 2D mask to indicies
if vars.ind
    out = beamline.mask2ind(mask);
else
    out = mask;
end                



% save to disk
if ~isempty(vars.file)
    valid_mask = out;
    try
        utils.savefast_safe(vars.file, 'valid_mask');
    catch
        fprintf('Failed to save mask to disk.');
    end
end


end

function orig_fig_slice_update(s)
val = s.figure1.UserData.ax.slider_handle.Value;
set(s.axes1.slider_handle, 'Value', val);
set(s.axes1.edit_handle, 'String', num2str(val));
s.axes1.update_fig(s.axes1);

end

