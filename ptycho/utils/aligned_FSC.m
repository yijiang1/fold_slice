%   [resolution stat] = aligned_FSC(file1,file2,param)
%
% Receives two filenames with path for ptychography reconstructions and a
% structure with parameters. The routine reads the reconstructions, matches
% the linear phase between them, registers the images, and returns the
% resolution estimates based on first and last crossing of the FSC with the
% threshold.
%
% References relevant to this code:
% For using this FSC code with ptychography: J. Vila-Comamala, et al., "Characterization of high-resolution diffractive X-ray optics by ptychographic coherent diffractive imaging," Opt. Express 19, 21333-21344 (2011).
% For subpixel alignment: M. Guizar-Sicairos, et al., "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156 (2008).
% For matching of phase ramp by approximate least squared error:  M. Guizar-Sicairos, et al., "Phase tomography from x-ray coherent diffractive imaging projections," Opt. Express 19, 21345-21357 (2011).
%
% Outputs:
%
% resolution    A two element variable that contains the resolution
%               obtained from first and last crossing of the FSC curve with
%               the threshold curve.
% stat          structure containing other statistics such as
%                  spectral signal to noise ratio (SSNR), average SNR and area under FSC curve
%
% Inputs:
%
% file1     Filename with path of reconstruction 1 or directly a 2D numerical array 
% file2     Filename with path of reconstruction 2 or directly a 2D numerical array 
% param    Structure with parameters as describred below
%
%  param.flipped_images    Flip one input image horizontally (= true or false).
%                           Useful when comparing 0 and 180 degree projections 
%                           in tomography (default = false).
%  param.crop      = '';       for using the default half size of the probe
%                   = 'manual'  for using GUI to select region. This will display the range, e.g. {600:800, 600:800}     
%                   = {600:800, 600:800} for custom vertical and horizontal cropping, respectively
%  param.GUIguess  To click for an initial alignment guess, if used it ignores 
%                   the values of param.guessx and param.guessy (default
%                   = false)
%  param.guessx
%  param.guessy        An intial guess for x and y alignment (default = [])
%  param.remove_ramp   Try to remove linear phase from whole image before initial 
%                       alignment (default = true)
%  param.image_prop    = 'complex' 
%                       = 'phasor' (phase with unit amplitude, default) 
%                       = 'phase'  (Note: phase should not be used if there is phase wrapping)
%  param.taper         = 20 (default)   Pixels to taper images - Increase until the FSC does not change anymore
%  param.plotting      Display plots (default = false)
%  param.dispfsc       Display FSC plot (default = true)
%  param.SNRt          SNR for FSC threshold curve   
%                       SNRt = 0.2071 for 1/2 bit threshold for resolution of the average of the 2 images
%                       SNRt = 0.5    for 1   bit threshold for resolution of each individual image (default)
%  param.thickring     Thickness of Fourier domain ring for FSC in pixels (default = 1)
%  param.freq_thr      To ignore the crossings before freq_thr for determining resolution (default 0.02)
%  param.out_fn        Filename of output of jpeg for FSC
%  param.pixel_size    Pixel size in the reconstruction, it is used only
%                       if file1 / file2 are not paths to the reconsturcted files  

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

function [resolution,stat] = aligned_FSC(file1,file2,param)
import utils.*
import math.*
import io.*
import plotting.*

file{1} = file1;
file{2} = file2;

for ii = 1:2
    if ischar(file{ii})
        [~,file_path{ii},~] =  fileparts(file{ii});
    else
        file_path{ii} = sprintf('image_id_%i', ii); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checks and defaults %%%
if isfield(param,'flag_imread')
    flag_imread = param.flag_imread;
else
    if ~(ischar(file1) && ischar(file2))  || (~isempty(regexpi(file1,'\.mat|\.h5')) && ~isempty(regexpi(file2,'\.mat|\.h5')))
        flag_imread = false;
    else
        flag_imread = true;
        warning('Files not in .mat: using ''imread'' for loading; image in real number (ignoring remove_ramp and image_prop)')
    end
end


% parse inputs
check_input = @(x) islogical(x) || isnumeric(x);
check_crop = @(x) isempty(x) || ischar(x) || iscell(x);
check_image_prop = @(x) assert(any(contains({'complex', 'phasor', 'phase', 'variation'}, x)), ...
    'image_prop must be either "complex", "phasor", "phase" or "variation".');
parse_param = inputParser;
parse_param.KeepUnmatched = true;

parse_param.addParameter('flipped_images', false, check_input)
parse_param.addParameter('crop', '', check_crop)
parse_param.addParameter('GUIguess', false, check_input)
parse_param.addParameter('guessx', [], check_input)
parse_param.addParameter('guessy', [], check_input)
parse_param.addParameter('plotting',false, check_input)
parse_param.addParameter('remove_ramp', false, check_input)
parse_param.addParameter('image_prop', 'phasor', check_image_prop)
parse_param.addParameter('SNRt', 0.5, @isnumeric)
parse_param.addParameter('thickring', 1, @isnumeric)
parse_param.addParameter('freq_thr', 0.02, @isnumeric)
parse_param.addParameter('prop_obj', [], check_input)
parse_param.addParameter('apod', [], check_input)
parse_param.addParameter('filter_FFT', [], check_input)
parse_param.addParameter('crop_factor', 1, @isnumeric)
parse_param.addParameter('crop_asize', [], @isnumeric)
parse_param.addParameter('z_lens', [], @isnumeric)
parse_param.addParameter('fourier_ptycho', false, check_input)
parse_param.addParameter('lambda', [], @isnumeric)
parse_param.addParameter('pixel_size', [], @isnumeric)
parse_param.addParameter('verbose_level', 3, @isnumeric)
parse_param.addParameter('fname', [], @iscell)
parse_param.addParameter('show_summary', true, check_input)

parse_param.parse(param)
param = parse_param.Results;

if isempty(param.crop)
    utils.verbose(3,'Cropping half size of the probe (default)')
end

if isfield(param,'taper')
    taper = param.taper;
else
    taper = 20;
    utils.verbose(3,'Using taper = 20 (default)')
end

if isempty(param.fname)
    param.fname = file_path;
end

utils.verbose(param.verbose_level);
utils.verbose(struct('prefix', {'FSC'}))

% set dispfsc value - used in utils.fourier_shell_corr_3D_2
param.dispfsc = (param.plotting > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(file{1})
    if ~flag_imread
        %%% Checking if file exist and loading %%%
        if exist(file{1},'file')
            param.fourier_ptycho = io.HDF.hdf5_load(file{1}, '/reconstruction/p/fourier_ptycho');
            utils.verbose(2, ['Loading: ' file{1}])
            recons{1} = load_ptycho_recons(file{1}, 'object');
            if param.fourier_ptycho
                recons{1} = FP_FSC_preprocess(recons{1}, file{1}, param);
            end
        else
            error(['Not found: ' file{1}])
        end
        if exist(file{2},'file')
            utils.verbose(2, ['Loading: ' file{2}])
            recons{2} = load_ptycho_recons(file{2}, 'object');
            if param.fourier_ptycho
                recons{2} = FP_FSC_preprocess(recons{2}, file{2}, param);
            end
        else
            error(['Not found: ' file{2}])
        end

        img1 = recons{1}.object;
        img2 = recons{2}.object;

        % load additional parameters
        asize = double(io.HDF.hdf5_load(file{1}, '/reconstruction/p/asize')); 
        if ~param.fourier_ptycho
            param.pixel_size = io.HDF.hdf5_load(file{1}, '/reconstruction/p/dx_spec');
        else
            param.pixel_size = recons{1}.p.objpix;
        end
        try
        if isempty(param.lambda) && io.HDF.hdf5_dset_exists(file{1}, 'lambda', '/reconstruction/p', true)
            param.lambda = io.HDF.hdf5_load(file{1}, '/reconstruction/p/lambda');
        end
        end
    else
        if ~isnumeric(file{1})
            img1 = imread(file{1});
        else
            img1 = file{1};  %input provides directly the numeric array 
        end
        if ~isnumeric(file{2})
            img2 = imread(file{2});
        else
            img2 = file{2};  %input provides directly the numeric array 
        end
        asize = [1 1];
        if isempty(param.pixel_size)
            param.pixel_size = 1e-6;
            warning('Using pixel size 1um')
        end
    end
else
    img1 = file{1};  %input provides directly the numeric array 
    img2 = file{2};  %input provides directly the numeric array 
    if ~isfield(param, 'asize')
        asize = [1 1];
    else
        asize = param.asize;
    end
    if isempty(param.pixel_size)
        param.pixel_size = 1e-6;
        warning('Using pixel size 1um')
    end
end


if param.flipped_images
    img2 = fliplr(img2);
end

% apply apodization

% check if apodization was used in the reconstruction
if isempty(param.apod)
    try 
        param.apod = io.HDF.hdf5_load(file{1}, '/reconstruction/p/plot/obj_apod');
    catch
        warning('Unable to load apodization parameter from reconstruction file.')
        param.apod = false;
    end
end
if param.apod
    img1 = apply_apod(img1, asize);
    img2 = apply_apod(img2, asize);
end

% propagate if needed
if isempty(param.prop_obj) || param.prop_obj
    if isempty(param.lambda)
        try
            param.lambda = io.HDF.hdf5_load(file{1}, '/reconstruction/p/lambda');
        catch
            try
                param.energy = io.HDF.hdf5_load(file{1}, '/reconstruction/p/energy');
                param.lambda = 12.4/param.energy*1e-10;
            catch
                error('Please specify your wavelength (param.lambda).')
            end
        end
    end
    if isempty(param.prop_obj) || islogical(param.prop_obj) && param.prop_obj
        % get values from file
        try
            param.prop_obj = io.HDF.hdf5_load(file{1}, '/reconstruction/p/prop_obj');
        catch
            error('Please specify the propagation distance or set it to "false" (param.prop_obj).')
        end
    end 
    
    img1 = utils.prop_free_nf(img1, param.lambda, param.prop_obj, param.pixel_size);
    img2 = utils.prop_free_nf(img2, param.lambda, param.prop_obj, param.pixel_size);
end


screensize = get( 0, 'Screensize' ); 
% Show phase images (not cropped)%
if param.plotting > 1
    plotting.smart_figure(21)
    set(gcf,'Outerposition',[1 screensize(4)-550 500 500])    %[left, bottom, width, height
    if ~isreal(img1)
        imagesc(angle(img1), math.sp_quantile(angle(img1),[1e-2, 1-1e-2],10));
        if ~param.fourier_ptycho
            rectangle('Position',[asize([2,1])/2, [size(img1,2),size(img1,1)]-asize([2,1])], 'EdgeColor', 'red')
        end
    else
        imagesc(img1,  math.sp_quantile((img1),[1e-2, 1-1e-2],10));
    end
    axis xy equal tight
    colormap bone
    colorbar
    if param.prop_obj
        [si_unit, val] = utils.get_unit_length(param.prop_obj);
        title(sprintf([param.fname{1} '\npropagated by %d %s'], val, si_unit),'interpreter','none')
    else
        title(param.fname{1},'interpreter','none')
    end
    plotting.smart_figure(22)
    if ~isreal(img1)
        imagesc(angle(img2), math.sp_quantile(angle(img2),[1e-2, 1-1e-2],10));
    else
        imagesc(img2, math.sp_quantile(img2,[1e-2, 1-1e-2],10));
    end
    if ~param.fourier_ptycho
        rectangle('Position',[asize([2,1])/2, [size(img2,2),size(img2,1)]-asize([2,1])], 'EdgeColor', 'red')
    end
    axis xy equal tight
    colormap bone
    colorbar
    if param.prop_obj
        [si_unit, val] = utils.get_unit_length(param.prop_obj);
        title(sprintf([param.fname{2} '\npropagated by %d %s'], val, si_unit),'interpreter','none')
    else
        title(param.fname{2},'interpreter','none')
    end
    set(gcf,'Outerposition',[500 screensize(4)-550 500 500])    %[left, bottom, width, height
end

% Crop images - default is half the size of the probe on each side plus
% whatever needed to make them of equal size
if isempty(param.crop)
    minsize = min(size(img1),size(img2))-asize;
    img1 = crop_pad(img1,minsize ); 
    img2 = crop_pad(img2,minsize );
elseif strcmpi(param.crop, 'manual')
    figure()
    imagesc(angle(img1), math.sp_quantile(angle(img1),[1e-2, 1-1e-2],10));
    if ~param.fourier_ptycho
        rectangle('Position',[asize([2,1])/2, [size(img1,2),size(img1,1)]-asize([2,1])], 'EdgeColor', 'red')
    end
    colormap bone 
    axis image xy 
    title('Select compared region')
    disp('Manually select the compared region ... ')
    rect = round(getrect);
    param.crop = {rect(2)+(1:rect(4)),rect(1)+(1:rect(3))};
    disp('===========================')
    fprintf('Selected region: {%i:%i,%i:%i}\n',rect(2), rect(2)+rect(4), rect(1), rect(1)+rect(3));
    disp('===========================')
    pause(1)
end
if ~isempty(param.crop)
    img1 = img1(param.crop{:});
    img2 = img2(param.crop{:});
end

if param.GUIguess
    plotting.smart_figure(21)
    disp(['Click on a feature on figure 2'])
    [xin yin] = ginput(1);
    plotting.smart_figure(22)
    disp(['Click on a feature on figure 3'])
    [xin2 yin2] = ginput(1);
    param.guessx = round(xin-xin2);
    param.guessy = round(yin-yin2);
end

if ~isempty(param.guessx)
    switch sign(param.guessx)
        case 1
            img1 = img1(:,1+param.guessx:end);
            img2 = img2(:,1:end-param.guessx);
        case -1
            img1 = img1(:,1:end+param.guessx);
            img2 = img2(:,1-param.guessx:end);
    end
end
if ~isempty(param.guessy)
    switch sign(param.guessy)
        case 1
            img1 = img1(1+param.guessy:end,:);
            img2 = img2(1:end-param.guessy,:);
        case -1
            img1 = img1(1:end+param.guessy,:);
            img2 = img2(1-param.guessy:end,:);
    end
end


% Remove ramp 
if param.remove_ramp
    utils.verbose(3,'Removing ramp for initial alignment')
    img1 = utils.stabilize_phase(img1,'binning', 4);
    img2 = utils.stabilize_phase(img2, img1, 'binning', 4);
end
    
if param.plotting >2 
    plotting.smart_figure(23)
    set(gcf,'Outerposition',[1 1 500 476])    %[left, bottom, width, height
    if ~isreal(img1)
        imagesc(angle(img1));
    else
        imagesc(img1); 
    end
    axis xy equal tight
    colormap bone; colorbar
    ax23 = gca; 
    title(param.fname{1},'interpreter','none')
    plotting.smart_figure(24); 
    if ~isreal(img2)
        imagesc(angle(img2));
    else
        imagesc(img2);
    end
    axis xy equal tight
    colormap bone; colorbar
    title(param.fname{2},'interpreter','none')
    set(gcf,'Outerposition',[500 1 500 476])    %[left, bottom, width, height]
    ax34 = gca; 
    linkaxes([ax23, ax34], 'xy')
end



% img1 = img1 - utils.imgaussfilt2_fft(img1, 20); 
% img2 = img2 - utils.imgaussfilt2_fft(img2, 20); 
% image_prop= 'complex';

%%% Initial alignment %%%
utils.verbose(2,'Initial alignment')  
if ~flag_imread
    switch lower(param.image_prop)
        case  'complex'
            imgalign1 = img1;
            imgalign2 = img2;
            utils.verbose(2,'Registering complex valued images')
        case 'phasor'
            imgalign1 = ones(size(img1)).*exp(1i*angle(img1));
            imgalign2 = ones(size(img2)).*exp(1i*angle(img2));
            utils.verbose(2,'Registering phasor of complex valued images')
        case 'phase'
            imgalign1 = angle(img1);
            imgalign2 = angle(img2);
            utils.verbose(2,'Registering phase of complex valued images')
        case 'variation'
            [dX,dY] = math.get_phase_gradient_2D(img1); 
            imgalign1 = sqrt(dX.^2+dY.^2); 
            [dX,dY] = math.get_phase_gradient_2D(img2); 
            imgalign2 = sqrt(dX.^2+dY.^2); 
            
    end
else
    imgalign1 = img1;
    imgalign2 = img2; 
end

upsamp = 100;
displ = utils.verbose>3;
W = 1;
x1 = [];%[1:150];
x2 = x1;
y1 = [];%[1:238];
y2 = y1;



% imgalign2 = shiftpp2(imgalign2,10,-10); % To test range adjustment
[subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(imgalign1,imgalign2, x1, y1, x2, y2, upsamp, displ,1);

%%% Fine alignment (second round) %%%

% Remove ramp for fine alignment
utils.verbose(2,'Removing ramp for fine alignment')
%%% A patch for deltafine large
if max(regionsout.y2+round(delta(1)))>size(img2,1)
    warning('First subpixel registration refinement found large values')
    regionsout.y2 = [min(regionsout.y2):size(img2,1)-round(delta(1))];
    regionsout.y1 = regionsout.y2;
end
if max(regionsout.x2+round(delta(2)))>size(img2,2)
    warning('First subpixel registration refinement found large values')
    regionsout.x2 = [min(regionsout.x2):size(img2,2)-round(delta(2))];
    regionsout.x1 = regionsout.x2;
end
%%%
subimg1 = img1(regionsout.y1,regionsout.x1);
subimg2 = img2(regionsout.y2+round(delta(1)),regionsout.x2+round(delta(2)));
if ~flag_imread
    subimg1 = remove_linearphase_v2(subimg1,ones(size(subimg1)),100);
    subimg2 = remove_linearphase_v2(subimg2,ones(size(subimg2)),100);
end


% Remove ramp 
if param.remove_ramp
    utils.verbose(2,'Removing ramp for initial alignment')
    subimg1 = utils.stabilize_phase(subimg1,'binning', 4);
    subimg2 = utils.stabilize_phase(subimg2, subimg1, 'binning', 4);
end

if ~flag_imread
    switch lower(param.image_prop)
        case  'complex'
            subimgalign1 = subimg1;
            subimgalign2 = subimg2;
            utils.verbose(2,'Registering complex valued images')
        case 'phasor'
            subimgalign1 = ones(size(subimg1)).*exp(1i*angle(subimg1));
            subimgalign2 = ones(size(subimg1)).*exp(1i*angle(subimg2));
            utils.verbose(2,'Registering phasor of complex valued images')
        case 'phase'
            subimgalign1 = angle(subimg1);
            subimgalign2 = angle(subimg2);
            utils.verbose(2,'Registering phase of complex valued images')
        case 'variation'
            [dX,dY] = math.get_phase_gradient_2D(subimg1); 
            subimgalign1 = sqrt(dX.^2+dY.^2); 
            [dX,dY] = math.get_phase_gradient_2D(subimg2); 
            subimgalign2 = sqrt(dX.^2+dY.^2); 
    end
else
   subimgalign1 = subimg1;
   subimgalign2 = subimg2; 
end

% Fine alignment %
utils.verbose(2,'Fine alignment')
[subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(subimgalign1,subimgalign2, x1, y1, x2, y2, upsamp, displ,1);

%%% propare images for FSC if variation was used for alignement 
if strcmpi(param.image_prop, 'variation')
    subim1 = subimg1(regionsout.y1, regionsout.x1); 
    subim2 = subimg2(regionsout.y2, regionsout.x2); 
    subim2 = shiftpp2(subim2,-deltafine(1), -deltafine(2)); %% Suboptimal, change to use a routine that receives FT data
    % convert images to phasor 
    subim1 = exp(1i*angle(subim1));
    subim2 = exp(1i*angle(subim2));
end

%%% Tapering %%%
filterx = fract_hanning_pad(size(subim1,2),size(subim1,2),size(subim1,2)-2*taper);
filterx = fftshift(filterx(1,:));
filtery = fract_hanning_pad(size(subim1,1),size(subim1,1),size(subim1,1)-2*taper);
filtery = fftshift(filtery(:,1));
filterxy = filterx.*filtery;

% Taper subimages %       
subim1 = subim1.*filterxy;% + (1-filterxy).*mean(subim1(:));
subim2 = subim2.*filterxy;% + (1-filterxy).*mean(subim2(:));

if param.plotting > 1
    plotting.smart_figure(23)
    set(gcf,'Outerposition',[1 1 500 476])    %[left, bottom, width, height
    if strcmpi(param.image_prop,'phase') || flag_imread
        imagesc(subim1);
    else
        imagesc(angle(subim1));
    end
    axis xy equal tight
    colormap bone; colorbar
    if param.prop_obj
        [si_unit, val] = utils.get_unit_length(param.prop_obj);
        title(sprintf([param.fname{1} '\npropagated by %d %s'], val, si_unit),'interpreter','none')
    else
        title(param.fname{1},'interpreter','none')
    end
    plotting.smart_figure(24)
    if strcmpi(param.image_prop,'phase') || flag_imread
        imagesc(real(subim2));
    else
        imagesc(angle(subim2));
    end
    axis xy equal tight
    colormap bone; colorbar
    if param.prop_obj
        [si_unit, val] = utils.get_unit_length(param.prop_obj);
        title(sprintf([param.fname{2} '\npropagated by %d %s'], val, si_unit),'interpreter','none')
    else
        title(param.fname{2},'interpreter','none')
    end
    set(gcf,'Outerposition',[500 1 500 476])    %[left, bottom, width, height
end
%% Computing the FSC
param.st_title = sprintf('%s\n %s\n flipped_images %d, taper %d',param.fname{1}, param.fname{2}, param.flipped_images, taper);
if flag_imread
    subim1 = real(subim1);
    subim2 = real(subim2);
    warning('Assuming images are in real number!');
end
    
subim1 = utils.stabilize_phase(subim1, subim2, 'binning', 4, 'fourier_guess', false); 
[resolution,FSC,T,freq,n,stat] = fourier_shell_corr_3D_2(subim1,subim2, param);


%% visually compare alignment quality 
if param.plotting>2
    plotting.smart_figure(4545)
    subplot(1,2,1)
    imagesc3D(angle(subim1 .* conj( subim2)))
    colormap bone
    axis off image
    colorbar
    title('Phase difference between aligned sub-images')

    subplot(1,2,2)
    imagesc3D(cat(3,angle(subim1),angle(subim2)))
    colormap bone
    axis off image
    title('Compare aligned sub-images')
    colorbar
    plotting.suptitle('Visually compare quality and verify alignement / drifts')
end

end


function img = apply_apod(img, asize)
ob_good_range = {asize(1)/2:size(img,1)-asize(1)/2, asize(2)/2:size(img,2)-asize(2)/2};
filt_size = [size(ob_good_range{1},2) size(ob_good_range{2},2)];
img = img.*fftshift(utils.filt2d_pad(size(img), max(1,filt_size), max(1,filt_size-min(floor(filt_size.*0.05)))));
end