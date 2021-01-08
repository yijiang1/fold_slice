function [img1_aligned, img2_aligned, delta_all] = fun_align_img(img1,img2,asize,pix)

%   [resolution] = aligned_FSC(file1,file2,params)
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
%
% Inputs:
%
% file1     Filename with path of reconstruction 1
% file2     Filename with path of reconstruction 2
% params    Structure with parameters as described below
%
%  params.flipped_images    Flip one input image horizontally (= true or false).
%                           Useful when comparing 0 and 180 degree projections 
%                           in tomography (default = false).
%  params.crop      = '';       for using the default half size of the probe
%                   = 'manual'  for using GUI to select region. This will display the range, e.g. {600:800, 600:800}     
%                   = {600:800, 600:800} for custom vertical and horizontal cropping, respectively
%  params.GUIguess  To click for an initial alignment guess, if used it ignores 
%                   the values of params.guessx and params.guessy (default
%                   = false)
%  params.guessx
%  params.guessy        An intial guess for x and y alignment (default = [])
%  params.remove_ramp   Try to remove linear phase from whole image before initial 
%                       alignment (default = true)
%  params.image_prop    = 'complex' 
%                       = 'phasor' (phase with unit amplitude, default) 
%                       = 'phase'  (Note: phase should not be used if there is phase wrapping)
%  params.taper         = 20 (default)   Pixels to taper images - Increase until the FSC does not change anymore
%  params.plotting      Display plots (default = false)
%  params.dispfsc       Display FSC plot (default = true)
%  params.SNRt          SNR for FSC threshold curve   
%                       SNRt = 0.2071 for 1/2 bit threshold for resolution of the average of the 2 images
%                       SNRt = 0.5    for 1   bit threshold for resolution of each individual image (default)
%  params.thickring     Thickness of Fourier domain ring for FSC in pixels (default = 1)
%  params.freq_thr      To ignore the crossings before freq_thr for determining resolution (default 0.02)
%  params.out_fn        Filename of output of jpeg for FSC

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|      Authors: CXS group
%|                                                                       |
%*-----------------------------------------------------------------------*
%
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form and additionally you should cite the references relevant
%   to this code.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checks and defaults %%%
flag_imread = 1;

GUIguess = false;
plotting = false;
crop = '';
guessx = [];
guessy = [];
remove_ramp = false;

image_prop = 'complex';
taper = 20;
dispfsc = true;
SNRt = 0.5;
thickring = 5;
freq_thr = 0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
img1_orig = img1;
img2_orig = img2;


screensize = get( 0, 'Screensize' ); 
% Show phase images (not cropped)%
if plotting
    figure(2)
    set(gcf,'Outerposition',[1 screensize(4)-550 500 500])    %[left, bottom, width, height
    if ~flag_imread
        imagesc(angle(img1));
    else
        imagesc(img1);
    end
    axis xy equal tight
    colormap bone
    if ~flag_imread 
        aux = angle(img1); %angle(img1(round(asize(1)/2):end-round(asize(1)/2),round(asize(2)/2):end-round(asize(2)/2)));
    else
        aux = img1;
    end
    caxis([min(aux(:)) max(aux(:))]); colorbar
    title(file{1},'interpreter','none')
    figure(3)
    if ~flag_imread
        imagesc(angle(img2));
    else
        imagesc(img2);
    end
    axis xy equal tight
    colormap bone
    if ~flag_imread 
        aux = angle(img2); %angle(img2(round(asize(1)/2):end-round(asize(1)/2),round(asize(2)/2):end-round(asize(2)/2)));
    else
        aux = img2;  
    end
    caxis([min(aux(:)) max(aux(:))]); colorbar
    title(file{2},'interpreter','none')
    set(gcf,'Outerposition',[500 screensize(4)-550 500 500])    %[left, bottom, width, height
end

% Crop images - default is half the size of the probe on each side plus
% whatever needed to make them of equal size
minsize1 = min(size(img1,1),size(img2,1));
minsize2 = min(size(img1,2),size(img2,2));
if isempty(crop)
    crop = {round(asize(1)/2):minsize1-round(asize(1)/2), ...
        round(asize(2)/2):minsize2-round(asize(1)/2)};
elseif strcmpi(crop, 'manual')
    figure()
    imagesc(angle(img1))
    colormap bone 
    axis image xy 
    title('Select compared region')
    rect = round(getrect);
    crop = {rect(2)+(1:rect(4)),rect(1)+(1:rect(3))};
    disp('===========================')
    fprintf('Selected region: {%i:%i,%i:%i}\n',rect(2), rect(2)+rect(4), rect(1), rect(1)+rect(3));
    disp('===========================')
    pause(1)
end
img1 = img1(crop{:});
img2 = img2(crop{:});


if GUIguess
    figure(2)
    disp(['Click on a feature on figure 2'])
    [xin yin] = ginput(1);
    figure(3)
    disp(['Click on a feature on figure 3'])
    [xin2 yin2] = ginput(1);
    guessx = round(xin-xin2);
    guessy = round(yin-yin2);
end

if ~isempty(guessx)
    switch sign(guessx)
        case 1
            img1 = img1(:,1+guessx:end);
            img2 = img2(:,1:end-guessx);
        case -1
            img1 = img1(:,1:end+guessx);
            img2 = img2(:,1-guessx:end);
    end
end
if ~isempty(guessy)
    switch sign(guessy)
        case 1
            img1 = img1(1+guessy:end,:);
            img2 = img2(1:end-guessy,:);
        case -1
            img1 = img1(1:end+guessy,:);
            img2 = img2(1-guessy:end,:);
    end
end

% Remove ramp 
if remove_ramp
    disp('Removing ramp for initial alignment')
    img1 = remove_linearphase_v2(img1,ones(size(img1)),100);
    img2 = remove_linearphase_v2(img2,ones(size(img2)),100);
end
    
if plotting
    figure(4)
    set(gcf,'Outerposition',[1 1 500 476])    %[left, bottom, width, height
    if ~flag_imread
        imagesc(angle(img1));
    else
        imagesc(img1); 
    end
    axis xy equal tight
    colormap bone; colorbar
    title(file{1},'interpreter','none')
    figure(5)
    if ~flag_imread
        imagesc(angle(img2));
    else
        imagesc(img2);
    end
    axis xy equal tight
    colormap bone; colorbar
    title(file{2},'interpreter','none')
    set(gcf,'Outerposition',[500 1 500 476])    %[left, bottom, width, height
end

%%% Initial alignment %%%
fprintf('\nInitial alignment\n')  
if ~flag_imread
    switch lower(image_prop)
        case  'complex'
            imgalign1 = img1;
            imgalign2 = img2;
            disp('Registering complex valued images')
        case 'phasor'
            imgalign1 = ones(size(img1)).*exp(1i*angle(img1));
            imgalign2 = ones(size(img1)).*exp(1i*angle(img2));
            disp('Registering phasor of complex valued images')
        case 'phase'
            imgalign1 = angle(img1);
            imgalign2 = angle(img2);
            disp('Registering phase of complex valued images')
    end
else
    imgalign1 = img1;
    imgalign2 = img2; 
end

upsamp = 100;
displ = 1;
W = 1;
x1 = [];%[1:150];
x2 = x1;
y1 = [];%[1:238];
y2 = y1;
% imgalign2 = shiftpp2(imgalign2,10,-10); % To test range adjustment
[subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(imgalign1,imgalign2, x1, y1, x2, y2, upsamp, displ,1);

%%% Fine alignment (second round) %%%

% Remove ramp for fine alignment
disp('Removing ramp for fine alignment')
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

if ~flag_imread
    switch lower(image_prop)
        case  'complex'
            subimgalign1 = subimg1;
            subimgalign2 = subimg2;
            disp('Registering complex valued images')
        case 'phasor'
            subimgalign1 = ones(size(subimg1)).*exp(1i*angle(subimg1));
            subimgalign2 = ones(size(subimg1)).*exp(1i*angle(subimg2));
            disp('Registering phasor of complex valued images')
        case 'phase'
            subimgalign1 = angle(subimg1);
            subimgalign2 = angle(subimg2);
            disp('Registering phase of complex valued images')
    end
else
   subimgalign1 = subimg1;
   subimgalign2 = subimg2; 
end

% Fine alignment %
disp(sprintf('\nFine alignment'))
[subim1, subim2, delta2, deltafine2, regionsout] = registersubimages_2(subimgalign1,subimgalign2, x1, y1, x2, y2, upsamp, displ,1);

%%% Tapering %%%
filterx = fract_hanning_pad(size(subim1,2),size(subim1,2),size(subim1,2)-2*taper);
filterx = fftshift(filterx(1,:));
filterx = repmat(filterx,[size(subim1,1) 1]);
filtery = fract_hanning_pad(size(subim1,1),size(subim1,1),size(subim1,1)-2*taper);
filtery = fftshift(filtery(:,1));
filtery = repmat(filtery,[1 size(subim1,2)]);
filterxy = filterx.*filtery;

% Taper subimages %
subim1 = subim1.*filterxy;% + (1-filterxy).*mean(subim1(:));
subim2 = subim2.*filterxy;% + (1-filterxy).*mean(subim2(:));

if plotting
    figure(4)
    set(gcf,'Outerposition',[1 1 500 476])    %[left, bottom, width, height
    if strcmpi(image_prop,'phase') || flag_imread
        imagesc(subim1);
    else
        imagesc(angle(subim1));
    end
    axis xy equal tight
    colormap bone; colorbar
    title(file{1},'interpreter','none')
    figure(5)
    if strcmpi(image_prop,'phase') || flag_imread
        imagesc(real(subim2));
    else
        imagesc(angle(subim2));
    end
    axis xy equal tight
    colormap bone; colorbar
    title(file{2},'interpreter','none')
    set(gcf,'Outerposition',[500 1 500 476])    %[left, bottom, width, height
end
%% Computing the FSC
param.st_title = sprintf('taper %d',taper);
param.pixel_size = pix; 

[resolution FSC T freq] = fourier_shell_corr_3D_2(subim1,subim2, param);

if 0
    img1_aligned = img1_orig;
    img1_aligned(asize(1)/2:asize(1)/2+size(subim1,1)-1, asize(2)/2:asize(2)/2+size(subim1,2)-1) = subim1;

    img2_aligned = img2_orig;
    img2_aligned(asize(1)/2:asize(1)/2+size(subim2,1)-1, asize(2)/2:asize(2)/2+size(subim2,2)-1) = subim2;
else
    img1_aligned = subim1;
    img2_aligned = subim2;
end

delta_all = round(delta) + delta2;

return