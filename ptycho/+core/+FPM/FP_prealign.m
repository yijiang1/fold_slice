%% FP_PREALIGN prealign Fourier ptychographic data
% [p] = FP_prealign(p)
% FP_prealign is an alignment routine for Fourier ptychographic
% measurements. Due to the changes in frequency content, a global
% registration is not sufficient. We therefore sort the frames such that
% images with similar frequency content can be aligned.
%
% ** p                  p structure    
% 
% *taken from p.prealign:*
% ** ctr_sh             shift the center before cropping the final dataset to asize
% ** crop_dft           crop the images by crop_dft before calculating the dftregistration
% ** axis               start the alignment procedure along specified axis or rotation (1 or 2)
% ** numiter            number of iterations
% ** rad_filt_min       discard positions below rad_filt_min
% ** rad_filt_max       discard positions beyond rad_filt_max
% ** mfiles             discard specific data points
% ** flat_corr          apply a flat-field correction
% ** filt_align         remove interpolation artifacts
% ** save_alignment     save final shifts / alignment
% ** load_alignment     overwrite alignment with previous alignment
% ** alignment_file     specify path+file 
% ** plot alignment     turn on/off plotting during the alignment
%
% returns:
% ++ p                  p structure
%
% see also: <base> utils.dftregistration
% 


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

function [p] = FP_prealign(p)

import utils.*

%% initial checks

par = p.prealign;

% parse inputs
check_input = @(x) islogical(x) || isnumeric(x);
parse_par = inputParser;
parse_par.KeepUnmatched = true;

parse_par.addParameter('ctr_sh', [0 0], @isnumeric)
parse_par.addParameter('crop_dft', 1, @isnumeric)
parse_par.addParameter('axis', 1, @isnumeric)
parse_par.addParameter('numiter', 3, @isnumeric)
parse_par.addParameter('rad_filt_min', 0, @isnumeric)
parse_par.addParameter('rad_filt_max',Inf, @isnumeric)
parse_par.addParameter('mfiles', [], @isnumeric)
parse_par.addParameter('save_alignment', false, check_input)
parse_par.addParameter('load_alignment', false, check_input)
parse_par.addParameter('plot_alignment', false, check_input)
parse_par.addParameter('sort_pos_radii', false, check_input)
parse_par.addParameter('mag_est', [], @isnumeric)

parse_par.parse(par);
par = utils.update_param(par, parse_par.Results);

if isempty(par.mag_est)
    warning('The magnification was not specified. I will assume mag_est=100.')
    par.mag_est = 100;
end

%% load data and optimize position arangement 

detStorage = p.detectors(p.scanID).detStorage;
detParams = p.detectors(p.scanID).params;

if p.scanID==1
    if ~isfield(detParams, 'detposmotor')
        error('Please specify the detector motor (det.detposmotor) in your detector template.')
    end
    p.positions_real = p.positions_orig;
    p.det_pos = p.positions_real(p.scanidxs{p.scanID},:);
    if (~isempty(p.spec.motor.coarse_motors))
        % adjust the positions in case of coarse stage movements
        coarse_pos = [];
        for ii = 1:length(p.scan_number)
            coarse_pos = [coarse_pos; [p.meta{ii}.spec.(p.spec.motor.coarse_motors{2}).*1e-3 p.meta{ii}.spec.(p.spec.motor.coarse_motors{1}).*1e-3]];
        end
        coarse_cen = [(max(coarse_pos(:,1)) - min(coarse_pos(:,1)))./2 (max(coarse_pos(:,2))-min(coarse_pos(:,2)))./2];
        for ii = 1:length(p.scan_number)
            p.positions_real(p.scanidxs{ii},1) = p.positions_real(p.scanidxs{ii},1)+p.meta{ii}.spec.(p.coarsey)*1e-3-min(coarse_pos(:,1))-coarse_cen(1);
            p.positions_real(p.scanidxs{ii},2) = p.positions_real(p.scanidxs{ii},2)+p.meta{ii}.spec.(p.coarsex)*1e-3-min(coarse_pos(:,2))-coarse_cen(2);
            p.det_pos(p.scanidxs{ii},1) = p.meta{ii}.spec.hy; 
            p.det_pos(p.scanidxs{ii},2) = p.meta{ii}.spec.hx; 
        end
        p.coarsex = []; 
        p.coarsey = [];
    else
        % load the detector positions
        for ii=1:length(p.scan_number)
            p.det_pos(p.scanidxs{ii},1) = p.meta{ii}.spec.(detParams.detposmotor{1});
            p.det_pos(p.scanidxs{ii},2) = p.meta{ii}.spec.(detParams.detposmotor{2});
        end
    end
    
end
par.det_pos = p.det_pos(p.scanidxs{p.scanID},:);
pos = p.positions_real(p.scanidxs{p.scanID},:);


data = detStorage.data;

% remove unwanted files
if ~isempty(par.mfiles)
    pos(par.mfiles,:) = [];
    data(:,:,par.mfiles) = [];
end

% remove files out of range
indx = [];
for ii=1:size(pos,1)
    if sqrt(pos(ii,1)^2 + pos(ii,2)^2)<par.rad_filt_min || sqrt(pos(ii,1)^2 + pos(ii,2)^2)>par.rad_filt_max %|| pos(ii,2)>35e-6
        indx = [indx ii];
    end
end

pos(indx,:) = [];
par.det_pos(indx,:) = [];
data(:,:,indx) = [];
% if p.scanID>1
%     p.positions_real(indx+sum(p.numpts(1:tmp.ii-1)),:) = [];
% end


par.pos = pos;
par.pos_orig = pos;

%%% update p values with new positions %%%
tmp_append_pos = p.positions_real([p.scanidxs{min(p.scanID+1, length(p.numpts)+1):end}],:);
tmp_append_det = p.det_pos([p.scanidxs{min(p.scanID+1, length(p.numpts)+1):end}],:);

p.numpts(p.scanID) = size(par.pos,1);
scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
for ii = 1:p.numscans
    p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
end
for ii = 1:p.numscans     
    p.scanidxs{ii} = p.scanindexrange(ii,1):p.scanindexrange(ii,2); 
end


% adjust positions container

p.positions_real([p.scanidxs{p.scanID:end}],:) = [];
p.positions_real(p.scanidxs{p.scanID},:) = par.pos;
p.positions_real = [p.positions_real; tmp_append_pos];

p.det_pos([p.scanidxs{p.scanID:end}],:) = [];
p.det_pos(p.scanidxs{p.scanID},:) = par.det_pos;
p.det_pos = [p.det_pos; tmp_append_det];




par.data = data;


par.sz = size(par.data);



fft_mask = ones(par.asize);
fft_mask(round(par.asize(1)/4)-5:round(par.asize(1)/4)+5,:) = 0;
fft_mask(:,round(par.asize(2)/4)-5:round(par.asize(2)/4)+5) = 0;
fft_mask(:,round(par.asize(2)*3/4)-5:round(par.asize(2)*3/4)+5) = 0;
fft_mask(round(par.asize(1)*3/4)-5:round(par.asize(1)*3/4)+5,:) = 0;

if ~isempty(detParams.mask_saturated_value)
    par.data = (abs(ifft2(fftshift(fft_mask).*fft2(par.data.*(par.data<detParams.mask_saturated_value)))));
else
    par.data = (abs(ifft2(fftshift(fft_mask).*fft2(par.data))));
end

clear data;

par.orig_data = (abs(ifft2(fftshift(fft_mask).*fft2(par.data))));%par.data;

if par.prealign_data
    par.sum_shift_total = zeros(par.sz(3), 2);
    
    for ii=1:par.numiter*length(par.type)
        par.iterii = ii;
        verbose(2, 'Iteration %d/%d', ii, par.numiter*length(par.type))
        
        % sort positions
        verbose(3, 'Sorting positions.')
        sort_type = par.type{mod(ii+1,length(par.type))+1};
        par = core.FPM.sort_pos(par, sort_type);
        
        verbose(4, 'Aligning data along sorted positions.')
        % align along sorted positions
        par = core.FPM.align_data(par);
        
        fig30 = plotting.smart_figure(30);
        clf;
        set(groot,'CurrentFigure',fig30);
        imagesc(mean(par.data,3));
        colormap(bone(256))
        title(sprintf('Alignment after %d iteration(s)', ii))
        drawnow()
        
        if ~mod(ii,length(par.type))
            if par.axis==1
                par.axis = 2;
            else
                par.axis = 1;
            end
        end
    end

end



if par.sort_pos_radii
    par = core.FPM.sort_pos_radii(par);
end



% use distortion matrix or load alignment from disk

if par.save_alignment
    sum_shift_total = par.sum_shift_total; 
    save(sprintf('alignment_S%05d.mat', p.scan_number(p.scanID)), 'sum_shift_total');
end

if par.use_distortion_corr && isempty(par.distortion_corr)
    par.distortion_corr = core.FPM.distortion_matrix(p, par.sum_shift_total, par.det_pos, par.mag_est);
elseif ischar(par.distortion_corr)
        f = io.load_ptycho_recons(par.distorion_corr);
        par.distortion_corr = f.p.prealign.distortion_corr;
end

if par.use_distortion_corr
    p.prealign.distortion_corr = par.distortion_corr;
    pos = core.FPM.get_positions(par.distortion_corr, par.pos.*1e3);
    sum_shift_total = (pos - par.det_pos)./p.ds./1e3;
elseif par.load_alignment
    if isempty(par.alignment_file) || ~exist(par.alignment_file, 'file')
        error('Could not load specified alignment file.')
    end
    f = load(par.alignment_file);
    sum_shift_total = f.sum_shift_total;
else
    sum_shift_total = par.sum_shift_total;
end


%%%%%%%%%%%%%%%%
% apply shifts %
%%%%%%%%%%%%%%%%

mask = zeros([size(detStorage.mask) p.numpts(p.scanID)]);
utils.verbose(2, 'Applying shifts to image stack and mask.');
for ii=1:p.numpts(p.scanID)
    data(:,:,ii) = ifftshift(utils.crop_pad(abs(utils.shiftpp2(par.orig_data(:,:,ii), sum_shift_total(ii,1), sum_shift_total(ii,2))), p.asize));
    mask(:,:,ii) = abs(utils.shiftpp2(detStorage.mask, round(sum_shift_total(ii,1)), round(sum_shift_total(ii,2))));
end


if utils.verbose > 2
        fig30 = plotting.smart_figure(30);
        clf;
        set(groot,'CurrentFigure',fig30);
        imagesc(fftshift(mean(data,3)));
        colormap(bone(256))
        title('Alignment')
        drawnow()
end

detStorage.data = data;
detStorage.mask = round(mask);

% if par.prealign_data
    p.positions_real(p.scanidxs{p.scanID},1) = p.positions_real(p.scanidxs{p.scanID},1);
    p.positions_real(p.scanidxs{p.scanID},2) = p.positions_real(p.scanidxs{p.scanID},2);
% end


if p.scanID==length(p.numpts)
    p = core.ptycho_adjust_positions(p);
    p = core.prepare_initial_guess(p);
end

end




