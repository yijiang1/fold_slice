%PREPARE_DATA_MATLAB prepares and normalizes data for matlab engines

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

function [ p ] = prep_data_matlab( p, data, fmask )
import utils.pshift
import utils.verbose
import utils.get_option
import math.fftshift_2D

scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number

% Region of interest [xmin xmax ymin ymax]
% Warning Not implemented for shared scans
% Need to update p.numpts, p.positions, data, fmask, indices
if isfield(p, 'scan') && isfield(p.scan, 'roi') &&  ~isempty(p.scan.roi)
    % Convert to p.positions centered on the object
    if p.share_object
       positions_centered = p.positions - p.object_size(1,:)/2 + p.asize(:)/2;
    else
        for ii = 1:p.numscans
            positions_centered(p.scanidxs{ii},:) = p.positions(p.scanidxs{ii},:) - p.object_size(ii,:)/2;
        end
    end
    xmin = p.scan.roi(1)/p.dx_spec(2);
    xmax = p.scan.roi(2)/p.dx_spec(2);
    ymin = p.scan.roi(3)/p.dx_spec(1);
    ymax = p.scan.roi(4)/p.dx_spec(1);
    % Quick check
    if (xmin>xmax)
        error('ROI is empty, xmax<xmin')
    elseif (ymin>ymax)
        error('ROI is empty, ymax<ymin')
    end
    % Do the comparison and update accordingly [xmin xmax ymin ymax]
    whichtokeep = find(  (positions_centered(:,1) > ymin) & ...
        (positions_centered(:,1)< ymax) & ...
        (positions_centered(:,2) > xmin) & ...
        (positions_centered(:,2) < xmax) );
    % update p.numpts, positions, data, fmask, indices
    data = data(:,:,whichtokeep);
    fmask = fmask(:,:,whichtokeep);
    p.positions = p.positions(whichtokeep,:);
    
    for ii = 1:p.numscans
        %         indaux = find(whichtokeep<scanfirstindex(ii+1),'last'); % Finds the last index of scan ii in the new variables
        indaux = find((whichtokeep>=scanfirstindex(ii))&(whichtokeep<scanfirstindex(ii+1))); % Finds indices of scan ii in the new variables
        if isempty(indaux)
            error('ROI: Scan %s has no points in the ROI. Change the ROI or remove this scan from the scan list',p.scan_str{ii})
        end
        p.numpts(ii) = numel(indaux);
    end
    scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
    for ii = 1:p.numscans
        p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
    end
    verbose(2, 'Region of interest reduced number of points from %d to %d', sum(p.numpts_orig), sum(p.numpts));
    
    % Redo a convenient offset
    p.positions = p.positions - min(p.positions) + 1;
    
   
    % Recompute object sizes
    p = update_object_size(p); 

   
end



% Skip some data point (for testing reduced dose) 
if isfield(p, 'skip_points')
    if ~isempty(p.skip_points) && p.skip_points>1
        offset = mod(p.scan_number, p.skip_points);
        %offset = mod([0 1 2], p.skip_points);
        whichtokeep = [1+offset(1):p.skip_points:p.numpts(1)];
        for idx = 2:length(p.numpts)
            last_idx = sum(p.numpts(1:idx-1));
            whichtokeep = [whichtokeep , (last_idx+1+offset(idx)):p.skip_points:(last_idx+p.numpts(idx))];
        end
        data = data(:,:,whichtokeep);
        fmask = fmask(:,:,whichtokeep);
        p.positions = p.positions(whichtokeep,:);
        for ii = 1:p.numscans
            %         indaux = find(whichtokeep<scanfirstindex(ii+1),'last'); % Finds the last index of scan ii in the new variables
            indaux = find((whichtokeep>=scanfirstindex(ii))&(whichtokeep<scanfirstindex(ii+1))); % Finds indices of scan ii in the new variables
            if isempty(indaux)
                error('ROI: Scan %s has no points in the ROI. Change the ROI or remove this scan from the scan list',p.scan_str{ii})
            end
            p.numpts(ii) = numel(indaux);
        end
        scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
        for ii = 1:p.numscans
            p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
        end
        verbose(2, 'Region of interest reduced number of points from %d to %d', sum(p.numpts_orig), sum(p.numpts));

        % Recompute object sizes
        if p.share_object
            p.object_size = p.asize + max(p.positions,[],1);
            verbose(3, 'Computed object size: %d x %d', p.object_size(1), p.object_size(2));
        else
            for ii = 1:p.numscans
                p.object_size(ii,:) = p.asize + max(p.positions(p.scanindexrange(ii,1):p.scanindexrange(ii,2),:),[],1);
                verbose(3, 'Computed object size: %d x %d', p.object_size(ii,1), p.object_size(ii,2));
            end
        end
    end
end

% % padding
% if any(datasize ~= p.asize)
%     newdata = zeros([p.asize, num_difpat]);
%     offset = floor(.5*(p.asize-datasize));
%     newdata(offset(1) + (1:datasize(1)), offset(2) + (1:datasize(2)), :) = data;
%     data = newdata;
%     clear newdata
%     newfmask = ones([p.asize, size(fmask,3)]);
%     newfmask(offset(1) + (1:datasize(1)), offset(2) + (1:datasize(2)), :) = fmask;
%     fmask = newfmask;
%     clear newfmask
% end


% 'auto_center_data' option is useful if detector shifts provided in
% template cannot be trusted, ie for more than 2 joined scans 
if get_option(p,'auto_center_data') && ~get_option(p,'get_artificial_data')
    for i = 1:p.numscans
        ind = p.scanidxs{i};
        [x,y] = math.center(mean(data(:,:,ind) .* fmask(:,:,ind),3));
        verbose(2,'Auto-shifting diffraction patterns by %i %i px', round(x), round(y))
        data(:,:,ind)  = utils.imshift_fast(data(:,:,ind), x, y, [], 'nearest');
        if size(fmask,3) == size(data,3)
            fmask(:,:,ind) = utils.imshift_fast(fmask(:,:,ind), x, y, [], 'nearest');
        elseif p.numscans == 1
            fmask = utils.imshift_fast(fmask, x, y, [], 'nearest');
        else
           error('Not implemented mask shifting option') 
        end
    end
end
    

p.fmask_per_scan = ndims(fmask) == 3;




% Prepare Fourier projections

% normalization ignores valid mask - should modify
for ii=1:p.numscans
    
    p.max_sum(ii) = max(sum(sum(data(:,:,p.scanidxs{ii}),1),2),[],3);
end
max_power = max(p.max_sum) / prod(p.asize);
p.renorm = sqrt(1/max_power);
p.Nphot = sum(data(:).*fmask(:)); % Number of photons for regularization normalization (ML optimization)

% store mask pre-fftshifted
p.fmask = fftshift_2D(fmask);
clear fmask 
% precalculate modulus of data, normalize and fftshift
p.fmag = fftshift_2D(sqrt(data)) * p.renorm;


end

