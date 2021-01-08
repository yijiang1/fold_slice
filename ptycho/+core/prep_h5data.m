%PREP_H5DATA prepare data and save it to disk
%   prep_h5data expects that fmask and fmag already exist, prepares them
%   for the C++ code and saves everything to disk. 
%
% ** p          p structure
%
% see also: core.ptycho_prepare_scans


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

function prep_h5data(p)
import utils.verbose
import io.HDF.save2hdf5


fmask = p.fmask;
fmag = p.fmag;

single_mask = true;
for ii=1:size(fmask,3)-1
    if ~isequaln(fmask(:,:,ii),fmask(:,:,ii+1))
        single_mask = false;
        break;
    end
end

%%%%% Prepare object and probe for hdf5 file %%%%%
for obnum = p.share_object_ID
    object_c(obnum) = struct('data',p.object{obnum}(:,:,1));
end
for prnum = p.share_probe_ID
    probe_c(prnum) = struct('data',p.probes(:,:,prnum));
end


%%%%% Prepare data for hdf5 file %%%%%
    bad_pixels = cell(p.numscans,1); 
    bad_pixels_index = cell(p.numscans,1);
    % Assumes that detector position is the same within a scan
    for ii = 1:p.numscans % loop over scans
        
        
        % Prepare structure for detector %
        % Here I detect the module gaps, can be later given by the
        % prepare_data function
        fmaski = fmask(:,:,p.scanindexrange(ii,1)); % First mask to detect modules (gaps)
        auxmodxo = any(fmaski,1);
        if auxmodxo(1)
            indbeginmodx = 1;
        else
            indbeginmodx = [];
        end
        auxmodx = diff(auxmodxo);
        indbeginmodx = [indbeginmodx find(auxmodx==1)+1];
        nummodx = length(indbeginmodx);
        indendmodx = find(auxmodx==-1);
        if auxmodxo(end)
            indendmodx = [indendmodx p.asize(2)];
        end
        
        auxmodyo = any(fmaski,2);
        if auxmodyo(1)
            indbeginmody = 1;
        else
            indbeginmody = [];
        end
        auxmody = diff(auxmodyo);
        indbeginmody = [indbeginmody find(auxmody==1).'+1];
        nummody = length(indbeginmody);
        indendmody = find(auxmody==-1).';
        if auxmodyo(end)
            indendmody = [indendmody p.asize(1)];
        end
        
        modulearray = zeros(nummody*nummodx,4);
        fmaskdet = zeros(p.asize); % Module mask for current detector position
        counter = 0;
        for kk = 1:nummody
            for jj = 1:nummodx
                counter = counter+1;
                numrows = indendmody(kk) - indbeginmody(kk) + 1;
                numcols = indendmodx(jj) - indbeginmodx(jj) + 1;
                fmaskdet(indbeginmody(kk):indendmody(kk),indbeginmodx(jj):indendmodx(jj))=1;
                modulearray(counter,:) = [numrows,numcols,indbeginmody(kk)-1,indbeginmodx(jj)-1]; %% Minus one to go to indexing convention in C
            end
        end
        
        verbose(3,['Detected ' num2str(nummodx*nummody) ' modules'])
        if verbose >= 3
            disp([num2str(modulearray)])
        end
        
        %%% Here there is the posibility to add bad pixels that are common
        %%% to all diffraction patterns. Could be identified in prepare
        %%% data
        %detector(ii) = struct('rows', uint32(192), 'columns', uint32(192), 'modules', transpose(uint32([192,192,0,0])), 'bad_pixels', transpose(uint32([9,10; 11,12; 13,14])));
        detector(ii) = struct('rows', uint32(p.asize(1)), 'columns', uint32(p.asize(2)),...
            'modules', transpose(uint32(modulearray)),'bad_pixels',uint32([]));
        
        if ~single_mask
            idx = 0;
            for jj = p.scanindexrange(ii,1):p.scanindexrange(ii,2) % loop over corresponding diffraction patterns
                [y, x] = find(1+fmask(:,:,jj)-fmaskdet==0);
                bps = transpose(reshape([y, x], length(x), 2))-1;
                idx = length(bps)+idx;
                bad_pixels_index{ii} = [bad_pixels_index{ii} idx];
                bad_pixels{ii} = [bad_pixels{ii} bps];
            end
        else
            [y, x] = find(1+fmaski-fmaskdet==0);
            bad_pixels{ii} = transpose(reshape([y, x], length(x), 2))-1;
        end
        
            % Prepare structure for measurement %
%         for jj = p.scanindexrange(ii,1):p.scanindexrange(ii,2) % loop over corresponding diffraction patterns
%             %%%%% Prepare bad pixels %%%%%
%             [y, x] = find(1+fmaski-fmaskdet==0);
%             bad_pixels(jj) = transpose(reshape([y, x], length(x), 2))-1;
%             measurement(jj) = struct('data', fmag(:,:,jj), 'position', uint32((p.positions(jj,:))),...
%                 'object', uint32(p.share_object_ID(ii)-1), 'probe', uint32(p.share_probe_ID(ii)-1),...
%                 'detector', uint32(ii-1));
%         end
    end
    
    
    %% prepare output
    
    h5_struc = [];
    h5_struc.Attributes.format = 2;
    for ii=1:size(probe_c,2)
        h5_struc.probes(:,ii) = uint64(size(probe_c(ii).data));
    end
    for ii=1:size(object_c,2)
        h5_struc.objects(:,ii) = uint64(size(object_c(ii).data));
    end
    
    
    %% detectors
    h5_struc.detector = [];
    for ii=1:p.numscans
        temp = detector(ii);
        h5_struc.detector.(['n' num2str(ii-1)]).Attributes.rows = temp.rows;
        h5_struc.detector.(['n' num2str(ii-1)]).Attributes.columns = temp.columns;
        if single_mask && ~isempty(bad_pixels{ii})
            h5_struc.detector.(['n' num2str(ii-1)]).bad_pixels = bad_pixels{ii};
        end
        h5_struc.detector.(['n' num2str(ii-1)]).modules = temp.modules;
        
    end
    
    %% measurements
    h5_struc.measurement = [];
    h5_struc.measurement.Attributes.max_power = 1/p.renorm^2;
    
    for ii=1:p.numscans
        % attributes
        h5_struc.measurement.(['n' num2str(ii-1)]).Attributes.detector = uint32(ii-1);
        h5_struc.measurement.(['n' num2str(ii-1)]).Attributes.probe = uint32(p.share_probe_ID(ii)-1);
        h5_struc.measurement.(['n' num2str(ii-1)]).Attributes.object = uint32(p.share_object_ID(ii)-1);
        h5_struc.measurement.(['n' num2str(ii-1)]).Attributes.max_sum = uint32(p.max_sum(ii));
        
        % datasets
        h5_struc.measurement.(['n' num2str(ii-1)]).positions = uint32(transpose(round(p.positions(p.scanidxs{ii},:))));
        h5_struc.measurement.(['n' num2str(ii-1)]).data = permute((fmag(:,:,p.scanidxs{ii})/p.renorm).^2, [2 1 3]);
        if ~single_mask
            h5_struc.measurement.(['n' num2str(ii-1)]).bad_pixels = bad_pixels{ii};
            h5_struc.measurement.(['n' num2str(ii-1)]).bad_pixels_index.Value = uint64(transpose(bad_pixels_index{ii}));
            h5_struc.measurement.(['n' num2str(ii-1)]).bad_pixels_index.Attributes.save2hdf5DataShape = size(uint64(transpose(bad_pixels_index{ii})),1);
        end

    end
    
    %% save to disk
    if ~exist(p.prepare_data_path, 'dir')
        mkdir(p.prepare_data_path)
    end
    verbose(2,'Writing H5 data file: %s',[p.prepare_data_path p.prepare_data_filename]);
    save2hdf5([p.prepare_data_path p.prepare_data_filename], h5_struc, 'overwrite', true, 'comp', p.io.data_compression);

end

