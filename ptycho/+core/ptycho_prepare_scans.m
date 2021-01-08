%PTYCHO_PREPARE_SCANS Load data from disk and prepare the scan for the
%ptychographic reconstruction.
%
% ** p          p structure
%
% returns:
% ++ p          p structure
% ++ status     status flag
%
%
% see also: core.initialize_ptycho
%
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

function [ p, status] = ptycho_prepare_scans( p )
import utils.verbose
import utils.get_option

status=1;

% legacy
if ~check_option(p.prepare, 'legacy')
    p.prepare.legacy = false;
end


%% check for lock file
for ii = 1:length(p.scan_number)
    % Write lock file
    if p.queue.lockfile
        lock_filename = [p.save_path{ii} '/' p.run_name '_lock'];
        if exist(lock_filename, 'file')
            verbose(1,sprintf('%s locked by other instance of this script. Continue with next scan.', p.scan_str{ii}));
            out = [];
            status = 0;
            p.getReport.completed = true;
            return
        else
            if ~exist(p.save_path{ii},'dir')
                mkdir(p.save_path{ii});
            end
            verbose(2, 'Creating lock file: %s', lock_filename);
            unix(['touch ' lock_filename]);
        end
    end
end

    
%% prepare data

% check prepared data file
if ~p.prepare.force_preparation_data
    if ~exist(fullfile(p.prepare_data_path, p.prepare_data_filename), 'file')
        verbose(1,'Missing prepared data %s, Forcing data preparation',fullfile(p.prepare_data_path, p.prepare_data_filename) )
        p.prepare.force_preparation_data = true;
    else
        try
            if core.check_prepared_data(p)
                verbose(2, 'Prepared data does not match reconstruction parameters. Forcing data preparation.')
                p.prepare.force_preparation_data = true;
            end
        catch ME
            verbose(2, [ME.getReport '\n']);
            verbose(2, 'Prepared data check failed. Forcing data preparation.')
            p.prepare.force_preparation_data = true;
        end
    end
end

% check if we need to prepare the data
if p.prepare.auto_prepare_data && (~exist(fullfile(p.prepare_data_path, p.prepare_data_filename),'file')||p.prepare.force_preparation_data)
    prepare_data_bool = true;
else
    prepare_data_bool = false;
end

% prepare initial object and probes
p = core.prepare_initial_guess(p);

if prepare_data_bool && ~p.prepare.legacy
    verbose(2, 'Loading raw data')
    p = core.run_data_preparator(p);
elseif p.prepare.legacy
    %% prepare.legacy mode to load prepared data from mat files
    % Prepare data path
    error('Loading from .mat files is not supported anymore.')
end

%% convert function handles to strings
for jj=1:length(p.detectors)
    fn = fieldnames(p.detectors(jj).params);
    for ii=1:length(fn)
        if isa(p.detectors(jj).params.(fn{ii}) , 'function_handle')
            p.detectors(jj).params.(fn{ii}) = func2str(p.detectors(jj).params.(fn{ii}));
        end
    end
    
    funcs_nm = fieldnames(p.detectors(jj).funcs);
    for ii=1:length(funcs_nm)
        p.detectors(jj).funcs.(funcs_nm{ii}) = func2str(p.detectors(jj).funcs.(funcs_nm{ii}));
    end
    
end

% cleanup detector storage
if isfield(p.detectors, 'detStorage')
    p.detectors = rmfield(p.detectors, 'detStorage');
end

p.detectors = struct2cell(p.detectors);



%% output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% save/load prepared data to/from disk %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.prepare.store_prepared_data && (prepare_data_bool || p.prepare.legacy) && ~strcmpi(p.prepare.data_preparator, 'libDetXR')
    core.prep_h5data(p);
end

% check if data needs to be loaded from disk
if  (~isfield(p,'fmag') || isempty(p.fmag) )&& (p.plot.prepared_data || ~p.external_engine0 || check_option(p, 'force_prepare_h5_files')) 
  
%% TO BE DELETED 
%   if (~p.external_engine0 && (~isfield(p,'fmag')) ) || p.external_engine0 &&  ...
%     (isfield(p.engines{1}, 'force_prepare_h5_files') && p.engines{1}.force_prepare_h5_files) || ...
%     (prepare_data_bool && strcmpi(p.prepare.data_preparator,'libDetXR') && p.plot.prepared_data && ~p.prepare.legacy)|| ...
%     (~prepare_data_bool && p.plot.prepared_data) || (p.model_object && strcmpi(p.model_object_type, 'prep_data') && ~prepare_data_bool)

    verbose(2, 'Loading already prepared data.');
    [p.fmag, p.fmask, pos, max_power, scanindexrange, p.max_sum] = io.load_prepared_data(fullfile(p.prepare_data_path ,p.prepare_data_filename));
    
    p.renorm = sqrt(1/max_power);
    p.Nphot = sum((p.fmag(:)/p.renorm).^2.*p.fmask(:));
    p.fmask_per_scan = (length(size(p.fmask)) == 3);   
    
    % shall the positions be overwritten?
    if p.io.load_prep_pos
        verbose(2,'Overwriting positions with values from prepared data.')
        p.positions = pos;
        p.scanindexrange = scanindexrange;
        p.numpts = diff(reshape(scanindexrange', 2, []),1);
        p.numpts(1) = p.numpts(1)+1; % scanindexrange starts at 1
        for ii = 1:p.numscans
            p.scanidxs{ii} = scanindexrange(ii,1):scanindexrange(ii,2);
        end
        
    end
    
end

% if not yet done, 
% apply binning on all relevant 
% parameters except data and mask that are already done 
if (p.detector.binning || p.detector.upsampling) && any(p.asize ~= (p.asize_nobin .* 2^-p.detector.binning * 2^p.detector.upsampling))
    if p.detector.binning
        p = core.apply_binning(p, 2^p.detector.binning); 
    end
    if p.detector.upsampling
        % just reverse operation to the binning 
        p = core.apply_binning(p, 2^(-p.detector.upsampling) ); 
    end
end


% initial guess for Fourier ptychography 
if p.model_object && strcmpi(p.model.object_type, 'amplitude')
    k = 2*pi/p.lambda;
    objpix = p.lambda*p.z_lens./(p.asize.*p.dx_spec);
    for ii=1:length(p.object)
        [Xp,Yp] = utils.get_grid(p.object_size(ii,:), objpix(1));
        pre_phase_factor = exp(-1i*k*((Xp./(p.object_size(ii,2)/(p.asize(2)))).^2+(Yp./(p.object_size(ii,1)/(p.asize(1)))).^2)/(2*p.z_lens));
        init_guess = rot90(mean(p.fmag,3),2);
        init_guess = init_guess./(max(max(init_guess)))./(p.asize(1).*p.asize(2)).*p.numpts;
        
        p.object{ii} = fftshift(fft2(ifft2(ifftshift(utils.crop_pad(fftshift(fft2(init_guess.*exp(-1j.*(init_guess./(max(max(init_guess))).*2*pi-pi)))), p.object_size(ii,:)))).*fftshift(pre_phase_factor)));

    end
end



end

