%libDetXR prepares and exports data for libDetXR/ptyhon data
%   preparation.
% 
% ** p              p structure
% 
% returns:
% ++ p              updated p structure
% ++ fdb            feedback structure
%
% see also: core.run_data_preparator


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

function [p, fdb] = libDetXR( p )
import utils.verbose
import io.HDF.save2hdf5
import io.mat2json
import beamline.identify_eaccount


verbose(1, 'Preparing data using json/python data preparation.')

if ~isfield(p, 'delete_temp_json_data')
    p.delete_temp_json_data = true;
end


% Parameters for reconstruction with external C code
if ~isfield(p, 'initial_conditions_file') || isempty(p.initial_conditions_file)%&&p.external_reconstruction
    if ~isempty(p.suffix)
        suffix = ['_' p.suffix];
    else
        suffix = '';
    end
    p.initial_conditions_file = [core.generate_scan_name(p) sprintf('_initial_conditions_%03dx%03d%s.h5', p.asize(1), p.asize(2), suffix)];
    verbose(3, 'C-code initial_conditions_file = %s', p.initial_conditions_file);
end


% path for saving temp data for C code
if isfield(p, 'initial_conditions_path') && ~isempty(p.initial_conditions_path)
    if p.prepare_data_path(1) == '~'
        p.initial_conditions_path = ['/sls/X12SA/Data10/' identify_eaccount p.prepare_data_path(9:end)];
    end
else
    p.initial_conditions_path = p.prepare_data_path;    
end

if ~isfield(p, 'initial_probe_temp_file')
    p.initial_probe_temp_file{1} = [];
end
if ~isfield(p, 'initial_object_temp_file')
    p.initial_object_temp_file{1} = [];
end

if isempty(p.initial_probe_temp_file{1})
    for ii=1:p.numscans
        p.initial_probe_temp_file{ii} = [core.generate_scan_name(p) sprintf('_temp_probe_S%05d_%03dx%03d.mat', p.scan_number(ii), p.asize(1), p.asize(2))];
    end
end
if isempty(p.initial_object_temp_file{1})
    for ii=1:p.numscans
        p.initial_object_temp_file{ii} = [core.generate_scan_name(p) sprintf('_temp_object_S%05d_%03dx%03d.mat', p.scan_number(ii), p.asize(1), p.asize(2))];
    end
end


s=struct;
s.glob          = struct;
s.detector              = struct;
s.measurement = struct;
s.initialCondition = struct;
s.initialCondition.probe    = struct;
s.initialCondition.object   = struct;
s.initialCondition.param    = struct;

for ii = 1:p.numscans
    p.scanID = ii;
    cid = sprintf('id%d',ii-1);
    
    % Structure begins
    s.glob.dst      = {'hdf', fullfile(p.prepare_data_path, p.prepare_data_filename)};
    s.glob.energy   = p.energy;
    s.glob.z        = p.z;
    s.glob.ds       = p.ds;
        
    s.detector.(cid)          = struct;   % id0 means "detector":{"0":
    s.detector.(cid).validMsk = p.detectors(ii).params.mask;  % Besides a matlab binary valid mask file
    % it accepts a definition based on
    % module and/or bad pixels
    s.detector.(cid).trfMsk = [0,0];     % Transformation for mask, [rot90,fliplr]
    
    % Mask for pilatus frames apparently not needed. Just datapath.
    % raw_data_filenamemask =fullfile(p.raw_data_path_full{1},sprintf('%s%05d_*.%s',p.detector.data_prefix,p.scan_number(ii),p.data_extension));

    s.measurement.(cid) = struct;
    s.measurement.(cid).detector  = ii-1; % NEEDS ADAPTATION FOR SHARING
    s.measurement.(cid).probe     = p.share_probe_ID(ii)-1; 
    s.measurement.(cid).object    = p.share_object_ID(ii)-1; 
    
    
    
    % extract h5location and reformat it 
    h5loc = [];
    for det_extraargs = 1:length(p.detectors(ii).params.image_read_extraargs)
        if strcmpi(p.detectors(ii).params.image_read_extraargs{det_extraargs}, 'H5Location')
            h5loc = p.detectors(ii).params.image_read_extraargs{det_extraargs+1}(2:end-1);
        end
    end
    
    % libDetXR uses 2 different readers: hdf for single h5 file and
    % oldEigerH5 for multiple files
    p = p.detectors(ii).params.get_filename(p);
    if strcmpi(p.detectors(ii).params.file_extension, 'h5')
        if numel(p.detectors(ii).detStorage.files)==1
            reader = 'hdf';
        else
            reader = 'oldEigerH5';
        end
    else
        reader = p.detectors(ii).params.file_extension;
    end
    
    % bug fix for reading cbfs
    if strcmpi(p.detectors(ii).params.file_extension, 'cbf')
        h5loc = 'eh5/images';
    end
    
    
    if numel(p.detectors(ii).detStorage.files)==1 && strcmpi(p.detectors(ii).params.file_extension, 'h5')
        fname = dir([p.raw_data_path_full{ii} '*.h5']);
        file_list = fullfile(fname.folder, fname.name);
    else
        file_list = p.raw_data_path_full{ii};
    end
    
    
    s.measurement.(cid).src = {reader, file_list, h5loc};

    % calculate transformation
    trfData = [0 0];
    if p.detectors(ii).params.orientation(1) == 1
        trfData = trfData+1;
    end
    if p.detectors(ii).params.orientation(2) == 1
        trfData(2) = trfData(2)+1;
    end
    if p.detectors(ii).params.orientation(3) == 1
        trfData = trfData+1;
        trfData(1) = trfData(1)+1;
    end
    
    % additional transpose for python data prep
    trfData = trfData - 1;
    
    % reduce number of fliplr
    trfData(2) = mod(trfData(2),2);
    
    s.measurement.(cid).trfData = trfData;
    assert(size(p.ctr,1) >= p.numscans, 'one line of p.ctrl is needed for each scan')
    s.measurement.(cid).roi       = [p.ctr(ii,2)-1, p.ctr(ii,1)-1, p.asize(1), p.asize(1)];
    s.measurement.(cid).pos       = round(p.positions(p.scanindexrange(ii,1):p.scanindexrange(ii,2),:));
    % Besides a list of points the positions
    % can be given as omni_pos dat file and in
    % that case it allows a transformation
    % matrix.
        
    s.initialCondition.dst      = {'hdf', fullfile(p.initial_conditions_path,p.initial_conditions_file)};
    %
    

    % Also supports hdf "0": ["hdf", "$INDIR/cSAXS_sxdm_2014_07_omny_commissioning/analysis/S00033/S00033_initial_conditions_400x400.h5","probes/probe_0"],
    % focusMdl and pinholeMdl are also
    % supported but not tested.

    
end

for ii=unique(p.share_probe_ID)
    cid = sprintf('id%d',ii-1);
    % Probe may be propagated or interpolated from initial guess, here we
    % have to save it to disk in order to be used by the Json ptycho
    % preparator
    
    % The probe is adapted in this code, interpolated, propagated or
    % created from a model. In order for the Python data preparer to access
    % these changes we save a temporary probe file to disk.

    verbose(3, 'Temporary probe file for python data prep = %s', p.initial_probe_temp_file{ii});
    temp_probe_full_file = fullfile(p.initial_conditions_path, p.initial_probe_temp_file{ii});
    %probe = squeeze(p.probes(:,:,1,:));
   
    probe = p.probes(:,:,ii,:);
    probe = permute(probe,[4 1 2 3]); % C-code expects probe mode index first.
    save(temp_probe_full_file,'probe');
    s.initialCondition.probe.(cid)    = {'mat', temp_probe_full_file, 'probe'};
    
end

for ii=unique(p.share_object_ID)
    cid = sprintf('id%d',ii-1);
    
    verbose(3, 'Temporary object file for python data prep = %s', p.initial_object_temp_file{ii});
    temp_object_full_file = fullfile(p.initial_conditions_path,p.initial_object_temp_file{ii});
    object = p.object{ii};
    object = squeeze(permute(object,[3 1 2])); % C-code expects probe mode index first.
    save(temp_object_full_file,'object', '-v6');
    s.initialCondition.object.(cid)   = {'mat', temp_object_full_file, 'object'};
end

% set parameters for c_solver engine; if only data preparation is needed,
% use dummy values
if p.external_engine0
    s.initialCondition.param.pfft_relaxation        = p.engines{1}.pfft_relaxation;
    s.initialCondition.param.probe_regularization   = p.engines{1}.probe_regularization;
    s.initialCondition.param.probe_radius           = p.engines{1}.probe_support_radius;
    s.initialCondition.param.diffmap_iterations     = p.engines{1}.number_iterations;
    s.initialCondition.param.max_mlh_iterations     = p.engines{1}.opt_iter;
else
    s.initialCondition.param.pfft_relaxation        = 0.05;
    s.initialCondition.param.probe_regularization   = 0.1;
    s.initialCondition.param.probe_radius           = 0.8;
    s.initialCondition.param.diffmap_iterations     = 200;
    s.initialCondition.param.max_mlh_iterations     = 100;
end

% Get a JSON string
t=mat2json(s);
% Write the file
if ~isfield(p, 'json_filename') || isempty(p.json_filename{1})
    p.json_filename{1} = [core.generate_scan_name(p) sprintf('_json_template_%03dx%03d.json', p.asize(1), p.asize(2))];
end
verbose(3, 'JSON template filename = %s', p.json_filename{1});
json_fullpath_filename = fullfile(p.initial_conditions_path,p.json_filename{1});

verbose(3,sprintf('Writting to file %s',json_fullpath_filename));
h_json = fopen(json_fullpath_filename,'w');
fprintf(h_json,t);
fclose(h_json);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Python prepare data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(p, 'ptychoPrep') || isempty(p.ptychoPrep)
    [status hostname] = system('hostname');
    fdb.status = core.engine_status(status);
    
    hostname = strsplit(hostname, '-');
    switch hostname{1}
        case 'ra'
            verbose(3, 'Loading python environment for the DaaS cluster.');
            p.ptychoPrep = fullfile(p.cSAXS_matlab_path, 'cSAXS_python_env/daas/lib/python2.7/site-packages/libDetXR/procPtycho.py');
        case 'x12sa'
            verbose(3, 'Loading python environment for the cSAXS beamline.');
            p.ptychoPrep = fullfile(p.cSAXS_matlab_path, 'cSAXS_python_env/x12sa/lib/python2.7/site-packages/libDetXR/procPtycho.py');
        otherwise
            error('Unknown host. Please specify your ptychoPrep or run it on the DaaS / beamline nodes.');
    end
end

if verbose > 4
    verb_ptychoPrep = ' -v1';
else
    verb_ptychoPrep = ' -v0';
end

python_call = ['python ' p.ptychoPrep ' --meta ' json_fullpath_filename verb_ptychoPrep ' --multiproc --cmpr s-zlib4'];

verbose(3,'Calling python prepare data:\n%s', python_call);
for ii = 1:3
    try
        [status, result] = system(python_call, '-echo');
        break
    catch ME
        warning('Data loading failed with error: %s\n Trying again ', ME.message)
    end
end
if status 
   error('Data loading failed') 
end

fdb.status = core.engine_status(status);

if p.delete_temp_json_data
    verbose(3,'Removing temporary probe and object files')
    for jj=1:length(p.initial_object_temp_file)
        if exist([p.initial_conditions_path p.initial_object_temp_file{jj}], 'file')
            [status, result] = system(['rm ' p.initial_conditions_path p.initial_object_temp_file{jj}], '-echo');
        end
    end
    for jj=1:length(p.initial_probe_temp_file)
        if exist([p.initial_conditions_path p.initial_probe_temp_file{jj}], 'file')
            [status, result] = system(['rm ' p.initial_conditions_path p.initial_probe_temp_file{jj}], '-echo');
        end
    end
else
    verbose(2,'Keeping H5 prepared data')
end



verbose(2,'Prepared data: %s', s.glob.dst{2});
verbose(2,'Initial conditions: %s', s.initialCondition.dst{2})


% bug fix for object sharing: procPtycho does not support object sharing,
% thus we have to overwrite the object size in the h5 data file. 
if any(p.share_object)
    h5_struc = [];
    for ii=1:size(p.object_size,1)
        h5_struc.objects(:,ii) = uint64(p.object_size(ii,:));
    end
    save2hdf5(fullfile(p.prepare_data_path, p.prepare_data_filename), h5_struc)
end


end

