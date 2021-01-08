%C_SOLVER external C++ code for DM and ML
% ** p      p structure
%
% returns:
% ++ p      p structure
% ++ fdb    feedback structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
%
%   Publications most relevant to the Difference-Map implementation
%       + P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, 
%       "High-Resolution Scanning X-ray Diffraction Microscopy," Science 321, 379-382 (2008)
%       + P. Thibault, M. Dierolf, O. Bunk, A. Menzel, F. Pfeiffer,
%       "Probe retrieval in ptychographic coherent diffractive imaging,"
%       Ultramicroscopy 109, 338–343 (2009)
%
%   Publications most relevant to the Maximum Likelihood refinement
%       + M. Guizar-Sicairos and J. R. Fienup, "Phase retrieval with transverse
%       translation diversity: a nonlinear optimization approach," Opt. Express 16, 7264-7278 (2008)
%       + P. Thibault and M. Guizar-Sicairos, "Maximum-likelihood refinement for
%       coherent diffractive imaging," New J. Phys. 14, 063004 (2012).

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

function [ p , fdb] = c_solver( p )
import utils.verbose
import utils.crop_pad
import io.HDF.*
import io.load_prepared_data
import beamline.identify_eaccount

fdb.status = [];
ob = p.object;

% make sure all the fields are defined

% make sure all the fields are defined
if ~isfield(p, 'probe_change_start')
    p.probe_change_start = 1;
end
if ~isfield(p, 'preshift_ML_probe')
    p.preshift_ML_probe = true;
end
if ~isfield(p, 'delta_z')
    p.delta_z = [];
end
N_layer = 1+length(p.delta_z);


% store the precise positions in the p-structure before padding
positions_float = p.positions;

if any(p.positions_pad~=0)
    verbose(2, 'Padding positions is not supported by the C++ code. I will go back to the original size...')
    
    for obnum=1:p.numobjs
        p.object_size(obnum,:) = p.object_size(obnum,:)-2*p.positions_pad;
        ob{obnum} = crop_pad(p.object{obnum},p.object_size(obnum,:));
    end
    p.positions = p.positions - p.positions_pad;
    p.positions = round(p.positions);
end


% Parameters for reconstruction with external C code
if isempty(p.initial_conditions_file)
    if ~isempty(p.suffix)
        suffix = ['_' p.suffix];
    else
        suffix = '';
    end
    p.initial_conditions_file = [core.generate_scan_name(p) sprintf('_initial_conditions_%03dx%03d%s.h5', p.asize(1), p.asize(2), suffix)];
    verbose(3, 'C-code initial_conditions_file = %s', p.initial_conditions_file);
end


if isempty(p.solution_file)
    for ii = 1:length(p.scan_number)
        p.solution_file{ii} = [p.run_name '_c.h5'];
        verbose(3, 'reconstruction filename = %s', p.solution_file{ii});
        if ~isempty(p.save_path{1})&&(p.save_path{ii}(1) == '~')
            p.save_path_c{ii} = ['/sls/X12SA/Data10/' identify_eaccount p.save_path{ii}(9:end)];
        else
            p.save_path_c{ii} = p.save_path{ii};
        end
    end
    
end

% path for saving temp data for C code
if isfield(p, 'initial_conditions_path') && ~isempty(p.initial_conditions_path)
    if p.prepare_data_path(1) == '~'
        p.initial_conditions_path = ['/sls/X12SA/Data10/' identify_eaccount p.prepare_data_path(9:end)];
    end
else
    p.initial_conditions_path = p.prepare_data_path;    
end

if p.current_engine_id > 1
    % recalculate object and object_size in case they have changed
    p.positions = round(p.positions);
    p.positions = p.positions - min(p.positions);
    p.object_size = p.asize + max(round(p.positions),[],1);
    for ii = unique(p.share_object_ID)
        p.object{ii} = crop_pad(double(p.object{ii}), p.object_size(ii,:));
    end
end

            
% Write hdf5 files if forced or if it is not the first engine 
if p.force_prepare_h5_files || p.current_engine_id > 1
    core.prep_h5data(p);
end


if ~ (  strcmpi(p.engines{1}.name, 'c_solver') && strcmpi(p.prepare.data_preparator, 'libDetXR'))
    % update initial conditions always except the case when c_solver is the
    % first engine and data preparator is libDetXR (python)
    engines.c_solver.prep_initial_conditions(p);
end
        
%%%%%
%%% Check if reconstruction name exist and append a number to avoid overwrite %%%
filename_with_path = fullfile(p.save_path_c{1}, p.solution_file{1});
if exist(filename_with_path, 'file')
    verbose(3,'File %s exists!', filename_with_path);
    alt_filename = filename_with_path;
    [~, fbase,f2] = fileparts(filename_with_path);
    append_number = 0;
    while exist(alt_filename, 'file')
        f1 = sprintf('%s_%02d', fbase, append_number);
        alt_filename = fullfile(p.save_path_c{1}, [f1 f2]);
        append_number = append_number + 1;
    end
    filename_with_path = alt_filename;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run external code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%


%%% prepare external call
[p, fdb] = engines.c_solver.prepare_external_call(p, fdb);

%%% call C++ code
ctic = tic;

if N_layer>1
    st_delta_z = num2str(p.delta_z(1)); % slice_distances
    for ii = 2:length(p.delta_z)
        st_delta_z = [st_delta_z ':' num2str(p.delta_z(ii))];
    end
else
    st_delta_z = '-';
end
        
if verbose < 3
    feedback_interval = 25;
elseif verbose == 3
    feedback_interval = 10;
else
    feedback_interval = 5;
end

if verbose > 1
    c_verbose = 1;
else
    c_verbose = 0;
end
if verbose > 4
    debug_flag = 8+64+2048;
elseif verbose >= 0
    debug_flag = 64;
else
    debug_flag = 0;
end



probe_support_string = ''; 
if check_option(p, 'probe_support_fft')
    % calculate optimal support size in fourier space if focused beam is used 
    if ~check_option(p.model, 'probe_focal_length') && ~check_option(p.model, 'probe_outer_zone_width')
        error('Missing  model.probe_focal_length and model.probe_outer_zone_width of Fresnel zone plate' )
    end        
    if ~check_option(p.model, 'probe_outer_zone_width')
        p.model.probe_outer_zone_width = p.lambda * p.model.probe_focal_length / p.model.probe_diameter; 
    end
    FZP_cone_diameter = p.lambda* p.z/(p.model.probe_outer_zone_width * p.ds); 
    % add some extra space 
    FZP_cone_diameter = FZP_cone_diameter * 1.2; 
    
    [cx, cy] = math.center(abs(fftshift(fft2(p.probes(:,:,1,1))))); 
 
    signal_radius = FZP_cone_diameter / p.asize(1); 
    signal_center_row = 0.5 + (cx / p.asize(1));
    signal_center_column = 0.5+ (cy /  p.asize(1)); 
    probe_support_string =  [' --signal_radius=',num2str(signal_radius), ' --signal_center_row=', num2str(signal_center_row), ' --signal_center_column=', num2str(signal_center_column)]; 
end

background_string = ''; 
if  check_option(p, 'background')
    if ~isfield(p, 'renorm')
       error('FIXME: background needs p.renorm which is know only when the data are prepared by matlab_ps or first engines is not external as c_solver')
    end
    background_string = [' --background_correction=' num2str(p.background * p.renorm^2)] ; 
end
c_propagator = '';
if check_option(p, 'propagator')
    c_propagator = [' --propagator=' p.propagator ' --detector_distance=' num2str(p.z)];
end

c_propagator = '';
if check_option(p, 'propagator')
    c_propagator = [' --propagator=' p.propagator ' --detector_distance=' num2str(p.z)];
end

external_call = [p.reconstruction_program ...
    ' --debug_flags=$((' num2str(debug_flag) '))' ... % 8+64+2048 (8: overall execution info, 64: timing info, 256: info on given arguments, 2048: max likelihood function values feedback).
    ' --feedback_interval=' num2str(feedback_interval) ...
    ' --verbose=' num2str(c_verbose) ...
    ' --diffmap_iterations=' num2str(p.number_iterations) ...
    ' --max_mlh_iterations=' num2str(p.opt_iter) ...
    ' --probe_modes=' num2str(p.probe_modes) ...
    ' --object_modes=' num2str(p.object_modes) ...
    ' --compress=' num2str(p.io.file_compression) ...
    ' --wavelength=' num2str(1.2398e-9/p.energy) ... % (required for multislice)
    ' --drow=' num2str(p.dx_spec(1)) ...    % (taken from initial_conditions file if present, command line takes precedence)   
    ' --dm_fixed_probe_iter=' num2str(p.probe_change_start) ...
    ' --slice_distances=' st_delta_z ... 
    ' --num_slices=' num2str(N_layer) ...
...  ' --mlh_distopt_interval=1 ' ' --mlh_distadj_interval=50 '... 
    background_string ... 
    probe_support_string ... 
    c_propagator ...
    ' ' p.prepare_data_path p.prepare_data_filename ' ' ...
    p.initial_conditions_path p.initial_conditions_file ' ' filename_with_path];

verbose(3, 'Calling external program\n%s', external_call);
[status, result] = system(external_call, '-echo');


fdb.status = core.engine_status(status);
if status ~= 0
    verbose(0,'External program reported an error!');
    return
end
p.recon_filename_c = filename_with_path;
ctoc = toc(ctic);
verbose(3, 'Elapsed time for external call: %f', ctoc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of external code %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the solution file
solution_data = hdf5_load(filename_with_path);

% Change to structure format used in the rest of the code
for obnum = 1:p.numobjs
    object_ptr = solution_data.objects.(['object_' num2str(obnum-1)]);
    p.object{obnum} = permute(complex(object_ptr.r,object_ptr.i), [2,1,3,4]);
end

p.probes = []; 
for prnum = 1:p.numprobs
    probe_ptr = solution_data.probes.(['probe_' num2str(prnum-1)]);
    p.probes(:,:,prnum,:) = permute(complex(probe_ptr.r,probe_ptr.i), [2,1,3,4]);
end


if N_layer > 1 && p.preshift_ML_probe
   % if multilayer extension is used, shift the probe to be
   % reconstructed at the center plane of the sample -> remove shift 
   probe_offset = +sum(p.delta_z)/2; 
   p.probes = utils.prop_free_nf(p.probes, p.lambda , probe_offset, p.dx_spec(1)) ;
end



if p.number_iterations == 0 && ~all(solution_data.feedback.max_likelihood.iteration==0)
    verbose(0,'ML solver ended after zero iteration')
end

if isfield(solution_data.feedback, 'difference_map') && isfield(solution_data.feedback, 'max_likelihood') && ~all(solution_data.feedback.max_likelihood.iteration==0)
    p.error_metric{1} = solution_data.feedback.difference_map;
    p.error_metric{2} = solution_data.feedback.max_likelihood;
    p.error_metric{1}.method = 'DM';
    p.error_metric{1}.err_metric = 'RMS';
    p.error_metric{2}.method = 'ML';
    p.error_metric{2}.err_metric = '-LogLik';
elseif isfield(solution_data.feedback, 'difference_map')
    p.error_metric = solution_data.feedback.difference_map;
    p.error_metric.method = 'DM';
    p.error_metric.err_metric = 'RMS';
elseif isfield(solution_data.feedback, 'max_likelihood') && ~all(solution_data.feedback.max_likelihood.iteration==0)
    p.error_metric = solution_data.feedback.max_likelihood;
    p.error_metric.method = 'ML';
    p.error_metric.err_metric = '-LogLik';
else
    verbose(2,'Missing feedback option, returning empty feedback')  % e.g. in case of too low verbosity or low number of iteration 
    p.error_metric.iteration = [];
    p.error_metric.value = [];
    p.error_metric.method = 'DM';
    p.error_metric.err_metric = 'RMS';
end


% load the data
if (p.external_engine0 && strcmpi(p.prepare.data_preparator, 'python')) || (p.external_engine0 && ~p.prepare.force_preparation_data)
    [p.fmag, p.fmask, p.positions, max_power] = io.load_prepared_data([p.prepare_data_path p.prepare_data_filename]);
    p.renorm = sqrt(1/max_power);
    p.Nphot = sum((p.fmag(:)/p.renorm).^2.*p.fmask(:));
    p.fmask_per_scan = (length(size(p.fmask)) == 3);
end

if any(p.positions_pad~=0)
    verbose(2, 'Reapplying padding...');
    for obnum=1:p.numobjs
        p.object_size(obnum,:) = p.object_size(obnum,:)+2*p.positions_pad;
        p.object{obnum} = crop_pad(p.object{obnum},p.object_size(obnum,:), 1e-5);
    end
end
p.positions = positions_float;  % keep the precise positions in the p-structure 


% delete h5 file if it is not needed anymore
if ~(p.current_engine_id == length(p.engines))  % check that it is the last engine 
    verbose(3, 'Removing h5 file...')
    delete(filename_with_path)
elseif strcmpi(p.save.output_file, 'h5') || strcmpi(p.save.output_file, 'cxs')
    % if the c_solver is the last engine, keep the h5 file but move the
    % data to group reconstruction and delete the attributes
    move_data_h5(filename_with_path);
    hdf5_rm_attr(filename_with_path, '/', {'max_mlh_iterations'; ...
        'probe_modes'; 'object_modes'; 'pfft_relaxation'; ...
        'probe_regularization'; 'diffmap_iterations'; 'probe_radius'});
    
    attr.MATLAB_class = 'cell';
    io.HDF.hdf5_append_attr(filename_with_path, attr, '/reconstruction/p/objects');
    attr.MATLAB_class = 'complex';
    for ii=0:p.numobjs-1
        io.HDF.hdf5_append_attr(filename_with_path, attr, ['/reconstruction/p/objects/object_' num2str(ii)]);
    end
    for ii=0:p.numprobs-1
        io.HDF.hdf5_append_attr(filename_with_path, attr, ['/reconstruction/p/probes/probe_' num2str(ii)]);
    end
    
end

end

function move_data_h5(filename)
    % move data in h5 file to reconstruction group
    io.HDF.hdf5_mv_data(filename, 'feedback', 'reconstruction/feedback');
    io.HDF.hdf5_mv_data(filename, 'objects', 'reconstruction/p/objects');
    io.HDF.hdf5_mv_data(filename, 'probes', 'reconstruction/p/probes');
end

