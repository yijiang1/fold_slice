%%  TEST TEMPLATE FOR FUNTIONALITY CHECK IN GPU ENGINES 
% 1) call standard template to get fresh settings defaults 
% 2) generate artificial data that should serve as a standart test "sample"
% 3) call GPU engine with different basic functionalities and test if all still works
%   !! THESE TEST ARE ONLY USEFUL TO FIND CRASHES IN THE CODE, QUALITY OF THE RECONSTRUCTION IS NOT EVALUATED !!

%% set shared parameters for all test scripts 
run(fullfile( fileparts(mfilename('fullpath')), 'init_test.m'))

plot_results = true; 

%% general settings
p.   artificial_data_file = 'tests/test_data_shifted_objects.m';     % artificial data parameters 
p.   asize = [128 128];                              % size of the reconstruction probe 

p.   scan.lx = 10e-6;                                                 % round_roi scan: width of the roi
p.   scan.ly = 10e-6; 
p.   preserve_previous_results = true;         % do not reset p-struct after the first reconstruction 


%% ENGINES
% External C++ code

% --------- GPU engines  ------------- 
if 1
    eng = struct();
    eng. name = 'GPU';    
    eng. use_gpu = true;            % if false, run CPU code, but it will get very slow 
    eng. keep_on_gpu = true;        % keep data + projections on GPU, false is useful for large data if DM is used
    eng. compress_data = true;      % use automatic online memory compression to limit meed of GPU memory
    eng. gpu_id = [];               % default GPU id, [] means choosen by matlab
    eng. check_gpu_load = true;     % check available GPU memory before starting GPU engines 
    
    %% general 
    eng. number_iterations = 20;   % number of iterations for selected method 
    %eng. share_probe = 1;           % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
    %eng. share_object = 0;          % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
    eng. method = 'MLs';            % choose GPU solver: DM, ePIE, pPIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. opt_errmetric = 'L1' ;     % optimization likelihood - poisson, L1
    eng. grouping = 50;            % size of processed blocks, larger blocks need more memory but they use GPU more effeciently
                                    % for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                    % for pPIE, MLc the convergence is similar 
                                    % for DM, RAAR is has no effect on convergence
    %eng. probe_modes  = 1;          % Number of coherent modes for probe
    eng. object_change_start = 1;    % Start updating object at this iteration number
    eng. probe_change_start = 1;     % Start updating probe at this iteration number
    eng.probe_support_fft = false;

    % regularizations
    eng. reg_mu = 0;                % Regularization constant ( = 0 for no regularization)
    eng. delta = 0;                 % press values to zero out of the illumination area, usually 1e-2 is enough 
    eng. positivity_constraint_object = 0; % enforce weak positivity in object, usually 1e-2 is already enough 

    eng. apply_multimodal_update = false; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
    eng. probe_backpropagate = 0;         % backpropagate the probe mask, inf == farfield 
    
    % basic recontruction parameters 
    % PIE / ML methods
    eng. beta_object = 1;           % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. beta_probe = 1;            % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. delta_p = 0.1;             % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, 
    % DM
    eng. pfft_relaxation = 0.1;     % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = 0.1;% Weigth factor for the probe update (inertia)

    
    % ADVANCED OPTIONS   
    % position refinement 
    eng. apply_subpix_shift = false;       % apply FFT-based subpixel shift, important for good position refinement but it is slow
    eng. probe_position_search = inf;           % reconstruct probe positions, from iteration == probe_pos_search, assume they are independed
    eng. probe_position_error_max = 10e-9; % max expected random position error of the stages 

    % wavefront refinement                  
    eng. probe_fourier_shift_search = inf; % refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
    eng. estimate_NF_distance = inf;       % try to estimate the nearfield propagation distance  
    eng. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan
    eng. variable_SVD_modes = 3;           % OPRP settings , number of SVD modes, apply only for PIE methods 
    eng. variable_probe_smooth = 0;        % OPRP settings , apply assumption of smooth evolution of the OPRP modes (ie slow drifts)
    eng. variable_intensity = false;       % account to changes in probe intensity

    

    %% test MLc/s codes 
    eng. method = 'MLs';            % choose GPU solver: DM, ePIE, pPIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng.grouping = 50; 
    eng.probe_support_radius  = []; 
    eng.number_iterations = 50;
    eng.accelerated_gradients_start = inf; 
    
    eng.mirror_objects = false;
    eng.share_object = false; 

    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    
    eng.number_iterations = 50;

    eng.share_object = true; 
    eng.share_probe = false; 
    eng. align_shared_objects = true;     % before merging multiple unshared objects into one shared, the object will be aligned and the probes shifted by the same distance -> use for alignement and shared reconstruction of drifting scans  

    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    
    
    
end

run(fullfile( fileparts(mfilename('fullpath')), 'run_test.m'))


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

