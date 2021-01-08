%% TEST TEMPLATE FOR FUNTIONALITY CHECK OF CPU ENGINES 
% 1) call standard template to get fresh settings defaults 
% 2) generate artificial data that should serve as a standart test "sample"
% 3) call GPU engine with different basic functionalities and test if all still works
%   !! THESE TEST ARE ONLY USEFUL TO FIND CRASHES IN THE CODE, QUALITY OF THE RECONSTRUCTION IS NOT EVALUATED !!


%% set shared parameters for all test scripts 
run(fullfile( fileparts(mfilename('fullpath')), 'init_test.m'))
close all 

plot_results = true;
plot_save_results = true; 


%% general settings
p.   artificial_data_file = 'tests/test_data_double_scan.m';     % artificial data parameters 
p.   asize = [128 128];                              % size of the reconstruction probe 
Nlayers = 2; 


%% Plot, save and analyze

p.   plot.prepared_data = false;                         % plot prepared data
p.   plot.interval = [];                                    % plot each interval-th iteration, does not work for c_solver code
p.   plot.log_scale = [1 1];                                % Plot on log scale for x and y
p.   plot.realaxes = true;                                  % Plots show scale in microns
p.   plot.remove_phase_ramp = true;                        % Remove phase ramp from the plotted / saved phase figures 
p.   plot.fov_box = true;                                   % Plot the scanning FOV box on the object (both phase and amplitude)
p.   plot.fov_box_color = 'r';                              % Color of the scanning FOV box
p.   plot.positions = true;                                 % Plot the scanning positions
p.   plot.mask_bool = true;                                 % Mask the noisy contour of the reconstructed object in plots
p.   plot.windowautopos = true;                             % First plotting will auto position windows
p.   plot.obj_apod = false;                                 % Apply apodization to the reconstructed object;
p.   plot.prop_obj = 1e-12;                                     % Distance to propagate reconstructed object before plotting [m]
p.   plot.show_layers = true;                              % show each layer in multilayer reconstruction 
p.   plot.show_layers_stack = true;                        % show each layer in multilayer reconstruction by imagesc3D
p.   plot.object_spectrum = true;                          % Plot propagated object (FFT for conventional ptycho); if empty then default is false if verbose_level < 3 and true otherwise
p.   plot.conjugate = true;                                % plot complex conjugate of the reconstruction 
p.   plot.horz_fact = 2.5;                                  % Scales the space that the ptycho figures take horizontally
p.   plot.FP_maskdim = 180e-6;                              % Filter the backpropagation (Fourier Ptychography)
p.   plot.calc_FSC = true;                                 % Calculate the Fourier Shell correlation for 2 scans or compare with model in case of artificial data tests 
p.   plot.show_FSC = true;                                 % Show the FSC plots, including the cropped FOV
p.   plot.residua = true;                                  % highlight phase-residua in the image of the reconstructed phase



%% ENGINES

if isunix
    % Please notice that you have to force data preparation (force_prepare_h5_files=true) if you have made any changes to 
    % the already prepared data (fmag, fmask, positions, sharing ...). 
    eng.  name = 'c_solver';
    eng.  method = 'DM+ML'; 
    eng.  number_iterations = 100;              % Total number of iterations
    eng.  opt_iter = 100;                       % Iterations for optimization     
    eng.  probe_regularization = .1;            % Weigth factor for the probe update; 
    eng.  probe_change_start = 1;               % Start updating probe at this iteration number
    eng.  probe_support_radius = 0.8;           % Normalized radius of circular support, = 1 for radius touching the window    
    eng.  pfft_relaxation = .05;                % Relaxation in the Fourier domain projection, = 0  for full projection    
    eng.  background = 0;                       % [PARTIALLY IMPLEMENTED (not fully optimized)] Add background to the ML model in form:  |Psi|^2+B, B is in average counts per frame and pixel
    eng.  probe_support_fft = false;            % [PARTIALLY IMPLEMENTED (not fully optimized)] Apply probe support in Fourier space, ! uses model zoneplate settings to estimate support size 

    eng.  N_layer =  Nlayers;                          % Number of virtual object layers (slices)
    eng.  delta_z = 10e-6 * ones(1, eng.N_layer-1); % Separation between object slices 
    %eng.  ms_init_ob_fraction = [1 0];
    if eng.  N_layer>1
        p.suffix = [p.suffix '_N' num2str(eng. N_layer)];
        eng.  number_iterations = 0; % highly recommended
    end
    
    eng.  single_prec = true;                   % single or double precision
    eng.  threads = 20;                         % number of threads for OMP
    eng.  beamline_nodes = [];                  % beamline nodes for the MPI/OMP hybrid, e.g. ['x12sa-cn-2'; 'x12sa-cn-3'];
    eng.  ra_nodes = 0;                         % number of nodes on ra cluster for the MPI/OMP hybrid; set to 0 for current node
    eng.  caller_suffix = '';                   % suffix for the external reconstruction program
    eng.  reconstruction_program = '';          % specify external reconstruction program that overwrites previous settings, e.g. 'OMP_NUM_THREADS=20 ./ptycho_single_OMP';
    eng.  check_cpu_load = true;                % check if specified nodes are already in use (only x12sa). Disable check if you are sure that the nodes are free.
    eng.  initial_conditions_path = '';         % path of the initial conditions file; default if empty (== prepare_data_path)
    eng.  initial_conditions_file = '';    		% Name of the initial conditions file, default if empty. Do not use ~ in the path
    eng.  measurements_file = '';				% Name of the measurements file, default if empty. Do not use ~ in the path
    eng.  solution_file = '';                   % Name of the solution file, default if empty. Do not use ~ in the path
    eng.  force_prepare_h5_files = 1;           % If true before running the C-code the data h5 file is created and the h5 file with initial object and probe too, regardless of whether it exists. It will use the matlab data preparator. 
    [p, ~] = core.append_engine(p, eng);        % Adds this engine to the reconstruction process
else
    eng = struct();
    eng. name = 'DM';
    eng. method = 'matlab'; 
    eng. number_iterations = 10;               % Total number of iterations
    eng. probe_change_start = 1;              % Start updating probe at this iteration number
    eng. average_start = 300;                 % Start averaging at this iteration number
    eng. average_interval = 5;                % Number of iterations between reconstruction estimates for average 
    eng. count_bound = 4e-2;                  % Relaxed Fourier projection parameter - average photons of change per pixel (= 0 no relaxation) 
    eng. pfft_relaxation = 0.05 ;               % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = .1;           % Weigth factor for the probe update
    eng. probe_mask_bool = true;              % If true, impose a support constraint to the probe
    eng. probe_mask_area = .9;                % Area ratio of the mask
    eng. probe_mask_use_auto = false;         % Use autocorrelation for probe_mask (if false: circular circle)
    eng. object_flat_region = [];             % Mask for enforcing a flat region in the object (to reduce artifacts)
    eng. remove_scaling_ambiguity = true;     % Remove ambiguity of the probe times object scalling by probe normalization
    eng. clip_object = true;                  % Clip the object transmission function
    eng. clip_max = 1.0;                      % Upper bound
    eng. clip_min = 0.0;                      % Lower bound
    eng. compute_rfact = false;               % If set to true, R-factor is computed at every iteration (large overhead!!!)
    eng. use_mex = [0,0,0]; 
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end


    
run(fullfile( ptycho_path, 'tests/run_test.m'))



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
