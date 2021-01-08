%PREPARE_DISTRIBUTED_DATA use current tomographic volume reconstruction  volData and
% previous projection reconstructions projData_rec to create a new initial
% guess for the ML ptychography
%
% [projData_model, projData_rec, shm_mem] = prepare_distributed_data(p, volData, projData_rec, layer_distance, par, init_level, save_init_guess)
% 
%  Inputs: 
%   **p               -  ptycho p structure 
%   **volData         - 3D array, linearized tomographic volume (ie tranmission == exp(sum(volData,1)) == prod(exp(volData))  )
%   **projData_rec    - structure that contains complex projection reconstructed in the previous iteration 
%   **layer_distance  - vector, distance between layers in meters 
%   **par             - parameter structure for ptychotomo
%   **init_level      - flag: -1 == very first initialization, 0 = no initialization, 1 == force ptychography cache to update  
%
%  *returns*
%   ++projData_model  - model projection structure
%   ++projData_rec    - reconstructed projection structure  
%   ++shm_mem         - @shm class object that contain reference to the shared memory 


% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
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
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
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



function  [projData_model, projData_rec, shm_mem] = prepare_distributed_data(p, volData, projData_rec, layer_distance, par, init_level, save_init_guess)

    try
        %% create initial guess for ptychography, clear object_c if number of layers in increased 
        if isempty(projData_rec.object)
            projData_rec.object = exp(projData_rec.object_c); 
            warning('Object was missing in projData_rec')
        end
        Npx_proj = size(projData_rec.object); 
        

        Nlayers = length(layer_distance)+1; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%% GET PROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        projData_model = get_forward_model(projData_rec, volData, Nlayers); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % move on GPU 
        projData_rec.object_c = utils.Garray(projData_rec.object_c); 

        expand_layers = size(projData_rec,3) < Nlayers; 
        if expand_layers
            % cheaper version but less exact 
            projData_rec.object_c = repmat(sum(projData_rec.object_c,3),1,1,Nlayers)/Nlayers; 
            projData_rec.object = [];
        end

        
        % find reliability region of the current model 
        Npx_proj_small = size(projData_model.object_c); 
        win = utils.Garray(single(tukeywin(Npx_proj_small(2), 0.2))'.*single(tukeywin(Npx_proj_small(1), 0.2))); 
        win = utils.crop_pad(win, Npx_proj(1:2));
        win_shifted = utils.imshift_fast(win, projData_rec.position_offset(1), projData_rec.position_offset(2)); 

      
        weight = utils.Garray(projData_rec.weight) ; 
        if init_level == -1 % first initialization 
            % estimate air region 
             W_tmp = weight .* exp(-abs(sum(projData_model.object_c,3))) .* win_shifted; 
             % remove offset between model and the reconstruction. 
             projData_rec = remove_offset_initial(projData_model,projData_rec, W_tmp, Npx_proj_small); 
            % run only the initialization  
            return 
        end

         
        %% merge phase and amplitude from tomo and measurements 
        projData_model.object_c =   win_shifted .*projData_model.object_c  + (1-win_shifted).* projData_rec.object_c;

        if ~isempty(weight)
            projData_rec.object_c   = projData_rec.object_c .*  (weight~=0); 
            projData_model.object_c = projData_model.object_c .*  (weight~=0); 
        end


    catch err
        err
        keyboard
    end
    
   
        

    
    %% get the complex transmission 
    object = exp(projData_model.object_c);
    projData_model.object = object; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STORE PREPARED PROJECTION FOR ASYNCHRONOUS SOLVERS %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    object = gather(object);
    
    if (init_level && par.force_initialization) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % precache the already initalized p structures to make the solver faster %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        p.ds = []; 
        p.use_gpu = true; 
        p.scan = []; 
        p.  scan_number = projData_rec.scan_id;
        p.  rotation_angle = projData_rec.angle; 
        p. proj_id = projData_rec.proj_id; 
        
        %% prepare data 
        p.  prepare_data_path = sprintf(par.prepare_data_path, p.scan_number);      
        
        p.object = []; 
        p.probe = []; 
        p.probes = []; 
        p.positions_real = []; 
        p.positions_orig = []; 
        %% store data into some fast cache 
        save('-v6', [p.prepare_data_path, '/prepared_data.mat'], 'p')
    end
    

    %% prepare the initial guess and store it 
    
    p = struct(); 
    p. proj_id = projData_rec.proj_id; 
    p.  scan_number = projData_rec.scan_id;
    p.  rotation_angle = projData_rec.angle; 
    p.  prepare_data_path = sprintf(par.prepare_data_path, p.scan_number); 
    p.  force_preparation_data = false; 
    p.  remove_object_ambiguity= false; 
    p.  initial_probe_file = ''; 
    p.  prefix = ''; 
    p.  run_name = '' ; 
    p.  save_path = ''; 
    
    % the input guess should be already optimally shifted 
    p.  position_offset = [0,0]; 
    
    
    if utils.verbose() 
        engine0.use_display = 1; 
        p.  use_display = 0; 
        p.  verbose_level = 3; 
        engine0.verbose_level = 3; 
    else
        p.  use_display = 0; 
        p.  verbose_level = -1; 
        engine0.verbose_level = -1; 
    end
    
    % initialize the reconstruction, load prepared data
    engine0.delta_z = layer_distance;
    engine0.N_layer = size(object,3);
    

    if init_level == 1
        engine0.initial_probe_rescaling = true;
    else
        engine0.initial_probe_rescaling = false;
    end
    
    
    
    %% load initial guesses 
    p.probes = single(projData_rec.probe); 
    p.positions = projData_rec.positions; 

    
        
    % ptychography expects the layers in opposite order
    object = object(:,:,end:-1:1); 
      
    p.object{1} = gather(reshape(object, [size(object,1), size(object,2),1,size(object,3)])); 
    
    p.object_size = [size(object,1), size(object,2)]; 
    
    ferr = nan;  % fourier error is initialized as nan 
    
    % prevent saving useless variables 
    try; p = rmfield(p, 'simulation'); end
    try; p = rmfield(p, 'probe'); end
    
    
    if save_init_guess

        if ~exist(p.prepare_data_path, 'dir')
            mkdir(p.  prepare_data_path); 
        end
             
        
        % upload to the shared memory 
        shm_mem = shm(false, p.proj_id);
        shm_mem.upload(p.object{1});
        shm_mem.detach
                
        % delete from p-struct , it is already in the shared memory 
        p.object = []; 
        % save only the empty p-structure 
        save([p.  prepare_data_path, '/prepared_initial_guess.mat'] , 'p', 'ferr', 'engine0', '-v6')
    else
        shm_mem = []; 
    end


    
end

function projData_model = get_forward_model(projData_rec, volData, Nlayers)


    Npx_vol = size(volData); 
    Nblock = ceil(prod(Npx_vol)*8 / 1e9);  % split into 2GB blocks 
    Npx_proj = size(projData_rec.object);
    
    % block by block get the projection model 
    projData_model = projData_rec; % create a copy of the input structure 
    projData_model.object_c = gpuArray.zeros([Npx_vol([3,1]),Nlayers], 'single'); 
    for ii = 1:Nblock
        ind = 1+(ii-1)*ceil(Npx_vol(3)/Nblock):min((ii)*ceil(Npx_vol(3)/Nblock), Npx_vol(3)); 

        volData_block = volData(:,:,ind); 
        volData_block = utils.Garray(volData_block); 

        [object_r] = ptychotomo.fwd_proj(volData_block, projData_rec.angle, Nlayers, 'real');
        [object_i] = ptychotomo.fwd_proj(volData_block, projData_rec.angle, Nlayers, 'imag');

        projData_model.object_c(ind,:,:) = complex(object_r, object_i); 

    end



    % pad the geprojData_modelnerated projection to the size of the measured
    % projections 
    projData_model.object_c = utils.crop_pad(projData_model.object_c,Npx_proj(1:2),0); 

     if any(any(round(projData_rec.position_offset) ~= projData_rec.position_offset))
        error('Noninterger shift may result in mixing real and imaginary data')
     end


    % apply shifts found by prealignement 
    projData_model.object_c  = utils.imshift_fast(projData_model.object_c, ...
                                                    projData_rec.position_offset(1), ...
                                                    projData_rec.position_offset(2)); 

        
end


function [projData_rec, offset] = remove_offset_initial(projData_model,projData_rec, weight, Npx_proj_small)
    import utils.*
    
    W = crop_pad(weight, Npx_proj_small) ; 
    diff = crop_pad(sum(projData_model.object_c,3), Npx_proj_small) - crop_pad(sum(projData_rec.object_c,3), Npx_proj_small); 
    % empirical estimation of weights to subtract the phase / amplitude offset

    offset   = mean(mean(mean(W.*diff,2))) ./ mean(mean(W,2));
    
    % subtract offset plane 
    projData_rec.object_c = projData_rec.object_c + offset  ; 
    
    verbose(0,'Offset removal  %3.2e+%3.2ei ', real(offset), imag(offset)); 
    
end
