%PTYCHO_SOLVER_DISTRIBUTED solve asynchronously ptychography, from
%"prepare_data_path" load parameters of the current reconstruction, data
%and from shared memory load the corresponding initial guess of the object 
% Finally, run ptychography and return results 
%
% [pout, ferr] = ptycho_solver_distributed(GPU_id, prepare_data_path)
% 
%  Inputs: 
%   **GPU_id                 - id of the used GPU 
%   **prepare_data_path      - path to the stored initialization files  
%
%  *returns*
%   ++pout                   - output from ptychography 
%   ++ferr                   - fourier error reported by ptychography

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


function [pout, ferr] = ptycho_solver_distributed(GPU_id, prepare_data_path)
  
    % fast loading to save time 
    %p = core.ptycho_recons(p, true);
    % IT IS 2X FASTER TO LOAD IT FROM CACHE
    load([prepare_data_path, '/prepared_data.mat'], 'p')  
    p.use_display = 0; 
    p.plot.interval = inf; 
    
    [p.fmag, p.fmask, p.positions] = io.load_prepared_data(fullfile( prepare_data_path, sprintf('S%05i_data_%ix%i.h5', p.scan_number, p.asize(1), p.asize(2)) )); 

          
    %% load initial guess + specific settings for this angles 
    init = load([prepare_data_path, '/prepared_initial_guess.mat']); 
    % get from the shared memory 
    shm_mem = shm(true, init.p.proj_id);
    [shm_mem, shm_data] = shm_mem.attach();
    init.p.object{1} = shm_data+eps;% force matlab to allocate new memory 

    
    % apply setttings from initial guess 
    for item = fieldnames(init.p)'
        if ~isempty(init.p.(item{1}))
            p.(item{1}) = init.p.(item{1}); 
        end
    end
    for item = fieldnames(init.engine0)'
        p.engines{1}.(item{1}) = init.engine0.(item{1}); 
    end
    
    
    p.gpu_id = GPU_id; 
    
    

    if isfield(p, 'gpu_id')
       gpu = gpuDevice();
       if p.gpu_id ~= gpu.Index
           gpuDevice(double(p.gpu_id)); 
       end
    end
    
    if strcmpi(p.engines{1}.name , 'GPU')
        p.object{1} = gpuArray(p.object{1}); 
        p.probes = gpuArray(p.probes);
    end
    
    
    %% run all engines
    utils.verbose(struct('prefix', {'ptycho'}))
    for ieng=1:length(p.engines)
        p.ieng = ieng;

        % engine call
        utils.verbose(1, 'Calling engine %s', p.engines{ieng}.name)
        utils.verbose(struct('prefix', {p.engines{ieng}.name}))
        
        %% run engine 
        
        try
            [pout, fdb] = core.run_engine(p,ieng);
        catch err 
            disp(err)
            keyboard
        end
        
        utils.verbose(struct('prefix', {'ptycho'}))
        
    end
    
    
        
    try
        ferr = pout.error_metric.value([1,end]); 
    catch
        ferr = nan(2,1); 
    end
    
    % delete input data, save time in storing  / loading 
    pout.fmag = []; 
    pout.fmask = []; 
    pout.probe = [];
    pout.simulation = []; 
   
    status = true; 
    for ii = 1:5
        try
            shm_mem.free;
            shm_mem.upload(pout.object{1})
            status = false;
            break
        catch err
            warning(err.message)
        end
        pause(0.1)
    end
    
    if status
        warning('Uploading results to shared mem failed')
        keyboard
    end
    
    % upload to shared memory
%      utils.add_to_3D_projection(shm_mem.data_handle{1}, real(pout.object{1}), [0,0], 1, false); 
%      utils.add_to_3D_projection(shm_mem.data_handle{2}, imag(pout.object{1}), [0,0], 1, false); 
     
    
    
    % delete the larger arrays from the p-structure
%     pout.probes = []; 
    pout.object = []; 
    
    % rewrite the structure with the initial guess and store there results 
    save([pout.  prepare_data_path, '/output_reconstruction.mat'] ,'ferr', 'pout', '-v6')
    utils.verbose(-1,'Save data %s', [pout.  prepare_data_path, '/output_reconstruction.mat'])
        
end
