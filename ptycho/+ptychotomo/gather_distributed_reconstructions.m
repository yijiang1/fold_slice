%GATHER_DISTRIBUTED_RECONSTRUCTIONS gather object reconstructions from the shared
% memory, calculate difference from the initial reconstruction, and apply
% update into the full tomographic volume. 
%
% [volData,projData_new, fourier_error, update_norm] = ...
%                     gather_distributed_reconstructions(volData, projData_model, ptycho_results, par)
% 
%  Inputs: 
%   **volData         - 3D array, linearized tomographic volume (ie tranmission == exp(sum(volData,1)) == prod(exp(volData))  )
%   **projData_model  - structure that contain complex projection, initial guess and other values related to the currently processed angle 
%   **ptycho_results  - output from ptychography, if empty, load it from the share memory array 
%   **par             - parameter structure for ptychotomo
%
%  *returns*
%   ++volData         - updated tomographic volume 
%   ++projData_new    - updated projection structure  
%   ++fourier_error   - vector, fourier error reported by ptychography 
%   ++update_norm     - scalar, difference between original and new tomo volume


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


function [volData,projData_new, fourier_error, update_norm] = ...
                    gather_distributed_reconstructions(volData, projData_model, ptycho_results, par)
    
                
    % prepare path 
    p.prepare_data_path = sprintf(par.prepare_data_path, projData_model.scan_id);
    Npx_proj = size(projData_model.object); 
    

    %% GATHER DATA FROM SHARED MEMORY 
    if ~isempty(ptycho_results) && isstruct(ptycho_results)
        % avoid loading from disk if the cached data are available 
        data.pout = ptycho_results; 
        data.ferr = nan; 
    else
        % load reoonstruction from disk, calculate optimal update for given angle, use it to calculate
        % optimal volume update 

        for ii = 1:3
            % try twice before failing 
            try
                data = load([p.  prepare_data_path, '/output_reconstruction.mat'], 'pout', 'ferr');
                assert(isfield(data, 'ferr'), 'Results not loeaded correctly')
                break
            end
            pause(0.5)
        end
        if ~exist('data', 'var')
            error('Reconstruction %s not found', p.  prepare_data_path)
        end
        for ii = 1:3
            try
                % gather to the shared memory 
                [ptycho_results, data_shm] = ptycho_results.attach();
                break
            end
            pause(0.5)
        end
        if isempty(ptycho_results)
            keyboard
        end
        if ~exist('data_shm', 'var')
            error('Shared data not found')
        end
        
       data.pout.object{1} = data_shm;% force matlab to allocate new memory 
       if numel(data_shm) ~= prod(Npx_proj)
           warning('Wrong object size returned')
           keyboard
       end
       if isreal(data_shm) || any(all(all(imag(data_shm)==0)))
          warning('Something went wrong, skipping projection update')
          % something wrong happened with the data in the shared memory, it
          % seems that the complex part was lost at least partly -> skip
          % this update and try next projection 
          projData_new = projData_model; 
          fourier_error= nan; 
          update_norm= nan; 
          return
       end
       
       clear  ptycho_results

    end

    %% GET UPDATE FROM THE COMPLEX VALUED PROJECTION 
    
    projData_new = projData_model;  % make a new structure by a compy of the previous one 
    % return sorted for ptychography
    projData_new.object = utils.Garray(squeeze(data.pout.object{1}(:,:,:,end:-1:1))); 
    projData_new.probe = squeeze(data.pout.probes); 
    projData_new.positions = data.pout.positions; 
    
    % enforce full transmission in the nonmeasured regions 
    missing = (projData_new.weight)==0; 
    projData_new.object = projData_new.object .* ~missing + missing; 

    % compare the model object_c and the updated object, use the difference to find new object_c 
    projData_new.object_c  = ptychotomo.prepare_projections(projData_new.object, Npx_proj, data.pout.asize, false, projData_model.object_c,  projData_model.weight);
    projData_new.object    = gather(projData_new.object); 
    projData_new.weight    = gather(projData_new.weight); 
        
    
    projData_c_upd = gpuArray(projData_new.object_c - projData_model.object_c); 
    
    
    % subtract potential offset in the global phase or aplitude  
    offset_tot = 0; 
    for ii = 1:3
        offset = mean(mean(mean((projData_c_upd),3) .*  projData_new.weight)) ./ mean(mean( projData_new.weight)); 
        projData_c_upd = (projData_c_upd-offset) .* projData_new.weight; 
        offset_tot = offset_tot + offset; 
    end
    utils.verbose(1,'Offset removal gather  %3.2e+%3.2ei', real(offset_tot), imag(offset_tot));
  
    

        
    if 2*sum(math.norm2(projData_c_upd)) > 0.5
        warning('Probably some convergence issue')
        keyboard
    end
    
    %%%%%%%%%%%%%%%%%% BACKPROJECT THE UPDATE TO THE VOLUME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [volData, update_norm] = ptychotomo.update_volume(volData, projData_c_upd, par.update_step, projData_model, par);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if update_norm > 0.1
        warning('Probably some convergence issue')
        keyboard
    end
    
    % store the filtered / regularized version for the next iteration  guess 
    projData_new.object_c = gather((projData_model.object_c + projData_c_upd));
    % already a version with corrections 
    projData_new.object = exp(projData_new.object_c); 
    
   
    
  
try
    fourier_error = data.ferr ; 
catch
    fourier_error = nan(2,1); 
end

   
              
            
end


