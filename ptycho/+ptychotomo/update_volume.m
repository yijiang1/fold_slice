%UPDATE_VOLUME apply the updated projection information to the shared
%tomographic volume using a gradient descent method 
%
% [volData, update_norm] = update_volume(volData, projData_c_upd, update_step, projData_model, par)
% 
% 
%  Inputs: 
%   **volData      - 3D array, linearized tomographic volume (ie tranmission == exp(sum(volData,1)) == prod(exp(volData))  )
%   **projData_upd - stacked 2D array, linearized difference of the complex-valued projections, ie log( P_new .* conj(P_orig) ./ |P_orig|^2)
%   **updated_step - complex scalar, update step length for the gradient descent method, should be < 1 
%   **projData_model - structure that contain information about the model projection at the current angle 
%   **par          - parameter structure for the ptychotomo method 
% *returns*
%   **volData      - 3D array, updated linearized tomographic volume (ie tranmission == exp(sum(volData,1)) == prod(exp(volData))  )
%   **update_norm  - scalar, relative tomo volume change since between the original and update volData 



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

function [volData, update_norm] = update_volume(volData, projData_c_upd, update_step, projData_model, par)
    % update the volData using the provided projection update
    % projData_c_upd and update_step

    Npx_vol = size(volData); 
    
    position_offset = -projData_model.position_offset; 
    if any(any(round(position_offset) ~= position_offset))
        error('Noninterger shift may result in mixing real and imaginary data')
    end
    
    % apply integer shift to correct for shifts of the projections 
    projData_c_upd = utils.imshift_fast(projData_c_upd, position_offset(1), position_offset(2)); 

    % crop back to projection size
    projData_c_upd = utils.crop_pad(projData_c_upd,[Npx_vol(3),Npx_vol(1)]);   
   

    Nblock = ceil(prod(Npx_vol)*8 / 1e9);  % split into 2GB blocks 
    update_norm = 0; 

    for ii = 1:Nblock
        ind = 1+(ii-1)*ceil(Npx_vol(3)/Nblock):min((ii)*ceil(Npx_vol(3)/Nblock), Npx_vol(3)); 
        
        %% update the reconstructed volume 
        args = {projData_model.angle, [Npx_vol(1:2), length(ind)]};
        [volData_upd_r] = ptychotomo.back_proj(real(projData_c_upd(ind,:,:)),args{:});
        [volData_upd_i] = ptychotomo.back_proj(imag(projData_c_upd(ind,:,:)),args{:});

        volData_upd = complex(volData_upd_r, volData_upd_i); 
        update_norm = update_norm + gather(norm(volData_upd(:))); 
        
        if update_norm /  par.norm_full  > 0.5
            keyboard
        end
        % update only the processed block 
        if update_step > 0 && update_norm /  par.norm_full  < 0.5
           %% update the full volume, vary the update speed in dependence on number of layers -> avoid instabilities 
            volData_block = volData(:,:,ind); 
            volData_block = arrayfun(@object_constraint_kernel,volData_block, par.support_mask, update_step,  volData_upd, par.max_value, par.min_value);
            %% call MEX function 
            set_to_block_gpu(volData, volData_block, uint16(ind));        
        end
    
    
    end

    update_norm = update_norm / par.norm_full ; 



    %utils.verbose(0, 'update_norm %g', update_norm)
end

  
 %% auxiliary function for GPU, use GPU arrayfun for fast inplace processing 

function volData = object_constraint_kernel(volData, support_mask, update_step, volData_update, max_value, min_value)

    % separated update_step / converegnce rate for phase and absorbtion 
    volData_update = complex(real(update_step) * real(volData_update), imag(update_step) * imag(volData_update)); 
    volData = volData .* support_mask;

       
    volData = volData + volData_update; 

    volData_r = real(volData);
    volData_i = imag(volData);
    
    % apply some clipping of the output values 
    volData_r = min(real(max_value), max(real(min_value), volData_r));
    volData_i = min(imag(max_value), max(imag(min_value), volData_i));

    volData = complex(volData_r, volData_i) ; 


end
