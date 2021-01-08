% MODULUS_CONSTRAINT universal fast relaxed modulus constraint for amplitude and poission likelihood
%
% [chi,R] = modulus_constraint(modF,aPsi ,Psi, mask, noise, relax_noise,likelihood , R_offset)
%
% ** modF           pre-fftshifted and sqrt-ed data
% ** aPsi           reciprocal amplitude model 
% ** Psi            where Psi is the propagated exitwave   
% ** mask           masked values, 1 = ignored, 0 = use this pixel 
% ** noise          estimated noise (STD) in each pixel after sqrt transform 
% ** relax_noise    multiplicative constant for the providede noise value 
% ** likelihood     L1 or poisson 
% ** R_offset       number to be subtracted from the modF / aPsi ratio 
% 
% returns:
% ++ chi        the updated wavefront or wavefront update, depends on the chosen R_offset
% ++ R          modF / aPsi ratio 

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
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
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
% 
%   

function [chi,R] = modulus_constraint(modF,aPsi ,Psi, mask, noise, par, R_offset)

    import engines.GPU.GPU_wrapper.*

    likelihood = lower(par.likelihood);
    relax_noise = par.relax_noise; 
    
    Nmodes = length(Psi); 
    R = [];
   
    if par.upsampling_data_factor      
        % calculate the convolution with delta functions 
        modF = utils.unbinning_2D(modF, 2^par.upsampling_data_factor );  
        aPsi = utils.unbinning_2D(aPsi, 2^par.upsampling_data_factor ); 
        if ~isempty(mask) && ~isscalar(mask)
           mask = utils.unbinning_2D(mask, 2^par.upsampling_data_factor );  
        end
    end
            
        
    %% write the most common cases as merged kernels 
    if nargout == 1
        if isempty(mask) && isempty(noise) && relax_noise == 0 && strcmpi(likelihood, 'l1') && Nmodes == 1  % common modulus constraint 
            chi{1} = Gfun(@modulus_non_relaxed,Psi{1},modF, aPsi, R_offset); 
            return
        end
        if ~isempty(mask) && isempty(noise) && strcmpi(likelihood, 'l1') && Nmodes == 1  % common modulus constraint 
            chi{1} = Gfun(@modulus_weight_relaxed,Psi{1},modF, aPsi,mask, R_offset); 
            return
        end
    end
    
    if isempty(mask) && isempty(noise) && relax_noise == 0  % common modulus constraint 
        switch likelihood
            case 'l1', R = Gfun(@non_relaxed, modF, aPsi, R_offset); 
            case 'poisson', R = Gfun(@poisson_noise_relaxed,modF, aPsi, R_offset);
        end
    elseif isempty(mask) && isempty(noise) && relax_noise > 0  % common modulus constraint 
        switch likelihood
            case 'l1', R = Gfun(@weight_relaxed, modF, aPsi, relax_noise, R_offset ); 
            case 'poisson', R = Gfun(@poisson_noise_weight_relaxed,modF, aPsi,relax_noise, R_offset); 
        end
    elseif  ~isempty(mask) && isempty(noise)
        switch likelihood
            case 'l1', R = Gfun(@weight_relaxed, modF, aPsi, max(mask, relax_noise), R_offset ); 
            case 'poisson', R = Gfun(@poisson_noise_weight_relaxed,modF, aPsi, max(mask, relax_noise), R_offset); 
        end
    elseif isempty(mask) && ~isempty(noise)
        R = Gfun(@noise_relaxed, modF, aPsi, noise, relax_noise, R_offset); 
    else
        R = Gfun(@noise_weight_relaxed, modF, aPsi, noise,relax_noise, mask, R_offset); 
    end
    
    for i = 1:Nmodes
        chi{i} = Psi{i} .* R;  % apply the constraint to the currect estimation 
    end
            
end

% classical modulus 
function chi = modulus_non_relaxed(Psi,modF, aPsi, R_offset)
    R = (modF./(aPsi+1e-9)- R_offset);
    chi = R .* Psi;
end
function chi = modulus_weight_relaxed(Psi,modF, aPsi,W, R_offset)
    R = (W+(1-W).*modF./(aPsi+1e-9))- R_offset; 
    chi = R .* Psi;
end
function R = non_relaxed(modF, aPsi, R_offset)
      R = modF./(aPsi+1e-9)- R_offset;
end

% modulus for noisy data with relaxation 

function R = weight_relaxed(modF, aPsi, W, R_offset)
      R = (W+(1-W).*modF./(aPsi+1e-9))- R_offset;
end
function R = noise_relaxed(modF, aPsi, noise, relax_noise, R_offset)
      W = 1 ./ (1+ ((aPsi - modF)./(noise .* relax_noise)).^2 ); 
      R = (W+(1-W).*modF./(aPsi+1e-9))- R_offset;
end
function R = noise_weight_relaxed(modF, aPsi, noise,relax_noise, W0, R_offset)
    % if mask == 1 then W has to be 0 
    W = 1-(1-W0) .* (1-1./ (1+ ((aPsi - modF)./(noise .* relax_noise) ).^2 )); 
    R = (W+(1-W).*modF./(aPsi+1e-9))- R_offset;
end
function R = poisson_noise_weight_relaxed(modF, aPsi, W, R_offset)
    % if mask == 1 then W has to be 0 
    % Maximum-likelihood refinement for coherent diffractive imaging
    R = (W+(1-W).*modF.^2 ./(aPsi.^2+1e-3))- R_offset;
end
function [R] = poisson_noise_relaxed(modF, aPsi, R_offset)
    % if mask == 1 then W has to be 0 
    % Maximum-likelihood refinement for coherent diffractive imaging
    R = modF.^2 ./(aPsi.^2+1e-3)- R_offset;
end
