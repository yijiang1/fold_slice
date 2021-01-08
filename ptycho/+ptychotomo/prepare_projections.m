%PREPARE_PROJECTIONS, use ptycho reconstruction object_0 and initial guess
%object_guess_c to calculate updated linearized complex object_c, i.e. object_0 = exp(object_c)
% -> projections needs to be unwrapped 
%
% object_c = prepare_projections(object_0, Npx_proj, probe_size, merge_layers, object_guess_c, weight)
% 
%  Inputs: 
%   **object_0        - stack of 2D reconstructed layers from ptychography
%   **Npx_proj        - size of the tomography projections 
%   **probe_size      - p.asize value 
%   **merge_layers    - bool, if true, merge all layers to one before unwrapping
%   **object_guess_c  - stack of 2D layers used as initial guess for from ptychography
%   **weight          - array denoting reliability for different regions of object_0 array 
%   
%  *returns*
%   ++object_c        - linearized complex  projection, ie object_0 = exp(object_c)


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


function object_c = prepare_projections(object_0, Npx_proj, probe_size, merge_layers, object_guess_c, weight)
        % get updated unwrapped object_c given the new object estiamte and the original unwrapped
        % object before the ptycho iteration 

        if nargin < 4
            merge_layers = false;
        end
        
        %% smoothly merge with the last reconstruction at this angle
        object_0 = utils.crop_pad(object_0,Npx_proj(1:2));
        
        if nargin > 4 && ~isempty(object_guess_c)
            % get only difference from the previous and then use this for
            % unwrapping 
            object_guess = exp(object_guess_c); 
            object = (object_0 .* conj(object_guess)) ./ (1e-4+abs(object_guess).^2);
            clear object_guess
        else
            object = object_0; 
        end


        %% multilayer expantion
        
        % find reliability region 
        if merge_layers
            object = prod(object,3); 
        end
        
        win = tukeywin(Npx_proj(2), probe_size(2)/Npx_proj(2)/2)'.*tukeywin(Npx_proj(1), probe_size(1)/Npx_proj(1)/2);

        %% empirical way to estimate the unwrapping region 
        weight = utils.Garray(weight .* win); 
        if ~isempty(object_guess_c)
            Niter_unwrap = 3; % assume that initial guess is good enough 
        else
            Niter_unwrap = 10; % try more iterations in the first step 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        unwrapped_phase = tomo.unwrap2D_fft2_split(object ,[],0,  weight, [], [], [], Niter_unwrap); 
        attenuation = log(max(1e-1,abs(object_0)));

        if ~isempty(object_guess_c)
            % make the fit exact, maybe not needed , useful when the update changes are small 
            unwrapped_phase = unwrapped_phase + angle(object .* exp(-1i*unwrapped_phase)) .* weight;
      
            
            % limit maximal absorbtion to 90% and tranmission to 100% 
            object_c = complex(attenuation, imag(object_guess_c)+unwrapped_phase); 
        else
            
            object_c = complex(attenuation , unwrapped_phase); 
        end

        % zero the empty regions 
        object_c = object_c .* (object_0 ~= 0); 
        
        
        if merge_layers
            object_c = repmat(object_c,1,1,Npx_proj(3))/Npx_proj(3); 
        end
            
end
