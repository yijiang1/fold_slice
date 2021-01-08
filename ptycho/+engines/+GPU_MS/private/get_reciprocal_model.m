% GET_RECIPROCAL_MODEL get estimate of the measured intensity from given electric field Psi
% 
% [aPsi, aPsi2, cache, self] = get_reciprocal_model(self, Psi, modF, mask,iter, g_ind, par, cache)
%
% ** self           structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** Psi            where Psi is the propagated exitwave   
% ** modF           pre-fftshifted and sqrt-ed data
% ** mask           masked values, 1 = ignored, 0 = use this pixel 
% ** iter           current iteration number 
% ** ind            processed indices 
% ** par            structure containing parameters for the engines 
% ** cache          structure with precalculated values 
% 
% returns:
% ++ aPsi           reciprocal amplitude model 
% ++ aPsi           reciprocal intensity model 
% ++ cache          structure with precalculated values 
% ++ self           structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..



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

function [aPsi, aPsi2, cache, self] = get_reciprocal_model(self, Psi, modF, mask,iter, g_ind, par, cache)
    import engines.GPU_MS.GPU_wrapper.*
    aPsi2 = [];

    if  par.probe_modes == 1 && isempty(self.background) && self.diff_pattern_blur == 0 ...
            && ~par.background_detection && strcmpi(par.likelihood, 'l1') &&  par.upsampling_data_factor == 0
        % or the simplest and fastest option: just get absolute value
        aPsi = abs(Psi{1});
        aPsi2 = [];
    elseif par.probe_modes == 1 && ~isempty(self.background) && self.diff_pattern_blur == 0 ...
            && ~par.background_detection && strcmpi(par.likelihood, 'l1')  &&  par.upsampling_data_factor == 0
        % second simplest option, abs + background
        aPsi = Gfun(@modulus_with_background,Psi{1}, self.background , cache.background_profile);
    else
        % apply corrected model and sum up all coherence modes
        aPsi2 = sumsq_cell(Psi);

        % assume that the data were upsampled by the utils.unbinning_2D function 
        if par.upsampling_data_factor
            aPsi2 = utils.binning_2D(aPsi2,  2^par.upsampling_data_factor);
        end

        %%%%%%%%%%%%%%%%  linear correction model %%%%%%%%%%%%%%%%%%%%%%%
        [aPsi2,cache, self] = get_linear_correction_model(self,par,cache,aPsi2,modF,mask,iter, g_ind );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aPsi = sqrt(aPsi2);
    end
end

function aPsi = modulus_with_background(Psi, background_value, background_shape)
    rPsi = real(Psi); 
    iPsi = imag(Psi); 
    aPsi2 = rPsi.^2 + iPsi.^2;    
    aPsi2 = aPsi2 + background_value .* background_shape;
    % sqrt is very slow ... 
    aPsi = sqrt(aPsi2);
    % aPsi = exp(0.5*log(aPsi2)); % log identity has exactly the same calculation time  
end

function y = sumsq_cell(x)
% Description: sum incoherently cells x, make it inplace and fast 

    N = length(x);
    if N <= 15 && builtin( 'isa', x{1}, 'gpuArray' )
       switch N
           case 1, fun = @sum_1; 
           case 2, fun = @sum_2; 
           case 3, fun = @sum_3; 
           case 4, fun = @sum_4; 
           case 5, fun = @sum_5; 
           case 6, fun = @sum_6; 
           case 7, fun = @sum_7; 
           case 8, fun = @sum_8; 
           case 9, fun = @sum_9; 
           case 10, fun = @sum_10; 
           case 11, fun = @sum_11; 
           case 12, fun = @sum_12; 
           case 13, fun = @sum_13; 
           case 14, fun = @sum_14; 
           case 15, fun = @sum_15; 
       end
        y = arrayfun(fun, x{:}); 
    else
        y = 0;
        for i = 1:N
             y = y + abs(x{i}).^2;
        end
    end
end

% !! using sqrt(imag(x)^2 + real(x)^2) is much slower !!!

% merged GPU kernels 
function y = sum_1(x)
    y = abs(x).^2;
end
function y = sum_2(x1,x2)
    y = abs(x1).^2+abs(x2).^2;
end
function y = sum_3(x1,x2,x3)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2;
end
function y = sum_4(x1,x2,x3,x4)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2;
end
function y = sum_5(x1,x2,x3,x4,x5)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2;
end
function y = sum_6(x1,x2,x3,x4,x5,x6)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2;
end
function y = sum_7(x1,x2,x3,x4,x5,x6,x7)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2;
end
function y = sum_8(x1,x2,x3,x4,x5,x6,x7,x8)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2;
end
function y = sum_9(x1,x2,x3,x4,x5,x6,x7,x8,x9)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2;
end
function y = sum_10(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2;
end
function y = sum_11(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2+abs(x11).^2;
end
function y = sum_12(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2+abs(x11).^2+abs(x12).^2;
end
function y = sum_13(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2+abs(x11).^2+abs(x12).^2+abs(x13).^2;
end
function y = sum_14(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2+abs(x11).^2+abs(x12).^2+abs(x13).^2+abs(x14).^2;
end
function y = sum_15(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
    y = abs(x1).^2+abs(x2).^2+abs(x3).^2+abs(x4).^2+abs(x5).^2+abs(x6).^2+abs(x7).^2+abs(x8).^2+abs(x9).^2+abs(x10).^2+abs(x11).^2+abs(x12).^2+abs(x13).^2+abs(x14).^2+abs(x15).^2;
end
function [aPsi2, cache, self] = get_linear_correction_model(self,par,cache,aPsi2,modF,mask, iter, ii )
       import engines.GPU_MS.GPU_wrapper.*
 
        if isempty(self.background) && self.diff_pattern_blur == 0 && strcmp(par.background_detection, 'none') 
            return  % nothing to be done, return 
        end
        
        %% add background 
        if ~isempty(self.background)  
            if ~isfield(cache, 'background_profile') || isscalar(cache.background_profile)
                aPsi2 = aPsi2 + self.background; 
            else
                aPsi2 = Gfun(@add_background, aPsi2, self.background, cache.background_profile,modF);
            end
        end
        
        if self.diff_pattern_blur > 0 
            %%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL CORRECTIONS START  %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%       
            if isempty(self.modes{1}.ASM_factor)  % is not nearfield
                aPsi2 = aPsi2(cache.fftshift_idx{:},:);
            end

            % apply blur correction 
            if self.diff_pattern_blur > 0
                % generate blurring kernel 
                x = [-1, 0,-1]/self.diff_pattern_blur;
                [X,Y] = meshgrid(x,x);
                blur_kernel = exp( -(X.^2 + Y.^2) );
                blur_kernel  = blur_kernel  / math.sum2( blur_kernel ); 
                aPsi2 = convn(aPsi2, blur_kernel, 'same');  
            end

            if isempty(self.modes{1}.ASM_factor)  % is not nearfield
                aPsi2 = aPsi2(cache.fftshift_idx{:},:);
            end
        end
        
        %  simple estimation of background 
        if par.background_detection && iter > par.background_detection
                if isempty(mask); mask = false; end
                % calculate the most optimal background update 
                [nom,denom] = Gfun(@get_background_estimate,modF, aPsi2, mask,  cache.background_profile_weight, cache.background_profile );
                update = sum2(nom)./sum2(denom); 
                if any(ii == 1)
                    fprintf('Background update: %3.3g curr value:%3.3g \n ', mean(update), self.background);
                end   
                
                
                
                self.background = posit(self.background + (par.grouping/self.Npos)*mean(update));               
         
%           %% Check if background is fitted well 
%             X = (-self.Np_p(1)/2:self.Np_p(1)/2-1);
%             Y = (-self.Np_p(2)/2:self.Np_p(2)/2-1);
%             [X,Y] = meshgrid(X,Y);
%  
%             R = (sqrt(X.^2 + Y.^2));
%             D = fftshift(single(modF.^2) - aPsi2);
%             for i = 1:mean(self.Np_p)/2
%                 progressbar(i, mean(self.Np_p)/2);
%                 mask = (R==i);
%                 mask = mask / sum2(mask);
%                 B(i) = median(sum2(bsxfun(@times, D, mask)));
%             end
%             plot(B)
%             ylim([-5,5])
%             drawnow 
%             


        end
end

function aPsi2 = add_background(aPsi2, background,background_profile,modF)
       aPsi2 = aPsi2 + background .* background_profile .* (modF > 0);  % leave empty pixels empty 
end

function [nom, denom] = get_background_estimate(modF, aPsi2, mask, distribution, background )
       W = ~mask .* distribution;
       nom =  W.* (modF.^2 - aPsi2).*background; 
       denom = W.*background.^2;
end





