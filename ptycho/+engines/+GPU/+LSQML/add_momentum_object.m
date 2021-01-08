% ADD_MOMENTUM_OBJECT speed up convergence by adding momentum / integral term into the
% update direction for the MLc method 
%
% [self, cache] = add_momentum_object(self, cache, par, object_upd_sum, iter, fourier_error) 
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** par       structure containing parameters for the engines 
% ** object_upd_sum     cell of arrays containing update direction from the LSQML method 
% ** iter               current iteation numebr 
% ** fourier_error      array [Npos,1] containing evolution of reconstruction error 
%
% returns:
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ cache     structure with precalculated values to avoid unnecessary overhead
%
%


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




function [self, cache] = add_momentum_object(self, cache, par, object_upd_sum, iter, fourier_error, beta_object)

   import utils.*
   import engines.GPU.GPU_wrapper.*
   import math.*

   verbose(1, 'Adding momentum to object')

               
   object_modes = size(self.object,1);

   % create velocity maps in the first iterations 
   if iter == par.object_change_start
        cache.object_upd_sum = cell(object_modes,1); 
        for ll = 1:object_modes
           for jj = 1:par.Nlayers
                cache.velocity_map_object{ll,jj} = Gzeros(self.Np_o);
                cache.object_upd_sum{ll,jj} = []; 
           end
        end
   end
   
   % store previous update directions 
   for ll = 1:object_modes
        for jj = 1:par.Nlayers
            upd = object_upd_sum{ll,jj}(cache.object_ROI{:})* mean(beta_object(:,jj)); 
            upd = upd / norm2(upd); 
            cache.object_upd_sum{ll,jj} = [cache.object_upd_sum{ll,jj}, {upd}]; 
        end
   end
   
   % how many steps are stored to calculate optimal friction 
   momentum_memory = 2; 

   % use Fourier error to avoid issues with convergence 
    ind_compare = find(all(~isnan(fourier_error),2),3,'last'); 
    if length(ind_compare) > 2
        merr = mean(fourier_error,2); 
        ferr_ok = max(merr(ind_compare([1,2]))) > min(merr(ind_compare([2,3]))); 
    else
        ferr_ok = true; 
    end
    
   if iter > momentum_memory+par.object_change_start
        for ll = 1:object_modes
            for jj = 1:par.Nlayers
               cache.object_upd_sum{ll,jj}(1) = []; % delete the oldest stored object update
               corr_level = nan; 
               
               if ferr_ok
                   % caculate correlation between updated to estimate optimal
                   % friction, !! NOTE that cache.object_upd_sum contains only the
                   % object_ROI region !!
                   % 
                   switch momentum_memory
                       case 2, [aux{1}, aux{2}] = compare_upd_directions_2(cache.object_upd_sum{ll,jj}{:});
                       case 3, [aux{1}, aux{2}, aux{3}] = compare_upd_directions_2(cache.object_upd_sum{ll,jj}{:}); 
                       othewise, error('Not implemented')
                   end
                   for kk = 1:momentum_memory
                       corr_level(kk) = real(Ggather(mean2(aux{kk}))); 
                   end
               end
               
               if ferr_ok && all(corr_level > 0 )


                    % estimate optimal friction from previous steps 
                    poly_fit = polyfit(0:momentum_memory,[0,log(corr_level)],1); 

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    gain = par.momentum;                           % smaller -> lower relative speed (less momentum)
                    friction =  0.5*max(-poly_fit(1),0);   % smaller -> longer memory, more momentum 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % update object by the momentum gradient 
                    [self.object{ll,jj}, cache.velocity_map_object{ll,jj}] = update_momentum(self.object{ll,jj}, cache.velocity_map_object{ll,jj}, object_upd_sum{ll,jj}, friction, gain, cache.illum_sum_0{ll}, cache.MAX_ILLUM(ll)); 
               else
                  % error was increased or updates are not well correlated, skip acceleration 
                   gain = 0; friction = inf; 
                   cache.velocity_map_object{ll,jj} =  cache.velocity_map_object{ll,jj}/2;
               end

               if verbose()> 1
                   norm_upd = norm2(object_upd_sum{ll,jj});
                   norm_vmap = norm2(gain.*cache.velocity_map_object{ll,jj});
                   verbose(2,['Object %i Corr=', repmat('1:%5.2f ',1,momentum_memory)],ll, corr_level )
                   verbose(0, 'Momentum: friction=%3.1e \tacceleration %1.1fx',friction,  1+norm_vmap/norm_upd) 
               end
           end
        end
       
        % plotting.smart_figure(25454)
        % img = cat(1,cat(2, cache.velocity_map_object{:}) .* cat(2,cache.illum_sum_0{:}),  cat(2, object_upd_sum{:}) .* cat(2,cache.illum_sum_0{:})); 
        % aimg = abs(img); 
        % img = min(aimg, quantile(aimg(:), 0.99)) .* img ./ (aimg+1e-3); 
        % plotting.imagesc3D(img(1:4:end, 1:4:end))
        % axis off image xy 
        % drawnow 
                
        try
           if verbose()> 2
               plotting.smart_figure(121)
               subplot(1,2,1)
               plotting.imagesc3D(cache.velocity_map_object{1,1})
               axis off image 
               title(sprintf('Velocity, iter=%i', iter))
               subplot(1,2,2)
               plotting.imagesc3D(object_upd_sum{1,1})
               axis off image 
               title('Gradient')
               drawnow 
           end
        catch
            keyboard
        end

   end
end

function  [object, Vmap] = update_momentum(object, Vmap, Vmap_upd, friction, gain, weight, w_max)
    % auxiliary function 
    weight = weight ./ (0.1*w_max+weight); 
    Vmap = (1-friction)*Vmap + Vmap_upd;
    object = object + weight.*gain.*Vmap; 
end

function  [out1, out2] = compare_upd_directions_2(upd1, upd2, upd3)
    % compare updates 
    upd3 = conj(upd3); 
    out2 = upd1 .* upd3; 
    out1 = upd2 .* upd3; 
end
function  [out1, out2, out3] = compare_upd_directions_3(upd1, upd2, upd3, upd4)
    % compare updates 
    upd4 = conj(upd4); 
    out3 = upd1 .* upd4; 
    out2 = upd2 .* upd4; 
    out1 = upd3 .* upd4; 
end
