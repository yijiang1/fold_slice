% ADD_MOMENTUM_PROBE speed up convergence by adding momentum / integral term into the
% update direction for the MLc method 
%
% [self, cache] = add_momentum_probe(self, cache, par, probe_upd, iter, fourier_error) 
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** par       structure containing parameters for the engines 
% ** probe_upd          cell of arrays containing update direction from the LSQML method 
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




function [self, cache] = add_momentum_probe(self, cache, par, probe_upd, iter, fourier_error, beta_probe)

   import utils.*
   import engines.GPU.GPU_wrapper.*
   import math.*

   verbose(1, 'Adding momentum to probe')

   %disp(size(self.probe))
       
   Nmodes = size(self.probe,1);

   % create velocity maps in the first iterations 
   if iter == par.probe_change_start
        cache.probe_upd_sum = cell(Nmodes,1); 
        for ll = 1:Nmodes
            cache.velocity_map_probe{ll} = Gzeros(self.Np_p);
            cache.probe_upd_sum{ll} = []; 
        end
   end
   
   % store previous update directions 
   for ll = 1:Nmodes
        probe_upd{ll} = probe_upd{ll}(:,:,:,1); % accelerate only the fundamental probe in the variable probe extension  
        probe_upd{ll} = probe_upd{ll} .* mean(beta_probe(:)); 
        upd = probe_upd{ll} ./ (norm2(probe_upd{ll})+eps); 
        cache.probe_upd_sum{ll} = [cache.probe_upd_sum{ll}, {upd}]; 
   end
   
   % how many steps are stored to calculate optimal friction 
   momentum_memory = 3; 

   if iter > momentum_memory+par.probe_change_start
        for ll = 1:Nmodes
           probe_modes= size(self.probe{ll},3);
           cache.probe_upd_sum{ll}(1) = []; % delete the oldest stored probe update
           % caculate correlation between updated to estimate optimal
           % friction, 
           
           [aux{1}, aux{2}, aux{3}] = compare_upd_directions(cache.probe_upd_sum{ll}{:}); 
           for kk = 1:momentum_memory
               corr_level(kk,:) = real(Ggather(mean2(aux{kk}))); 
           end


            % use Fourier error to avoid issues with convergence 
            ind_compare = find(all(~isnan(fourier_error),2),3,'last'); 
            if length(ind_compare) > 2
                merr = mean(fourier_error,2); 
                ferr_ok = max(merr(ind_compare([1,2]))) > min(merr(ind_compare([2,3]))); 
            else
                ferr_ok = true; 
            end
                        
            if all(corr_level(:) > 0 ) && ferr_ok

try
                % estimate optimal friction from previous steps 
                poly_fit = polyfit(repmat(0:momentum_memory,probe_modes,1),[zeros(1,probe_modes);log(corr_level)]',1); 
catch
    keyboard
end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                gain = par.momentum;                           % smaller -> lower relative speed (less momentum)
                friction =  0.5*max(-poly_fit(1),0);   % smaller -> longer memory, more momentum 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % update probe by the momentum gradient ,  accelerate only the fundamental probe in the variable probe extension  
                [self.probe{ll}(:,:,:,1), cache.velocity_map_probe{ll}] = update_momentum(self.probe{ll}(:,:,:,1), cache.velocity_map_probe{ll}, probe_upd{ll}, friction, gain); 
           else
              % error was increased or updates are not well correlated, skip acceleration 
               gain = 0; friction = inf; 
               cache.velocity_map_probe{ll} =  cache.velocity_map_probe{ll}/2;
           end

           if verbose()  > 0
               norm_upd = mean(norm2(probe_upd{ll}));
               norm_vmap = eps+mean(norm2(gain.*cache.velocity_map_probe{ll}));
               verbose(2,['Probe %i Corr=', repmat('1:%5.2f ',1,momentum_memory)],ll, corr_level(:,1) )
               verbose(0, 'Momentum: friction=%3.1e \tacceleration %1.1fx',friction(:,1),  1+norm_vmap./norm_upd(1)) 
           end
        end
       
        % plotting.smart_figure(25454)
        % img = cat(1,cat(2, cache.velocity_map_probe{:}) .* cat(2,cache.illum_sum_0{:}),  cat(2, probe_upd_sum{:}) .* cat(2,cache.illum_sum_0{:})); 
        % aimg = abs(img); 
        % img = min(aimg, quantile(aimg(:), 0.99)) .* img ./ (aimg+1e-3); 
        % plotting.imagesc3D(img(1:4:end, 1:4:end))
        % axis off image xy 
        % drawnow 
                
        try
           if verbose()  > 2
               plotting.smart_figure(1213)
               subplot(1,2,1)
               plotting.imagesc3D(cache.velocity_map_probe{1,1})
               axis off image 
               title(sprintf('Velocity, iter=%i, norm=%g', iter, norm2(cache.velocity_map_probe{1,1})))
               subplot(1,2,2)
               plotting.imagesc3D(probe_upd{1})
               axis off image 
               title(sprintf('Gradient, norm=%g',norm2(probe_upd{1})))
               drawnow 
           end
        catch
            keyboard
        end

   end
end

function  [probe, Vmap] = update_momentum(probe, Vmap, Vmap_upd, friction, gain)
    % auxiliary function 
    Vmap = (1-friction)*Vmap + Vmap_upd;
    probe = probe + gain.*Vmap; 
end
function  [out1, out2, out3] = compare_upd_directions(upd1, upd2, upd3, upd4)
    % compare updates 
    upd4 = conj(upd4); 
    out3 = upd1 .* upd4; 
    out2 = upd2 .* upd4; 
    out1 = upd3 .* upd4; 
end
