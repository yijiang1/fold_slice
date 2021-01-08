% UPDATE_PROBE calculate improved update direction
% apply "overlap" constraint to get better estimate  of the update direction
% 
%  probe = update_probe(probe, m_probe_update, par, p_ind, g_ind, beta_probe, Nind)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** object_update_proj  [Nx, Ny, N] array, estimate of the object update for each scan position, ie conj(P)*chi
% ** object_upd_sum      cell of object sized arrays containg previous optimal updates 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** beta_probe  (scalar) relaxation parameter of the update step 
% ** Nind        (int) number of groups that are solved serially

% returns:
% ++ probe      cell of the updated probes 
%
%
% see also: engines.GPU_MS.LSQML 


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



function  probe = update_probe(probe, m_probe_update, par, p_ind, g_ind, beta_probe, Nind)
    import engines.GPU_MS.shared.*
    import engines.GPU_MS.GPU_wrapper.*
    import utils.verbose 
    
    %% update probe 
    beta_probe = beta_probe(g_ind);
    probe_ids = unique(p_ind); 
    if is_method(par, 'MLc')
        beta_probe = beta_probe / Nind; % in order to make the compact version closer to original ML method, the accumulated probe step per iteration should be 1 
    end
    
    if (par.share_probe || par.Nscans == 1) && size(probe,3)==1
        % most simple case, no multiprobe needed 
        probe = probe + m_probe_update .* mean(beta_probe);
    elseif length(probe_ids) == 1
        % variable probe extension with shared probe 
        probe(:,:,probe_ids,1) = probe(:,:,probe_ids,1) + m_probe_update .* mean(beta_probe);
    else % unshared probe 
        % update each of the probes separately 
        for id = probe_ids(:)'
            ind = p_ind == id;
            beta = mean(beta_probe(ind));
            probe(:,:,id,1) = Gfun(@upd_probe_Gfun, probe(:,:,id,1),m_probe_update(:,:,min(end,id)), beta);
        end     
    end
        
    if verbose()> 3
        % show applied subsets (update amplitude) and probe update
        % amplitude 
        plotting.smart_figure(11)
        probe_modes = size(m_probe_update,3); 
        for ll = 1:probe_modes
            subplot(probe_modes,2,2+probe_modes*(ll-1))
            p = fftshift(fft2(fftshift(m_probe_update(:,:,ll)))); 
            p = min(abs(p), quantile(abs(p(:)), 0.999)) .* p ./ abs(p); 
            plotting.imagesc3D(p);
            axis off image xy 
        end
        title('Probe update')
        colormap bone 
        drawnow 
    end
end

function probe = upd_probe_Gfun(probe,probe_update, alpha_p)
    probe =  probe + alpha_p.*probe_update;
end
