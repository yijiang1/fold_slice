% GET_PARALLEL_BLOCKS Find the optimal groups to be solved in parallel on GPU/CPU
% 
%[cache, par] = get_parallel_blocks(self, par, cache)
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ par       structure containing parameters for the engines 



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


function [cache, par] = get_parallel_blocks(self, par, cache)
    import engines.GPU_MS.GPU_wrapper.*
    import utils.*
    import engines.GPU_MS.shared.*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIND MAXIMAL GROUP SIZE IF GPU IS USED %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grouping_0 = par.grouping; 

    par.grouping = round(min(self.Npos, par.grouping));

    global gpu 

    if par.use_gpu
        %% !! very empirical estimation of the GPU memory requirements !!!
        % use the estimated memory requirements to prevent low memory issues 

        [required_mem_constant] = estimate_req_memory(self, par, 0); 
        [required_mem] = estimate_req_memory(self, par, 1); 
        % precheck size of the block and try to optimize size of the groups
        % allowed for given GPU 
        max_group_size = floor( (gpu.AvailableMemory - required_mem_constant) ./ (required_mem - required_mem_constant));
        max_group_size = min(self.Npos, max_group_size); 
        verbose(1,'Maximal possible grouping %i', max_group_size);

        % if group size was set to infinity, is the maximal group size possible 
        if isinf(grouping_0)
            par.grouping= max_group_size; 
        else
            % otherwise use max_group_size as a top limit 
            par.grouping = min(par.grouping, max_group_size); 
        end
        
        % adjust grouping to minimize overhead -> make the group sizes more
        % equal 
        if is_method(par, {'ML', 'PIE'})
            % allows to calculate several scans together 
            par.grouping = ceil(self.Npos/ceil(self.Npos/par.grouping));
        else
            % consider each scan separatelly
            Npos_scan = cellfun(@length, self.reconstruct_ind); 
            par.grouping = max(ceil(Npos_scan./ceil(Npos_scan./par.grouping)));
        end

        if par.grouping ~= grouping_0 
            verbose(1,'Optimal grouping was changed from %i to %i ', grouping_0, par.grouping);
        end
        if par.grouping < 1
            error('Too low memory, use smaller dataset or try ePIE')
        end

        verbose(1,'Selected grouping %i', par.grouping);

    else
        if is_method(par, {'DM', 'ML'})
            par.grouping = self.Npos;
        end
    end

    % precalculate distance matrix for pseudo ePIE / hPIE / MLs to get
    % least overlapping indices 
    if is_method(par, {'ML', 'PIE'})
        if self.Npos/par.Nscans < 1e3 %added by YJ to save memory
            for ll = 1:par.Nscans
                dist_mat =  single(distmat(self.probe_positions_0(self.reconstruct_ind{ll},:))); 
                dist_mat(dist_mat==0 |  dist_mat > max(self.Np_p)/2) = inf;
                cache.distances_matrix{ll} = dist_mat;
            end
        end
    end
    if is_method(par, 'MLc')
        % get higly overlapping subsets of indices for PIE / ML
         [cache.preloaded_indices_compact{1}.indices,cache.preloaded_indices_compact{1}.scan_ids] = ...
             get_close_indices(self, cache, par );     
    elseif is_method(par, {'MLs', 'PIE'})
        % preload order of indices , generate several of them to add randomness
        for i = 1:min(par.number_iterations,10)  
            [cache.preloaded_indices_sparse{i}.indices,cache.preloaded_indices_sparse{i}.scan_ids] = ...
                get_nonoverlapping_indices(self, cache, par );
        end
    end
    % get just some predefined sets of indices - RAAR, DM , !! order
    % does not matter 
    [cache.preloaded_indices_simple{1}.indices,cache.preloaded_indices_simple{1}.scan_ids] = ...
         get_scanning_indices(self, cache, par );
     

end
