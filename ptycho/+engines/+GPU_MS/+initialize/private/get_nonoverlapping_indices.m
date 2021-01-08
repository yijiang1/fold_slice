%  GET_NONOVERLAPPING_INDICES a heuristic based method to select pseudorandom indices of non overlapping regions   
% Note: It can be slow for large number of scanning positions 
%
% [indices_out, scan_ids_out] = get_nonoverlapping_indices(self, cache, par )      
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns: 
% ++ indices_out   cell of arrays, contain indices of positions processed in parallel 
% ++ scan_ids_out  cell of arrays, contain scan numbers for each position 


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


function [indices_out, scan_ids_out] = get_nonoverlapping_indices(self, cache, par )      

% find groups accross the scans in order to further minimize overlap
group_across_scans = true; %need to be true for sharing object amongs scans

if group_across_scans
% divide the grouping equally over all the scans 
    grouping = ceil(par.grouping/par.Nscans);
else
    grouping = par.grouping;
end
max_groups = 0;

for kk = 1:par.Nscans
        
    % setdiff sort the indices by size 
    indices_0 = self.reconstruct_ind{kk};  % remove unwanted from the decision process 
    ind_start(kk) = min(indices_0)-1;
    indices_0 = indices_0 - ind_start(kk);  % 
    Npos_tmp=length(indices_0);
    % randomly permutate the indices    
    indices_0 = indices_0(randperm(Npos_tmp)); 
    max_groups = max(max_groups, ceil(Npos_tmp/grouping));
    % fill it with some initial random guess 
    for ii = 1:ceil(Npos_tmp/grouping) 
         indices{kk}{ii} = indices_0(1+(ii-1)*grouping : min(Npos_tmp,ii*grouping));
    end
   
    % no need for this method ot it calculation would be too long -> use
    % just the random initial guess 
    if (self.Npos/par.Nscans > 1e3 ) || (grouping == 1) || ~isfield(cache, 'distances_matrix')
        %%%for ii = 1:length(indices{1}) %why length(indices{1})? Bug?
        for ii = 1:length(indices{kk}) %modified by YJ to prevent error when different scans have differernt number of positions
            scan_ids{kk}{ii} = ones(1,length(indices{kk}{ii}))*kk;  % note their scan origin
        end
        continue 
    end  % hope that for large number of positions the random statistics will be enough
    
    try

     update_score = 0;
     for i = 1:ceil(Npos_tmp/grouping)-1
         id = indices{kk}{i};
         dist_mat_small = cache.distances_matrix{kk}(id,id);
         for ii = 0:2*length(indices{kk}{i+1})  % go twice through all positions 
             j = 1+mod(ii, length(indices{kk}{i+1}));
             min_dist = 1./sum(1./dist_mat_small.^2); % find the shortest distance between the probes
             if all(isinf(min_dist)) % all(isnan(min_dist))
                 break
             end
             [~,min_dist_ind] = min(min_dist);
            % make a swap with the j position in i+1 index array
             tmp = indices{kk}{i+1}(j);
             indices{kk}{i+1}(j) = indices{kk}{i}(min_dist_ind);
             indices{kk}{i}(min_dist_ind) = tmp;
             
             
            % update distance matrix 
            dist_mat_small_update = cache.distances_matrix{kk}(tmp,indices{kk}{i});
            dist_mat_small(min_dist_ind,:) = dist_mat_small_update';
            dist_mat_small(:,min_dist_ind) = dist_mat_small_update;
         end
         update_score = update_score +j;
     end
    catch
        keyboard
    end
        

     % fill the last group by the skip indieces but do not expand it 
     skip_ind = cache.skip_ind(randperm(length(cache.skip_ind)));
     indices{kk}{end} = [indices{kk}{end}, skip_ind(1:min(end, grouping-length(indices{kk}{end})))]; % join skip_ind back to the last (smallest) set    
          
     for ii = 1:length(indices{kk})
        scan_ids{kk}{ii} = ones(1,length(indices{kk}{ii}))*kk;  % note their scan origin
     end
end


if group_across_scans
    indices_out = cell(max_groups,1);
    scan_ids_out = cell(max_groups,1);
    %% merge groups from difference scans into larger chunks if required 
    for ii = 1:max_groups
        indices_out{ii} = [];
        scan_ids_out{ii} = [];
        % from each scan add one group
        for kk = 1:par.Nscans
            if ii <= length(indices{kk})
                indices_out{ii} = [indices_out{ii}, indices{kk}{ii}+ind_start(kk)];
                scan_ids_out{ii} = [scan_ids_out{ii}, scan_ids{kk}{ii}];
            end
        end
        if length(scan_ids_out) > 1 && length(scan_ids_out{end}) < grouping / 10
            % if the a group is too small, merge it with the previous to
            % reduce the overhead 
            indices_out{end-1} = [indices_out{end-1}, indices_out{end}]; 
            scan_ids_out{end-1} = [scan_ids_out{end-1}, scan_ids_out{end}];
            scan_ids_out(end) = []; indices_out(end) = []; 
        end
    end
else
    indices_out = {};
    for ii = 1:par.Nscans
        for kk = 1:length(indices{ii})
            indices_out = [indices_out, indices{ii}{kk}+ind_start(ii)];
        end
    end
    scan_ids_out = [scan_ids{:}]';

end

%% sort them to minimize allocation of new projection matrices
Nitems = cellfun(@length, indices_out); 

if all(max(Nitems) - min(Nitems) <= 1) && all(Nitems > 100)
    % just neglect one scanning position to keep the bunches with the same
    % size -> faster run on GPU 
    for i = 1:length(indices_out)
        indices_out{i} = indices_out{i}(1:min(Nitems));
        scan_ids_out{i} = scan_ids_out{i}(1:min(Nitems));
    end
else
    [~,ind] = sort(Nitems(:),1,'descend' );
    indices_out = indices_out(ind);
    scan_ids_out = scan_ids_out(ind);
end


    
end



