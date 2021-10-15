% GET_CLOSE_INDICES simple based method to select indices for DM 
% !! GPU needs the sets to be with similar , ideally the same sizes !!! 
%
% [indices_out, scan_ids_out] = get_close_indices(self, cache, par )          
% 
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



function [indices_out, scan_ids_out] = get_close_indices(self, cache, par )      

import math.*
import utils.*

grouping = par.grouping;

% in case of a shared scan join together all positions to find the optimal groups 
group_across_scans = true; 


if par.share_object && group_across_scans
    Nsets = 1; 
else
    Nsets = par.Nscans; 
end
cluster_refinement_time = 0;
cluster_time = 0; 

% in simplest case process all positions together 
if Nsets == 1 && grouping >= self.Npos
    indices_out = {[self.reconstruct_ind{:}]}; 
    scan_ids_out{1} = []; 
    for ii = 1:length(self.reconstruct_ind)
        scan_ids_out{1} = [scan_ids_out{1}; ii*ones(length(self.reconstruct_ind{ii}),1)];
    end
    return
end

%rng default

for kk = 1:Nsets
    % take them sequentially but with random offset 
    if par.share_object  && group_across_scans
        % join all indices into one large set if the object is shared 
        indices_0 = [self.reconstruct_ind{:}];
        for ii = 1:length(self.reconstruct_ind)
            scans_0(self.reconstruct_ind{ii}) = ii;
        end
    else
        indices_0 = self.reconstruct_ind{kk};
        scans_0 = kk * ones(size(indices_0)); % scan number 
    end
    N = length(indices_0);
    Ngroups=ceil(N/grouping);
    positions = self.probe_positions_0(indices_0,:);
    Npos = length(positions);

    % get initial set  distribution 
    [groups, C, sum_D, D] = get_best_kmeans(positions, Ngroups); 

    iter = 0;
    t0 = tic;
    
    while true 
        iter= iter +1;
        [nbins,bins] = hist(groups, unique(groups));
        % if less than 2 types of groups are present, finish 
        Ngroups_sizes = length(unique(nbins));
        % try to find distribution with most similar sets sizes, if not
        % easy, end with suboptimal distribution  after 50 iterations 
        if ( Ngroups_sizes <= max(2, ceil(iter/1e3)) && (Ngroups*grouping ~= N || iter > 1e3 )) ...
                ||  Ngroups_sizes == 1 % choose suboptimal solution if better is not found soon
            break
        end

        % find group with lowest number of members , add new points into
        % this group 
        min_group = bins(argmin(nbins));
        large_groups = bins(nbins>grouping);
        if isempty(large_groups) || any(ismember(min_group, large_groups)) ; break; end 
        % choose closest position from the largest group to be moved to the
        % smallest group 
        ind_large = (D(:,min_group) == min(D(ismember(groups, large_groups), min_group)));

        groups(ind_large) = min_group;
       
    end

    % remove empty groups 
    ugroups = unique(groups); 
    Ngroups = length(ugroups); 
    groups = sum((1:Ngroups) .*(groups == ugroups'),2);
    
    
    
    cluster_time = cluster_time + toc(t0); 

    
    for ii = 1:Ngroups
        C(ii,:) = median(positions(groups == ii,:));
    end
    for ii = 1:Ngroups
        D(:,ii) = (sum((positions - C(ii,:)).^2,2));
    end

        
    t0 = tic;
    %% find more compact refinement 
    % find the most distanced points 
    [~,sind] = sort(D,2);
    % positions to be improved  -> find the best matching group 
    optimal_group = sind(:,1); 

    nonoptimal_ratio_0 = 1; 

    for iter = 1:10
        ind_switch  = (groups ~=  optimal_group);
        nonoptimal_ratio = sum(ind_switch) / numel(ind_switch);
        if nonoptimal_ratio > 0
            verbose(0, 'Indexes to be switched:  %3.2g%% positions', nonoptimal_ratio * 100)
        end
        
        if nonoptimal_ratio >= nonoptimal_ratio_0
            break
        end
        nonoptimal_ratio_0 = nonoptimal_ratio; 

        max_dist_0 = inf;
        for i = 1:sum(ind_switch)
            % calculate distance for each point from its group center 
            center_dist = (D(sub2ind(size(D), (1:Npos)', groups)));
            max_dist_0 = max(center_dist(ind_switch));
            % start from the worst case
            ind_worse = find(max(center_dist(ind_switch)) == center_dist, 1, 'first');
            % initial group 
            group_old = groups(ind_worse);
            % better fitting group 
            group_new = optimal_group(ind_worse);
            % position to be switched in the new group 
            ind_new = find(D(:,group_old) == min(D(groups == group_new, group_old)), 1, 'first');
            % switch the group members  
            groups(ind_worse) = group_new;
            groups(ind_new) = group_old;
            ind_switch([ind_worse, ind_new]) = 0;
            if all(ind_switch == 0)
                break
            end
        end
      
%          ind_switch = (groups ~=  sind(:,1));
%          for ii = Ngroups
%              clf
%              hold all;
%              ind = groups == ii;
%              ax = plot(self.probe_positions_0(ind & ind_switch, 1), self.probe_positions_0(ind & ind_switch, 2), 'o');
%              ax2 = plot(self.probe_positions_0(ind & ~ind_switch, 1), self.probe_positions_0(ind & ~ind_switch, 2), 'x');
%              try; ax2.Color = ax.Color; end
%              plot(C(ii,1),C(ii,2),'x','Linewidth', 2)
% %              drawnow 
% %             pause(1)
%          end
%          title(num2str(iter))
%          axis tight equal
%          pause(1)
%        
    end
    
    cluster_refinement_time = cluster_refinement_time + toc(t0); 


    %% optimally sort the indices to help GPU
    [nbins,bins] = hist(groups, unique(groups));
    [~,ind] = sort(nbins,2,'descend');
    for ii = 1:length(bins)
        indices{kk}{ii} = indices_0((groups == bins(ind(ii))));
        scan_ids{kk}{ii} = scans_0((groups == bins(ind(ii))));
    end 
    verbose(2,'=== Number of cluster sizes %i', length(unique(nbins)))
end

verbose(0,'=== Position clusters found in %i iterations in %3.2gs', iter, cluster_time)
verbose(0,'=== Position clusters refined in %i iterations in %3.2gs', iter, cluster_refinement_time)


%rng shuffle 

if verbose()  > 1 && Ngroups_sizes > 1
    warning('Unequal group sizes, it may cause slower calculation')
end


indices_out = horzcat(indices{:});
scan_ids_out =  horzcat(scan_ids{:});


if Ngroups  == 1 && Nsets == 1
    %% merge groups from multiple scans into larger chunks if grouping is too large 
    indices_out = {horzcat(indices_out{:})};
    scan_ids_out =  {horzcat(scan_ids_out{:})};
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


function [groups, C, sum_D, D] = get_best_kmeans(positions, Ngroups)
    % make several guesses to get better Kmean distribution 
    warning('off','stats:kmeans:FailedToConverge')
    for i = 1:10
        [groups{i}, C{i}, sum_D{i}, D{i}] = kmeans(positions, Ngroups);
        nbins = hist(groups{i}, unique(groups{i}));
        score(i) = std(nbins);
    end
    best = math.argmin(score); 
    groups = groups{best};
    C = C{best};
    sum_D = sum_D{best};
    D = D{best};
end

