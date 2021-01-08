% GET_SCANNING_INDICES simple based method to select indices for DM 
% 
% [indices_out, scan_ids_out] = get_scanning_indices(self, cache, par ) 
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
% 
%       

function [indices_out, scan_ids_out] = get_scanning_indices(self, cache, par )      

import engines.GPU.GPU_wrapper.*
import engines.GPU.shared.*

grouping = par.grouping;
max_groups = 0;

if self.Npos == grouping && par.Nscans == 1
    indices_out = self.reconstruct_ind;
    scan_ids_out = {ones(self.Npos,1)};
    return
end


for kk = 1:par.Nscans
     N = length(self.reconstruct_ind{kk});

     % !! indices ordering has to be always the same for DM !!! 
      
     indices_0 = self.reconstruct_ind{kk};

     max_groups = max(max_groups, ceil(N/grouping));

     for ii = 1:ceil(N/grouping)
         indices{kk}{ii} = indices_0(1+(ii-1)*grouping : min(end,ii*grouping));
     end     
     % fill the last group by the skip indices but do not expand it 
     skip_ind = cache.skip_ind(randperm(length(cache.skip_ind)));
     indices{kk}{end} = [indices{kk}{end}, skip_ind(1:min(end, grouping-length(indices{kk}{end})))]; % join skip_ind back to the last (smallest) set    
     for ii = 1:length(indices{kk})
        scan_ids{kk}{ii} = kk * ones(1,length(indices{kk}{ii}));  % note their scan origin
     end 
end

% how many scans should be merged to reach the desired grouping 
Njoin = ceil(par.grouping / (self.Npos/par.Nscans));


if Njoin > 1 && par.Nscans > 1 && is_method(par, {'PIE', 'ML'})
    % join several scan to improve performance 
    indices_out = cell(ceil(par.Nscans/Njoin),1);
    scan_ids_out = cell(ceil(par.Nscans/Njoin),1);
    %% merge groups from difference scans into larger chunks 
    for kk = 1:ceil(par.Nscans/Njoin)
        indices_out{kk} = [];
        scan_ids_out{kk} = [];
        for ii = 1:Njoin
            if kk+(ii-1)*ceil(par.Nscans/Njoin) <= par.Nscans
            indices_out{kk} = [indices_out{kk}, indices{kk+(ii-1)*ceil(par.Nscans/Njoin)}{1}];
            scan_ids_out{kk} = [scan_ids_out{kk}, scan_ids{kk+(ii-1)*ceil(par.Nscans/Njoin)}{1}];
            end
        end
    end
else
    indices_out = horzcat(indices{:});
    scan_ids_out =  horzcat(scan_ids{:});
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