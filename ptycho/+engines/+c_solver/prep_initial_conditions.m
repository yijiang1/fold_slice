%PREP_INITIAL_CONDITIONS prepares the initial_conditions file for the
% reconstruction with external C++ code
% save2hdf5 is used to save the required datasets and attributes to an h5 file
%
% ** p      p structure
%
% see also: engines.c_solver
%

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

function prep_initial_conditions( p )
import utils.verbose
import io.HDF.save2hdf5

% define structure and append attributes
h5_struc = [];

h5_struc.Attributes.pfft_relaxation = p.pfft_relaxation;
h5_struc.Attributes.probe_regularization = p.probe_regularization;
h5_struc.Attributes.probe_radius = p.probe_support_radius;
h5_struc.Attributes.diffmap_iterations = uint32(p.number_iterations);
h5_struc.Attributes.max_mlh_iterations = uint32(p.opt_iter);

% for multi-slice c-code:
% currently only allow the same initial guess for all slices

%if isfield(p,'N_layer') && strcmpi(p.engines{1}.name, 'c_solver') && (isfield(p,'object') && var(p.object{1}(:))>1e-5) %&& strcmpi(p.initial_iterate_object,'file')
if isfield(p,'N_layer') && strcmpi(p.engines{1}.name, 'c_solver') && ~isempty(p.initial_iterate_object_file{1})

    for ii=p.share_object_ID
%% Simple initial guess 
        if size(p.object{ii},4) < p.N_layer
            extra_layers = p.N_layer-size(p.object{ii},4); 
            p.object{ii} = cat(4, p.object{ii}, ones([p.object_size(ii,:),p.object_modes,extra_layers]));
            % get the empty layers on top / bottom
            p.object{ii} = circshift(p.object{ii}, floor(extra_layers/2),4);
        end
%% Esther version, unrealiable .... 
%         if ~isfield(p,'ms_init_ob_fraction') || sum(p.ms_init_ob_fraction)~=1 || length(p.ms_init_ob_fraction)~=p.N_layer
%             verbose(3, '## p.ms_init_ob_fraction: using default 1/N_layer');
%             p.ms_init_ob_fraction = ones(1,p.N_layer) / p.N_layer;
%         end
%         if sum(p.ms_init_ob_fraction==1)==0  
%             for ii=p.share_object_ID
%                 ob_phase = engines.ML_MS.fun_ramp_unwrap(p.object{ii}, p.asize);
%                 ob_abs = abs(p.object{ii});
%                 for n = 1:p.N_layer
%                     p.object{ii}(:,:,1,n) = ob_abs.^(p.ms_init_ob_fraction(1)) .* exp(1i*ob_phase.*p.ms_init_ob_fraction(n));
%                 end
%             end	      
%         elseif size(p.object{1},4)==1 % No p.object_layers loaded; No unwrapping
%             for ii=p.share_object_ID
%                 ob = p.object{ii};
%                 for n = 1:p.N_layer
%                     if p.ms_init_ob_fraction(n)==1
%                         p.object{ii}(:,:,1,n) = ob;
%                     else
%                         p.object{ii}(:,:,1,n) = ones(size(ob))*(1+1i*1e-8);
%                     end
%                 end
%             end	
%         end        
%         if 1
%             figure(100);
%             for n = 1:p.N_layer
%                 subplot(1,p.N_layer,n); 
%                 imagesc(angle(p.object{1}(:,:,1,n))); axis xy equal tight; colorbar
%                 title(sprintf('Initial condition (phase), layer %d',n));
%             end
%             drawnow
%         end

    end
end

if ~isempty(p.delta_z) && p.preshift_ML_probe
   % if multilayer extension is used, assume that the probe is
   % reconstructed at the center plane of the sample -> apply shift 
   probe_offset = -sum(p.delta_z)/2; 
   p.probes = utils.prop_free_nf(p.probes, p.lambda , probe_offset, p.dx_spec(1)) ;
end



% define object and probe datasets
% axes have to be permuted
for ii=p.share_object_ID
    h5_struc.objects.(['object_' num2str(ii-1)]) = permute(p.object{ii}, [2 1 3 4]);
end
for ii=p.share_probe_ID
    h5_struc.probes.(['probe_' num2str(ii-1)]) = permute(squeeze(p.probes(:,:,ii,:)), [2 1 3]);
end

% write data to disk, using a compression level as defined in
% template_ptycho
verbose(3,'Writing H5 initial conditions file: %s',[p.initial_conditions_path p.initial_conditions_file]);
save2hdf5([p.initial_conditions_path p.initial_conditions_file], h5_struc, 'overwrite', true, 'comp', p.io.file_compression);


% update object size defined in the prepared data
verbose(3,'Updating measurement file: %s',[p.prepare_data_path p.prepare_data_filename]);
for ii = 1:p.numobjs
    objects(:,ii) = uint64([size(p.object{ii},1),size(p.object{ii},2)]);
end 
save2hdf5([p.prepare_data_path p.prepare_data_filename], objects, 'data_name', 'objects', 'overwrite', false);

% save updated positions into the prepared data 
% save2hdf5([p.prepare_data_path p.prepare_data_filename], uint32(p.positions'), 'data_name', 'measurement/n0/positions', 'overwrite', false);




end

