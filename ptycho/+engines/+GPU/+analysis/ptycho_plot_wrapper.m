% PTYCHO_PLOT_WRAPPER wrapper around the default ptychoshelves plotting routine 
%
% ptycho_plot_wrapper(self, par, fourier_error)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 

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


   
%% PLOTTING 
function ptycho_plot_wrapper(self, par, fourier_error)
    %% wrapper to the cSAXS default plotting function 
    import engines.GPU.GPU_wrapper.*

    p = par.p;
    p.object_size = ceil(p.object_size .*  ( self.Np_p ./ p.asize));  % modify the size in case of presolver with different probe size 
    p.asize = self.Np_p;
    
    p.numobjs = size(self.object,1);
    Nlayers   = size(self.object,2);
    p.object = {};
    for ii = 1:p.numobjs
        p.object_size(ii,:) = self.Np_o;
        p.object{ii} = []; 
        for jj = 1:Nlayers
            p.object{ii}(:,:,1,jj) = Ggather(utils.crop_pad(self.object{ii,jj}, p.object_size));
        end
    end
    p.object_modes = par.object_modes;
    p.probe_modes = par.probe_modes;
    p.probes = [];
    for ii = 1:p.probe_modes
        p.probes(:,:,:,ii) = Ggather(self.probe{ii}(:,:,:,1));
    end
    p.dx_spec=[self.pixel_size]/self.relative_pixel_scale;
    p.engines = {struct()};

    iterations = Ggather(find(any(~isnan(fourier_error),2)));
    p.engines{1}.error_metric_final = struct();
    p.engines{1}.error_metric_final.iteration=iterations;
    p.engines{1}.error_metric_final.value = Ggather(fourier_error( iterations,:));
    p.engines{1}.error_metric_final.method = par.method;
    p.engines{1}.error_metric_final.err_metric = par.likelihood;

    position_offset = 1+floor((p.object_size-self.Np_p)/2);
    for ii = 1:p.numscans
        ind = p.scanidxs{ii}; 
        p.positions(ind,:) = self.modes{1}.probe_positions(ind,[2,1]) + position_offset(p.share_object_ID(ii),:); 
    end

    p.plot.extratitlestring = '';
    p.plot.show_only_FOV = true;
    p.plot.mask_bool = false; 
    p.plot.log_scale = [1 1];
    p.plot.subplwinobj_dir = 'vertical';
    p.plot.show_layers = true;
    p.plot.residua = true; 
        
    if isempty(p.plot.obtitlestring)
        p.plot.obtitlestring =  [core.generate_scan_name(p) ' '];
    end
    if isempty(p.plot.prtitlestring)
        p.plot.prtitlestring = [core.generate_scan_name(p) ' '];
    end
    
    if par.share_object
        p.share_object_ID = ones(p.numobjs,1);
    else
        p.share_object_ID = 1:p.numobjs;
    end
    core.analysis.plot_results(p, 'final', true)
     
end
