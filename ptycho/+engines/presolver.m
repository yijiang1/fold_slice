%C_SOLVER external C++ code for DM in lower resolution 

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
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
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

function [ p , fdb] = presolver( p )
    import utils.interpolateFT
    
    scale = p.asize ./ p.asize_presolve; 
    assert(scale(1) == scale(2), 'Nonequal probe size is not supported'); 
    assert(~p.prepare.legacy, 'Presolver with prepare.legacy option is not implemented')
    
    p_small = p; 
    
    p_small.   asize = p.asize_presolve;   % size of the small reconstruction 
    p_small.   prepare_data_filename = []; 
    p_small.   prefix = 'presolver_';
    p_small.   run_name = [p_small.prefix core.generate_scan_name(p) '_' num2str(p_small.asize(1)) 'x' num2str(p_small.asize(2)) '_' p.suffix];  % If empty: automatically generated
    p_small.   plot.prepared_data = false; % plot prepared data
    p_small.   save.store_images = false;               % Write nice jpegs in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
    p_small.   plot.positions = false; 			% Plot the scanning positions
    p_small.   plot.mask_bool = false;          % Mask the noisy contour of the reconstructed object in plots
    p_small = rmfield(p_small, 'ds');           % avoid complains of the initialization code 
    p_small = rmfield(p_small, 'object');       % remove before new initalization 
    p_small = rmfield(p_small, 'detectors');       % remove before new initalization 

    
%     p_small.fmag = utils.crop_pad(p_small.fmag, p.asize_presolve);
%     p_small.fmask = utils.crop_pad(p_small.fmask, p.asize_presolve);
%     p_small.prepare.auto_prepare_data = false; 
      p_small.fmag = [];
      p_small.fmask = [];
      
    if isfield(p_small, 'prepare.legacy_filename')
        p_small.prepare.legacy_filename = replace(p_small.prepare.legacy_filename, sprintf('%ix%i', p.asize), sprintf('%ix%i', p_small.asize));
    end
    % load data nad other settings for small probe 
    [p_small, status] = core.initialize_ptycho(p_small);
   
    % rescale the input estimates to presolver resolution 
    p_small.probes = interpolateFT(p_small.probes,p_small.asize) ;   
    for ii = 1:length(p_small.object)
        p_small.object{ii} = interpolateFT(p.object{1}, ceil(p_small.object_size(ii,:) ) );
    end

    % avoid ML code  !!!
    p_small.opt_iter = 0; 
    
    % run low resolution presolver 
    [p_small, fdb] = engines.c_solver(p_small);
    
    % rescale the input estimates to final resolution 
    p.probes = interpolateFT(p_small.probes,p.asize) ;     
    for ii = 1:length(p_small.object)
        p.object{ii} = interpolateFT(p_small.object{ii},p.object_size(ii,:));
    end
    
    %% renormalize the probe and objects 
    pnorm = math.norm2(p.probes); 
    if p.share_object
        pnorm = mean(pnorm);
    end
    p.probes = p.probes / pnorm; 
    for ii = 1:length(p.object)
        p.object{ii} = p.object{ii} * pnorm(p.share_object_ID(ii));
    end
    
    
    
    if isfield(p_small, 'error_metric')
        p.error_metric = p_small.error_metric;
    else
        p.error_metric = []; 
    end


   p = rmfield(p,  'asize_presolve');


end