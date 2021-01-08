% INITIAL_CHECKS check if the inputs are valid or try to correct them 
% 
% [self,par] = initial_checks(self, par)
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
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

function [self,par] = check_inputs(self, par)
    import engines.GPU_MS.shared.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*

[self.Np_o(1),self.Np_o(2),~] = size(self.object{1});
[self.Np_p(1),self.Np_p(2),~] = size(self.probe{1});
par.Nrec  = 1;
par.Nscans = length(self.reconstruct_ind);


if ischar(par.extension)
    par.extension = {par.extension};
end

for ii = 1:numel(self.object)
   assert(all(isfinite(self.object{ii}(:))), 'Provided object contains nan / inf')
end

for ii = 1:numel(self.probe)
   assert(all(isfinite(self.probe{ii}(:))), 'Provided probes contain nan / inf')
end

Np_d =  size(self.diffraction);
if any(self.Np_p ~= Np_d(1:2)) % && isempty(self.modF_ROI)
    error('Size of probe and data is different')
end


tmp = self.diffraction(1:self.Np_p(1)*7:end);  % get some small sample 
tmp = tmp * 2^(par.upsampling_data_factor*2); % remove upsampling effects 
if par.compress_data && max(abs((tmp - round(tmp)))) > 0.2
    verbose(1,'Data are not integers, cannot use compression')
    par.compress_data = false;
end

%%%%%%%%%%%%% accelerated solver %%%%%%%%%%%%%%%%%%%

if  par.accelerated_gradients_start < par.number_iterations && ~is_method(par, 'MLs')
    verbose(3, 'accelerated_gradients_start < number_iterations	 is supported only for MLc engine')
    par. accelerated_gradients_start = inf; 
end


% if  par.accelerated_gradients_start < par.number_iterations && par.momentum > 0 && is_method(par, 'ML')
%     error('accelerated_gradients_start < inf cannot be used if momemtum > 0 ')
% end


%%%%%%%%%%%%% variable probe %%%%%%%%%%%%%%%%%%%%%%%%

if ~par.variable_probe
    par.variable_probe_modes = 0; 
end

if par.variable_probe && par.variable_probe_modes  > 0 && ~is_method(par, {'PIE', 'ML'})
    warning('Variable probe implemented only for PIE and ML')
    par.variable_probe = false;
end


if par.variable_probe &&  par.variable_probe_modes == 0 
    error('Choose more than 0 variable_probe_modes for OPRP')
    par.variable_probe_modes = 1; 
end

if par.variable_probe && ~par.share_probe && is_method(par, 'PIE')
    par.share_probe = true;
    % variable probe means automatically shared variable probe 
end


if ~is_method(par, {'PIE', 'ML'}) && strcmpi(par.likelihood, 'poisson')
    warning('Poisson likelihood supported only for PIE methods')
    par.likelihood = 'L1';
end
if ~ismember(lower( par.likelihood), {'l1','poisson'})
    error('Unsupported error estimation')
end

   
%%%%%%%%%%%%%% check if position correction is allowed 
if ~ is_method(par, {'PIE', 'ML'}) && par.probe_position_search < par.number_iterations
    verbose(2, 'Position correction supported only for PIE/ML methods ')
    par.probe_position_search = inf; 
end

if any(~ismember(par.probe_geometry_model, {'scale', 'asymmetry', 'rotation', 'shear'}))
   missing_option = setdiff(par.probe_geometry_model, {'scale', 'asymmetry', 'rotation', 'shear'} ); 
   error('Unsupported geometry model option:  "%s"', missing_option{1}) 
end

if par.probe_position_search < par.number_iterations && par.detector_scale_search < par.number_iterations && any(ismember(par.probe_geometry_model,'scale')) 
    error('Do not use probe_position_search with probe_geometry_model==''scale'' and detector_scale_search together')
end

%%%%%% checks for the multilayer method  %%%%%%%%%%%%%%%%
par.Nlayers = length(self.z_distance); %Note: z_distance is initialied in load_from_p.m
% Added by ZC: exclude last inf layer for multisluce
if par.Nlayers > 1 && isinf(self.z_distance(end)) 
    if isfield(par,'rmvac') && par.rmvac
        par.Nlayers = par.Nlayers - 1;
    end
end

assert(sum(~isfinite(self.z_distance)) <= 1, 'Provided distanced of layers are not possible to be used')

if par.Nlayers > 1 && ~is_method(par, {'PIE', 'ML'})
    error('Multilayer extension is supported only for PIE/ML methods')
end
% Modified by ZC: multislice now works with multiple probe modes
% if par.Nlayers > 1 && par.probe_modes > 1
%     error('Multilayer extension is not supported with incoherent modes')
% end
if par.Nlayers > 1 && par.Nscans > 1
	error('Multilayer extension is not supported with multiple scans')
end

%%%%%%%%%%  fast scanning %%%%%%%%%%%%%%%%%%%%%%%%%
if is_used(par, 'fly_scan') && ~is_method(par, {'PIE', 'ML'})
    error('Fly scan is  supported only for PIE/ML methods')
end
if is_used(par, 'fly_scan')
    if par.Nmodes == 1
        warning('Flyscan has no effect with a single mode')
        par.extension = setdiff(par.extension, 'fly_scan');
        par.apply_subpix_shift= true;
    end
    par.Nrec = par.Nmodes;
%     par.apply_multimodal_update = true;
end

%%%%%%%% nearfield %%%%%%%%%%%%%%%%%%%%%%%%%%%

if par.estimate_NF_distance < par.number_iterations && isinf(self.z_distance(end))
    error('estimate_NF_distance valid only for nearfield mode')
end


%%%%%%%%%%%%%%% OTHER %%%%%%%%%%%%%%%%%%%%
if strcmpi(par.likelihood, 'poisson') && par.background_detection && ~isinf(par.background_detection)
    error('Background detection does not work well with Poisson likelihood')
end

if prod(self.Np_p) *self.Npos >  intmax('int32') && par.keep_on_gpu && is_method(par, {'MLs', 'ePIE'})
    warning('Dataset as more than 2147483647 elements (max of int32), in case of problems try par.keep_on_gpu = false')
end



if any(self.noise(:) == 0) && par.relax_noise
    warning('Some values of expected noise are 0')
    self.noise = max(0.5, self.noise); 
end


if par.Nrec > max([par.Nmodes, par.probe_modes , par.object_modes])
    warning('Number of modes is too high')
end

if length(self.probe_positions) ~= self.Npos
   self.probe_positions = [];
end


if par.mirror_objects && par.Nscans ~= 2 
    error('Object mirroring is supported only for two scans')
end

%%%%%% position correction %%%%% 
if ~is_method(par, {'PIE', 'ML'}) && par.probe_position_search < par.number_iterations
    warning('Position corrections works only for PIE/ML methods')
end
if is_method(par, {'PIE', 'ML'}) && par.probe_position_search< par.number_iterations && ~(par.apply_subpix_shift || is_used(par,'fly_scan'))
   verbose(2,'Subpixel shifting is strongly recommended for position refinement => enforcing par.apply_subpix_shift = true') 
   par.apply_subpix_shift = true; 
end

    

end
