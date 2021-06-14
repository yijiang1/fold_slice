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
%Note: self.z_distance is first initialized in load_from_p.m, where a
%vacuum layer is appended: self.z_distance = [p.delta_z, inf] for far-field 
par.Nlayers = length(self.z_distance);
if par.Nlayers > 1 && isinf(self.z_distance(end)) 
    % Added by ZC: exclude the last vacuum (inf) layer for multisluce
    par.Nlayers = par.Nlayers - 1;
end

assert(sum(~isfinite(self.z_distance)) <= 1, 'Provided distanced of layers are not possible to be used')

if par.Nlayers > 1 && ~is_method(par, {'PIE', 'ML'})
    error('Multilayer extension is supported only for PIE/ML methods')
end

% Added by ZC. allow user to specify the layer used for position correction
if ~isfield(par,'layer4pos') || isempty(par.layer4pos) 
    par.layer4pos = ceil(par.Nlayers/2);
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
