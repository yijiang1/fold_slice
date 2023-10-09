function [outputDir, paramInfo] = generateResultDir(param, resultDir, extra)
%Generate output directory for GPU and GPU_MS engines
%   Written by YJ

paramInfo = strcat(param.method,'_',param.opt_errmetric,'_p',num2str(param.probe_modes),'_g',num2str(param.grouping));

if isfield(param,'asize_presolve') && length(param.asize_presolve)==2
    paramInfo = strcat(paramInfo,'_Ndp',num2str(param.asize_presolve(1)));
end

if strcmp(param.method, 'MLc') && param.accelerated_gradients_start < param.number_iterations    
    paramInfo = strcat(paramInfo,'_ag',num2str(param.accelerated_gradients_start));
end
if strcmp(param.method, 'MLc') && param.momentum ~=0  
    paramInfo = strcat(paramInfo,'_m',num2str(param.momentum));
end

if param.beta_object < 1
    paramInfo = strcat(paramInfo,'_betaO',num2str(param.beta_object));
end
if param.beta_probe < 1
    paramInfo = strcat(paramInfo,'_betaP',num2str(param.beta_probe));
end

%{
if isfield(param,'beta_LSQ')
    paramInfo = strcat(paramInfo,'_betaLSQ',num2str(param.beta_LSQ));
end

if param.delta_p ~= 0.1
    paramInfo = strcat(paramInfo,'_LSQdamping',num2str(param.delta_p));
end
%}

if param.probe_position_search < param.number_iterations
    paramInfo = strcat(paramInfo,'_pc',num2str(param.probe_position_search));
    if isfield(param,'apply_relaxed_position_constraint') && ~param.apply_relaxed_position_constraint
        paramInfo = strcat(paramInfo,'_noModel');    
    elseif ~isempty(param.probe_geometry_model)
        if ismember('scale', param.probe_geometry_model)
            paramInfo = strcat(paramInfo,'_scale');
        end
        if ismember('asymmetry', param.probe_geometry_model)
            paramInfo = strcat(paramInfo,'_asym');
        end
        if ismember('rotation', param.probe_geometry_model)
            paramInfo = strcat(paramInfo,'_rot');
        end
        if ismember('shear', param.probe_geometry_model)
            paramInfo = strcat(paramInfo,'_shear');
        end
        %for i=1:length(param.probe_geometry_model)
        %    paramInfo = strcat(paramInfo,'_',param.probe_geometry_model{i});
        %end
    end
    if isfield(param,'update_pos_weight_every') && param.update_pos_weight_every< param.number_iterations
        paramInfo = strcat(paramInfo,'_updW',num2str(param.update_pos_weight_every));
    end
    if param.probe_position_error_max < inf
        %paramInfo = strcat(paramInfo,'_maxError',num2str(param.probe_position_error_max/1e-9),'nm');
        paramInfo = strcat(paramInfo,'_maxError',num2str(param.probe_position_error_max));
    end
    if isfield(param,'max_pos_update_shift') && param.max_pos_update_shift~=0.1
        paramInfo = strcat(paramInfo,'_maxUpdShift',num2str(param.max_pos_update_shift));
    end
    if isfield(param,'probe_position_search_momentum') && param.probe_position_search_momentum>0
        paramInfo = strcat(paramInfo,'_m',num2str(param.probe_position_search_momentum));
    end
end

if param.detector_scale_search < param.number_iterations
    paramInfo = strcat(paramInfo,'_detScaleSearch',num2str(param.detector_scale_search));
end

if param.probe_fourier_shift_search < param.number_iterations
    paramInfo = strcat(paramInfo,'_fpc',num2str(param.probe_fourier_shift_search));
end

if param.background>0
    paramInfo = strcat(paramInfo,'_bg',num2str(param.background));
end
if param.delta>0
    paramInfo = strcat(paramInfo,'_delta',num2str(param.delta));
end
if param.reg_mu>0
    paramInfo = strcat(paramInfo,'_regSmooth',num2str(param.reg_mu));
end 
if isfield(param,'TV_lambda') && param.TV_lambda>0
    paramInfo = strcat(paramInfo,'_TV',num2str(param.TV_lambda));
end 

if param.positivity_constraint_object>0
    paramInfo = strcat(paramInfo,'_posObj',num2str(param.positivity_constraint_object));
end

if param.variable_probe
    paramInfo = strcat(paramInfo,'_vp',num2str(param.variable_probe_modes));
    if param.variable_probe_smooth>0
        paramInfo = strcat(paramInfo,'_smooth',num2str(param.variable_probe_smooth));
    end
end
if param.variable_intensity
    paramInfo = strcat(paramInfo,'_vi');
end
if param.apply_multimodal_update
    paramInfo = strcat(paramInfo,'_mm');
end

if is_used(param, 'fly_scan')
    paramInfo = strcat(paramInfo,'_apFly',num2str(param.Nmodes));
end

if ~isempty(param.delta_z)
    paramInfo = strcat(paramInfo,'_Ns',num2str(length(param.delta_z)));
    paramInfo = strcat(paramInfo,'_dz',num2str(mean(param.delta_z)));
    if param.regularize_layers~=0
        paramInfo = strcat(paramInfo,'_reg',num2str(param.regularize_layers));
    end
    if param.preshift_ML_probe
        paramInfo = strcat(paramInfo,'_centerProbe');
    end
end

if isfield(param,'diff_pattern_blur') && param.diff_pattern_blur>0
    paramInfo = strcat(paramInfo,'_dpBlur',num2str(param.diff_pattern_blur));
end

if any(param.custom_data_flip)
    paramInfo = strcat(paramInfo,'_dpFlip');
    if param.custom_data_flip(1)==1
        paramInfo = strcat(paramInfo,'_lr');
    end
    if param.custom_data_flip(2)==1
        paramInfo = strcat(paramInfo,'_ud');
    end
    if param.custom_data_flip(3)==1
        paramInfo = strcat(paramInfo,'_T');
    end
end
%{
if any(param.det_bad_pixels(:))
    output = strcat(output,'_badPixels');
end
%}
if nargin==3
    paramInfo = strcat(paramInfo,extra);
end

outputDir = strcat(resultDir,'/',paramInfo,'/');
end
