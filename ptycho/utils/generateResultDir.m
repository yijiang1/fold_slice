function [outputDir, paramInfo] = generateResultDir(param, resultDir, extra)
%Generate output directory for GPU engine
%   Written by YJ.

paramInfo = strcat(param.method,'_',param.opt_errmetric,'_p',num2str(param.probe_modes),'_g',num2str(param.grouping));
if isfield(param,'asize_presolve') && length(param.asize_presolve)==2
    paramInfo = strcat(paramInfo,'_Ndp',num2str(param.asize_presolve(1)));
end

if strcmp(param.method, 'MLc') && param.accelerated_gradients_start < param.number_iterations    
    paramInfo = strcat(paramInfo,'_accGrad',num2str(param.accelerated_gradients_start));
end

if strcmp(param.method, 'MLc') && param.momentum ~=0  
    paramInfo = strcat(paramInfo,'_mom',num2str(param.momentum));
end

if param.beta_object < 1
    paramInfo = strcat(paramInfo,'_betaObj',num2str(param.beta_object));
end
if param.beta_probe < 1
    paramInfo = strcat(paramInfo,'_betaProb',num2str(param.beta_probe));
end

if isfield(param,'beta_LSQ')
    paramInfo = strcat(paramInfo,'_betaLSQ',num2str(param.beta_LSQ));
end
if param.delta_p ~= 0.1
    paramInfo = strcat(paramInfo,'_LSQdamping',num2str(param.delta_p));
end

if param.probe_position_search < param.number_iterations
    paramInfo = strcat(paramInfo,'_pc',num2str(param.probe_position_search));
    if ~isempty(param.probe_geometry_model)
        paramInfo = strcat(paramInfo,'_model');
        for i=1:length(param.probe_geometry_model)
            paramInfo = strcat(paramInfo,'_',param.probe_geometry_model{i});
        end
    end
    if param.probe_position_error_max < inf
        paramInfo = strcat(paramInfo,'_maxPosError',num2str(param.probe_position_error_max/1e-9),'nm');
    end
    if isfield(param,'apply_relaxed_position_constraint') && ~param.apply_relaxed_position_constraint
        paramInfo = strcat(paramInfo,'_noModelCon');
    end
    if isfield(param,'update_pos_weight_every') && param.update_pos_weight_every< param.number_iterations
        paramInfo = strcat(paramInfo,'_updPosWeight',num2str(param.update_pos_weight_every));
    end
    if isfield(param,'max_pos_update_shift') && param.max_pos_update_shift~=0.1
        paramInfo = strcat(paramInfo,'_maxPosUpdShift',num2str(param.max_pos_update_shift));
    end
    if isfield(param,'probe_position_search_momentum') && param.probe_position_search_momentum>0
        paramInfo = strcat(paramInfo,'_mom',num2str(param.probe_position_search_momentum));
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
    paramInfo = strcat(paramInfo,'_posConObj',num2str(param.positivity_constraint_object));
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
%if isfield(param, 'extension') && param.extension{} = {'fly_scan'}
if is_used(param, 'fly_scan')
    paramInfo = strcat(paramInfo,'_afs',num2str(param.Nmodes));
end

if ~isempty(param.delta_z)
    %paramInfo = strcat(paramInfo,'_Nlayers',num2str(length(param.delta_z)+1));
    paramInfo = strcat(paramInfo,'_Nlayers',num2str(length(param.delta_z)));
    paramInfo = strcat(paramInfo,'_avg_dz',num2str(mean(param.delta_z)));
    if param.regularize_layers~=0
        paramInfo = strcat(paramInfo,'_reg',num2str(param.regularize_layers));
    end
    if param.preshift_ML_probe
        paramInfo = strcat(paramInfo,'_centerProbe');
    end
end

if param.custom_data_flip(1)==1
    paramInfo = strcat(paramInfo,'_fliplr');
end
if param.custom_data_flip(2)==1
    paramInfo = strcat(paramInfo,'_flipud');
end
if param.custom_data_flip(3)==1
    paramInfo = strcat(paramInfo,'_transpose');
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

