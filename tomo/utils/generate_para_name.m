function [paraName] = generate_para_name(par, extra)
%Generate output names based on alignment parameters
%   Written by YJ.

paraName = strcat('_hps',num2str(par.high_pass_filter));
if isfield(par,'angle_correction') && par.angle_correction~=0
    paraName = strcat(paraName,'_angle_corr',num2str(par.angle_correction));
end
if par.binning~=1
    paraName  = strcat(paraName,'_bin',num2str(par.binning));
end
if par.use_localTV
    paraName = strcat(paraName,'_TV',num2str(par.localTV_lambda));
end
if par.apply_positivity
    paraName = strcat(paraName,'_pos');
end
if par.position_update_smoothing>0
    paraName = strcat(paraName,'_smooth',num2str(par.position_update_smoothing));
end
paraName = strcat(paraName,'_sc',num2str(par.min_step_size));

if nargin==2
    paraName = strcat(paraName,extra);
end

end

