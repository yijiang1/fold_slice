clear variables
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)
utils.ccc
%%
alpha_max = optimizableVariable('alpha_max',[1, 30]); %mrad    
df = optimizableVariable('df',[-1000, 1000]); %nm
rot_ang = optimizableVariable('rot_ang',[-180, 180]); %nm
cs = optimizableVariable('cs',[0, 1]); %mm
dp_transpose = optimizableVariable('dp_transpose',[0, 1], 'Type', 'integer'); 
scan_step_size = optimizableVariable('scan_step_size',[0.4, 10]); %angstrom            
N_dp_factor = optimizableVariable('N_dp_factor',[1, 4], 'Type', 'integer'); 

verbose = 2;
%plot_funcs = {};
plot_funcs = {@plotObjectiveModel, @plotMinObjective};

fun = @(x)ptycho_recon_bilayer_MoSe2(x.df, 30, 0);
        results = bayesopt(fun, [df],...
            'Verbose',verbose,...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'IsObjectiveDeterministic',false,...
            'MaxObjectiveEvaluations', 25,...
            'NumSeedPoints',3,...
            'PlotFcn',plot_funcs);
        
%% line search
df_s = -1000:5:1000;
data_errors = zeros(1,length(df_s));
for i=1:length(df_s)
    disp(df_s(i))
    data_errors(i) = ptycho_recon_bilayer_MoSe2(df_s(i), 30, 0);
end
%%
close all
plot(df_s,data_errors, '.', 'MarkerSize', 9)
ylabel('Data error [a.u.]')
xlabel('defocus [angstrom]')
