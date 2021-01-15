%simulate CBED for FSC analysis
%% parameters
addpath(fullfile(pwd,'utils_electron'))
FOV = 60; %fixed FOV in angstrom
N = 128;   % size of diffraction pattern in pixels. only square dp allowed 
scan_step_zie = 4; %scan step size in angstrom
N_scans_x = round(FOV/scan_step_zie); % number of scan positions along horizontal direction
N_scans_y = round(FOV/scan_step_zie); % number of scan positions along vertical direction

maxPosError = 0; %largest randrom position error
dose = 2e3; %total electron dose (e/A^2)
Nc_avg = dose*scan_step_zie^2/N^2; %average electron count per detector pixel. For poisson noise, SNR = sqrt(Nc_avg);
base_dir = '/home/beams2/YJIANG/research/algorithm/simulation/FSC_study/electron_ptycho_temp/';

%% load test object
disp('Loading test object...')
load(fullfile(pwd,'utils_electron','CuPcCl.mat'))
%pad object in case of large FOV
phase_true = padarray(phase_psi,[6400,6400],'circular','post');
r = 4; %resample phase
phase_true = imresize(phase_true, 1/r);
%create a complex object
object_true = ones(size(phase_true)).*exp(1i*phase_true);
dx = dx*r; %real-space pixel size in angstrom

N_obj = size(object_true,1); %only square object allowed
ind_obj_center = floor(N_obj/2)+1;

%% generate probe function
disp('Generating probe function...')
df = 600; %defocus in angstrom
cs = 0;
voltage = 300; %beam voltage in keV
alpha_max = 18; %convergence angle in mrad
show_figure = true;
[probe_true, mask, ~] = generateProbeFunction(dx,N, 0,0,df ,cs, 1, voltage, alpha_max, show_figure );
%calculate rbf
lambda = 12.398/sqrt((2*511.0+voltage).*voltage); %angstrom
dk = 1/dx/N; %fourier-space pixel size in 1/A
rbf = alpha_max/1e3/lambda/dk;
%% save initial probe and parameters
probe = probe_true;
data_dir = strcat('CuPcCl_ss',num2str(scan_step_zie),'_a',num2str(alpha_max),'_df',num2str(df),'_dose',num2str(dose),'/');
mkdir(fullfile(base_dir,data_dir))
p = {};
p.binning = false;
p.detector.binning = false;
p.dk = dk;
p.voltage = voltage;
p.df = df;
p.alpha_max = alpha_max;
p.cs = cs;
p.N_scans_x = N_scans_x;
p.N_scans_y = N_scans_y;
save(strcat(base_dir,data_dir,'init_probe'),'probe','p')
%% generate two datasets (required by FSC)
for j=1:2
    disp(j)
    % generate scan positions
    close all
    disp('Generating scan positions...')

    pos_x = (1 + (0:N_scans_x-1) *scan_step_zie);
    pos_y = (1 + (0:N_scans_y-1) *scan_step_zie);
    % centre this
    pos_x  = pos_x - (mean(pos_x));
    pos_y  = pos_y - (mean(pos_y));
    [Y,X] = meshgrid(pos_x, pos_y);
    ppX = X(:); % true posoition
    ppY = Y(:);

    ppX_recon_init = ppX;
    ppY_recon_init = ppY;

    %add random position errors - to simulate scan noise
    ppX_recon_init = ppX_recon_init + maxPosError*(rand(size(ppX_recon_init))*2-1);
    ppY_recon_init = ppY_recon_init + maxPosError*(rand(size(ppY_recon_init))*2-1);

    %
    disp('Generating diffraction patterns...')
    %calculate indicies for all scans
    N_scan = length(ppX);
    %position = pi(integer) + pf(fraction)
    py_i = round(ppY/dx);
    py_f = ppY - py_i*dx;
    px_i = round(ppX/dx);
    px_f = ppX - px_i*dx;

    ind_x_lb = px_i - floor(N/2) + ind_obj_center;
    ind_x_ub = px_i + ceil(N/2) -1 + ind_obj_center;
    ind_y_lb = py_i - floor(N/2) + ind_obj_center;
    ind_y_ub = py_i + ceil(N/2) -1 + ind_obj_center;

    dp = zeros(N,N,N_scan);
    dp_true = zeros(N,N,N_scan);

    snr = ones(N_scan,1)*inf; %signal-to-noise ratio of each diffraction pattern

    f = waitbar(0,'1','Name','Simulating diffraction patterns...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);

    for i=1:N_scan
        % Check for clicked Cancel button
        if getappdata(f,'canceling')
            break
        end
        % Update waitbar and message
        waitbar(i/N_scan,f,sprintf('No.%d/%d',i,N_scan))

        probe_s = shift(probe_true, dx, dx, px_f(i), py_f(i));
        obj_roi = object_true(ind_y_lb(i):ind_y_ub(i),ind_x_lb(i):ind_x_ub(i));
        psi =  obj_roi .* probe_s;

        %FFT to get diffraction pattern
        dp_true(:,:,i) = abs(fftshift(fft2(ifftshift(psi)))).^2;
        dp(:,:,i) = dp_true(:,:,i);

        %Add poisson noise
        if Nc_avg<inf
            dp_true_temp = dp_true(:,:,i);
            dp_temp = dp_true_temp/sum(dp_true_temp(:))*(N^2*Nc_avg);
            dp_temp = poissrnd(dp_temp);
            dp_temp = dp_temp*sum(dp_true_temp(:))/(N^2*Nc_avg);
            snr(i) = mean((dp_true_temp(:)))/std(dp_true_temp(:) - dp_temp(:));
            dp(:,:,i) = dp_temp;
        end
    end
    delete(f)
    disp('Generating diffraction patterns...done')

    % save cbed 
    disp('Saving diffraction patterns...')
    save_dir = fullfile(base_dir,data_dir,strcat('data',num2str(j)));
    mkdir(save_dir)
    
    save_name = strcat('data_roi0_dp.hdf5'); %save diffraction patterns
    h5create(fullfile(save_dir,save_name), '/dp', size(dp),'ChunkSize',[size(dp,1) size(dp,1), 1],'Deflate',4)
    h5write(fullfile(save_dir,save_name), '/dp', dp*100)
    save_name = strcat('data_roi0_para.hdf5'); %save scan positions
    hdf5write(fullfile(save_dir,save_name), '/ppX', ppX)
    hdf5write(fullfile(save_dir,save_name), '/ppY', ppY,'WriteMode','append')
    disp('done')
end

%% run ptychosheleves script to prepare reconstruction parameters
template = 'ptycho_electron_simulation_template';
% you can adjust more parameters in the template
base_path_ptycho = fullfile(base_dir,data_dir); %base path needed by ptychoshelves
grouping = N_scans_x; %adjust group size based on total # of diffraction patterns
run(template)

%% Run the reconstruction
tic
out = core.ptycho_recons(p);
toc

%% To exam the final FSC score, load the .mat file generated by PtychoSheleves
%for example:
load(fullfile(eng.fout,'Niter200.mat'))
fsc = p.dx_spec(1)/outputs.fsc_score{end}.resolution; %unit: angstrom
disp(fsc)