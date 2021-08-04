% pout = ptycho_model_probe(p)

function pout = ptycho_model_probe(p)
import utils.*
import io.*

% Define often-used variables
lambda = p.lambda;
asize = p.asize; % Diffr. patt. array size  
% added by YJ for up-sampled diffraction patterns.
if p.detector.upsampling >0
	asize = asize*2^p.detector.upsampling;
end
dx_spec = p.dx_spec;
a2 = prod(asize);

if check_option(p, 'prop_regime', 'nearfield')
    % for this task use original values of the pixel sizes
    dx_spec =  p.lambda*p.z*p.nearfield_magnification ./ (p.asize*p.ds);
end

% Prepare probe
if p.model_probe
    % STEM probe: based on Eq.(2.10) in Advanced Computing in Electron 
    % Microscopy (2nd edition) by Dr.Kirkland
    if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
        df = p.model.probe_df;
        alpha_max = p.model.probe_alpha_max;
        amax = alpha_max*1e-3; %in rad
        amin = 0;

        klimitmax = amax/lambda;
        klimitmin = amin/lambda;
        N = asize(1);
        dk = 1/(dx_spec(1)*N);

        kx = linspace(-floor(N/2),ceil(N/2)-1,N);
        [kX,kY] = meshgrid(kx,kx);

        kX = kX.*dk;
        kY = kY.*dk;
        kR = sqrt(kX.^2+kY.^2);
        theta = atan2(kY,kX);
        
        mask = single(kR<=klimitmax).*single(kR>=klimitmin);
        chi = -pi*lambda*kR.^2*df;
        %third-order spherical aberration in angstrom
        if isfield(p.model,'probe_c3') && p.model.probe_c3~=0
            chi = chi + pi/2*p.model.probe_c3*lambda^3*kR.^4;
        end
        %fifth-order spherical aberration in angstrom
        if isfield(p.model,'probe_c5') && p.model.probe_c5~=0
            chi = chi + pi/3*p.model.probe_c5*lambda^5*kR.^6;
        end
        %seventh-order spherical aberration in angstrom
        if isfield(p.model,'probe_c7') && p.model.probe_c7~=0
            chi = chi + pi/4*p.model.probe_c7*lambda^7*kR.^8;
        end
        %twofold astigmatism in angstrom & azimuthal orientation in radian
        if isfield(p.model,'probe_f_a2') && isfield(p.model,'probe_theta_a2') && p.model.probe_f_a2~=0
        	chi = chi + pi*p.model.probe_f_a2*lambda*kR.^2*sin(2*(theta-p.model.probe_theta_a2));
        end
        %threefold astigmatism in angstrom & azimuthal orientation in radian
        if isfield(p.model,'probe_f_a3') && isfield(p.model,'probe_theta_a3') && p.model.probe_f_a3~=0
        	chi = chi + 2*pi/3*p.model.probe_f_a3*lambda^2*kR.^3*sin(3*(theta-p.model.probe_theta_a3));
        end
        %coma in angstrom & azimuthal orientation in radian
        if isfield(p.model,'probe_f_c3') && isfield(p.model,'probe_theta_c3') && p.model.probe_f_c3~=0
        	chi = chi + 2*pi/3*p.model.probe_f_c3*lambda^2*kR.^3*sin(theta-p.model.probe_theta_c3);
        end
        
        probe = mask.*exp(-1i.*chi);
        probe = fftshift(ifft2(ifftshift(probe)));
        probe = probe/sum(sum(abs(probe)));
    else  %X-ray probe  
        if p.model.probe_is_focused
            verbose(2, 'Using focused probe as initial model.');
            if asize(1) ~= asize(2)
                error('Focused probe modeling is only implemented for square arrays (please feel free to change that).');
            end

            if isempty(p.model.probe_zone_plate_diameter) || isempty(p.model.probe_outer_zone_width)
                zp_f = p.model.probe_focal_length;
                verbose(3, 'Using model.probe_focal_length for modeled probe.');
            else
                zp_f = p.model.probe_zone_plate_diameter * p.model.probe_outer_zone_width / lambda;
            end

            % The probe is generated in a larger array to avoid aliasing
            upsample = p.model.probe_upsample;

            defocus = p.model.probe_propagation_dist;
            Nprobe = upsample*asize(1);                       % Array dimension for the simulation
            dx = (zp_f+defocus)*lambda/(Nprobe*dx_spec(1));   % pixel size in the pupil plane
            r1_pix = p.model.probe_diameter / dx;              % size in pixels of first pinhole
            r2_pix = p.model.probe_central_stop_diameter / dx;       % size in pixels of central stop


            % Pupil
            [x,y] = meshgrid(-Nprobe/2:floor((Nprobe-1)/2),-Nprobe/2:floor((Nprobe-1)/2));
            r2 = x.^2 + y.^2;
    %         w = (r2 < (r1_pix)^2);
            if upsample*asize(1) < round(r1_pix)-5
                error(sprintf('For this experimental parameters asize must be at least %d in order for the lens to fit in the window.',ceil((round(r1_pix)-5)./upsample+1)))
            end
            w = fftshift(filt2d_pad(upsample*asize(1), round(r1_pix)+5, round(r1_pix)-5, 'circ'));
            if p.model.probe_central_stop
                w = w .*(1-fftshift(filt2d_pad(upsample*asize(1), round(r2_pix)+2, round(r2_pix-2), 'circ')));
            end
            if isfield(p.model,'probe_structured_illum_power') && p.model.probe_structured_illum_power
                rng default
                r = utils.imgaussfilt2_fft(randn(upsample*p.asize),upsample*2); 
                r = r / math.norm2(r); 
                r = exp(1i*r*p.model.probe_structured_illum_power);
                w = imgaussfilt(w,upsample/2).*r; 
            end

            % Propagation
            probe_hr = prop_free_ff(w .* exp(-1i * pi * r2 * dx^2 / (lambda * zp_f)), lambda, zp_f + defocus, dx);

            % Cropping back to field of view
            probe = crop_pad(probe_hr, asize); 

            % prevent unreal sharp edges from the cropped tails in probe        
            [probe] = utils.apply_3D_apodization(probe, 0); 

            probe = probe .* sqrt(1e5/sum(sum(abs(probe).^2)));
            clear x y r2 w probe_hr  

        else
            verbose(2, 'Using circular pinhole as initial model.');
            [x1,x2] = ndgrid(-asize(1)/2:floor((asize(1)-1)/2),-asize(2)/2:floor((asize(2)-1)/2));
            probe = ( (x1 * dx_spec(1)).^2 + (x2 * dx_spec(2)).^2 < (p.model.probe_diameter/2)^2);
            probe = prop_free_nf(double(probe), lambda, p.model.probe_propagation_dist, dx_spec);
            clear x1 x2
        end
    end
    verbose(3, 'Successfully generated model probe.');
else
    if ~isfield(p,'probe_file_propagation')
        p.probe_file_propagation = [];
    end
    verbose(2, 'Using previous run as initial probe.');
    
    % if string allows it, fill in the scan numbers 
    p.initial_probe_file = sprintf(replace(p.initial_probe_file,'\','\\'), p.scan_number(1)); 

    for searchpath = {'', p.ptycho_matlab_path}
        fpath = dir(fullfile(searchpath{1},p.initial_probe_file)); 
        % check if only one unique file is found
        if length(fpath) > 1
            error('Too many paths corresponding to patterns %s were found', p.initial_probe_file)
        elseif length(fpath) == 1
            p.initial_probe_file = fullfile(fpath.folder, fpath.name);
            break
        end
    end
    if isempty(fpath)
        error(['Did not find initial probe file: ' p.initial_probe_file])
    end
    
    fileokflag = 0;
    while ~fileokflag
        try
            S = load_ptycho_recons(p.initial_probe_file, 'probe'); % avoid object loading when it is not needed 
            probe = S.probe;
            probe = probe(:,:,:,1); %%added by YJ. Force to ignore the 4-th dimension (used for storing OPR modes)
            %disp(size(probe))
            S = load_ptycho_recons(p.initial_probe_file, 'p');
            fileokflag = 1; 
            verbose(2, 'Loaded probe from: %s',p.initial_probe_file );

            %% check if the loaded probe was binned or no
            if isfield(S, 'p') && isfield(S.p, 'binning')
                binning = S.p.binning;
            elseif isfield(S, 'p') && isfield(S.p, 'detector') && isfield(S.p.detector, 'binning')
                binning = S.p.detector.binning;
            else
                if verbose() > 0
                    binning = []; 
                    while isempty(binning)
                        binning = str2num(input('Define binning factor 2^x for loaded initial probe (i.e. 0 for no binning):','s'));
                    end
                    % save provided binning option to the loaded probe file
                    S.p.detector.binning = binning; 
                    save(p.initial_probe_file, '-append', '-struct', 'S')
                else
                    % prevent stopping code if automatic reconstructions are running
                    verbose(0, 'Initial probe binning could not be determined, assuming no binning')
                    binning = 0; 
                end
            end
            % modify the loaded probe into a nonbinned version
            probe = crop_pad(probe, [size(probe,1),size(probe,2)]*2^binning);            
        catch err
            disp(['File corrupt: ' p.initial_probe_file])
            disp(err.message)
            disp('Retrying')
            
            pause(1)
        end
    end
    if ndims(probe)==3
        sz_pr = size(probe);
        probe = reshape(probe, [sz_pr(1) sz_pr(2) 1 sz_pr(3)]);
    end
    verbose(3, 'File %s loaded successfully.', p.initial_probe_file);

    if ~all([size(probe,1) size(probe,2)] == asize)
        verbose(2,'Loaded probe has the wrong size.');
        if isfield(p,'crop_pad_init_probe') && p.crop_pad_init_probe %added by YJ
            verbose(2,'Crop/pad probe in file %s, from (%d,%d) to (%d,%d).', p.initial_probe_file,size(probe,1),size(probe,2),asize(1),asize(2));
            probe = crop_pad(probe, asize);            
        else
            verbose(2,'Interpolating probe in file %s, from (%d,%d) to (%d,%d).', p.initial_probe_file,size(probe,1),size(probe,2),asize(1),asize(2));
            probe = interpolateFT(probe,asize);
        end
    end
    if ~isempty(p.probe_file_propagation) && any(p.probe_file_propagation ~= 0)
        verbose(2,'Propagating probe from file by %f mm',p.probe_file_propagation*1e3);
        probe= prop_free_nf(double(probe), lambda, p.probe_file_propagation, dx_spec);
    end

end
if isfield(p,'normalize_init_probe') %%added by YJ
    if p.normalize_init_probe
        probe = probe .* sqrt(a2 ./ sum(sum(abs(probe).^2)));
    end
else
	probe = probe .* sqrt(a2 ./ sum(sum(abs(probe).^2)));
end
pout = p;
pout.probe_initial = probe;
