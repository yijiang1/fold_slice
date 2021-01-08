% Prepare probe - Only single file supported, either the file has a 3D matrix
% of many probes or the probe will be repeated for the number of probes
% needed for reconstruction

function [ p ] = prepare_initial_probes( p )

    p = core.ptycho_model_probe(p);
    
    p.probes = double(p.probe_initial);
    if size(p.probes,3) ~= p.numprobs        
        p.probes = repmat(p.probes,[1 1 p.numprobs]);
    end
    % If prepared without modes, but reconstruction needs modes
    % allocate initial guess.
    if size(p.probes,4) ~= p.probe_modes
        % Determine mode energies
        Emod = zeros(p.probe_modes,1);
        for jj = 1:numel(Emod)-1
            if jj <= numel(p.mode_start_pow)
                Emod(jj+1) = p.mode_start_pow(jj);
            else
                Emod(jj+1) = p.mode_start_pow(end);
            end
        end
        %             if (numel(p.mode_start_pow) == 1)||(numel(p.mode_start_pow) == p.probe_modes-1)
        %                 Emod(2:end) = p.mode_start_pow;
        if sum(Emod) > 1
            error('Energy distribution between modes exceeds 1, see p.mode_start_pow')
        else
            Emod(1) = 1-sum(Emod);
        end
        Emod_init = Emod;
        %disp(p.numprobs)
        for prnum = 1:p.numprobs
            % Determine the total energy of the probe first mode
            aux = p.probes(:,:,prnum,1);
            Etot = sum(abs(aux(:)).^2);
            Emod = Emod_init*Etot; % Now Emod has really the expected total sum
            
            if strcmpi(p.mode_start,'rand')
                for prmode = 2:p.probe_modes
                    p.probes(:,:,prnum,prmode) = p.probes(:,:,prnum,1).*(2*rand(p.asize)-1);
                end
            elseif strfind(p.mode_start,'herm')
                if strcmpi(p.mode_start,'herm')
                    M = ceil(sqrt(p.probe_modes))-1;
                    N = ceil(p.probe_modes/(M+1))-1;
                elseif strcmpi(p.mode_start,'hermver')
                    M = 0;
                    N = p.probe_modes-1;
                elseif strcmpi(p.mode_start,'hermhor')
                    M = p.probe_modes-1;
                    N = 0;
                else
                    error('Unknown p.mode_start')
                end
                x = [1:size(p.probes,2)]-size(p.probes,2)/2;
                y = [1:size(p.probes,1)]-size(p.probes,1)/2;
                [X Y] = meshgrid(x,y);
                
                H = core.hermite_like(squeeze(p.probes(:,:,prnum,1)),X,Y,M,N);
                if prnum == 1
                    p.probes(:,:,:,2:p.probe_modes) = 0;
                end
                p.probes(:,:,prnum,2:p.probe_modes) = reshape(H(:,:,2:p.probe_modes),size(p.probes(:,:,prnum,2:p.probe_modes)));
            else
                error('Undefined p.mode_start')
            end
            % Normalization
            for prmode = 1:p.probe_modes
                p.probes(:,:,prnum,prmode) = p.probes(:,:,prnum,prmode)*sqrt(Emod(prmode)/(sum(sum(abs(p.probes(:,:,prnum,prmode)).^2))));
            end
            
            %  p.   mode_start_pow = [0.02] ;   % Integrated intensity on modes. Can be a number (all modes equal) or a vector
            %  p.   mode_start = 'rand'         % (for probe) = 'rand', = 'her' (Hermitian-like base), = 'herver' (vertical modes only), = 'herhor' (horizontal modes only)
            %  p.   mode_her_ord = [];          % (for probe) Specify a 2xn vector with the (m,n) starting orders, leave = [] for default
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*(1+0.01.*rand(p.asize))*0.5;
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*exp(1i*0.1*pi.*rand(p.asize))*0.1;
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*(2*rand(p.asize)-1)*5;
        end
    end
    
    %Added by YJ. Force orthogonalization of initial probes
    if isfield(p,'ortho_init_probes') && p.ortho_init_probes
        % Determine mode energies
        Emod = zeros(p.probe_modes,1);
        for jj = 1:numel(Emod)-1
            if jj <= numel(p.mode_start_pow)
                Emod(jj+1) = p.mode_start_pow(jj);
            else
                Emod(jj+1) = p.mode_start_pow(end);
            end
        end
        %             if (numel(p.mode_start_pow) == 1)||(numel(p.mode_start_pow) == p.probe_modes-1)
        %                 Emod(2:end) = p.mode_start_pow;
        if sum(Emod) > 1
            error('Energy distribution between modes exceeds 1, see p.mode_start_pow')
        else
            Emod(1) = 1-sum(Emod);
        end
        Emod_init = Emod;
        for prnum = 1:p.numprobs
            % Determine the total energy of the probe first mode
            aux = p.probes(:,:,prnum,1);
            Etot = sum(abs(aux(:)).^2);
            Emod = Emod_init*Etot; % Now Emod has really the expected total sum
            
            if strcmpi(p.mode_start,'rand')
                for prmode = 2:p.probe_modes
                    p.probes(:,:,prnum,prmode) = p.probes(:,:,prnum,1).*(2*rand(p.asize)-1);
                end
            elseif strcmpi(p.mode_start,'zeros') %added by YJ
                for prmode = 2:p.probe_modes
                    p.probes(:,:,prnum,prmode) = ones(p.asize)*eps;
                end
            elseif strfind(p.mode_start,'herm')
                if strcmpi(p.mode_start,'herm')
                    M = ceil(sqrt(p.probe_modes))-1;
                    N = ceil(p.probe_modes/(M+1))-1;
                elseif strcmpi(p.mode_start,'hermver')
                    M = 0;
                    N = p.probe_modes-1;
                elseif strcmpi(p.mode_start,'hermhor')
                    M = p.probe_modes-1;
                    N = 0;
                else
                    error('Unknown p.mode_start')
                end
                x = [1:size(p.probes,2)]-size(p.probes,2)/2;
                y = [1:size(p.probes,1)]-size(p.probes,1)/2;
                [X Y] = meshgrid(x,y);
                H = core.hermite_like(squeeze(p.probes(:,:,prnum,1)),X,Y,M,N);
                %disp(size(H))

                if prnum == 1
                    p.probes(:,:,:,2:p.probe_modes) = 0;
                end
                p.probes(:,:,prnum,2:p.probe_modes) = reshape(H(:,:,2:p.probe_modes),size(p.probes(:,:,prnum,2:p.probe_modes)));
            else
                error('Undefined p.mode_start')
            end
            % Normalization
            for prmode = 1:p.probe_modes
                p.probes(:,:,prnum,prmode) = p.probes(:,:,prnum,prmode)*sqrt(Emod(prmode)/(sum(sum(abs(p.probes(:,:,prnum,prmode)).^2))));
            end
            
            %  p.   mode_start_pow = [0.02] ;   % Integrated intensity on modes. Can be a number (all modes equal) or a vector
            %  p.   mode_start = 'rand'         % (for probe) = 'rand', = 'her' (Hermitian-like base), = 'herver' (vertical modes only), = 'herhor' (horizontal modes only)
            %  p.   mode_her_ord = [];          % (for probe) Specify a 2xn vector with the (m,n) starting orders, leave = [] for default
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*(1+0.01.*rand(p.asize))*0.5;
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*exp(1i*0.1*pi.*rand(p.asize))*0.1;
            %   p.probes(:,:,:,prmode) = p.probes(:,:,:,1).*(2*rand(p.asize)-1)*5;
        end
    end

end

