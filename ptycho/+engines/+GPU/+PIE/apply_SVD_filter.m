% APPLY_SVD_FILTER The core of the variable probe (OPRP) code: 
% find SVD decomposition and limit the probe into several orthogonal
% modes 
% additional prior knowledge can be also included 
%
% [probe, probe_evolution] = apply_SVD_filter(probe, Nmodes, mode)
%
% ** probe     [Nx,Ny,N] variable probe for each position 
% ** Nmodes     (int)   number of variable modes 
% ** mode       structure containing parameters for selected probe mode 
% returns
% ++ probe      [Nx,Ny,variable_modes] variable probe modes
% ++ probe_evolution   [Npos,variable_modes] updated array containing evolution of the varaible modes for each position 
%
% see also: engines.GPU.PIE 

function [probe, probe_evolution] = apply_SVD_filter(probe, Nmodes, mode)

    import engines.GPU.shared.*
    import engines.GPU.GPU_wrapper.*
    import engines.GPU.shared.*

    import math.*
    import utils.* 
    
    
    Np = size(probe);
         
    [U,S,V] = fsvd(reshape((probe),[],Np(3)) ,Nmodes);
%      if any(diag(S.^2)/sum(diag(S.^2)) < 1e-3)  % 2e-3 is the weakest that FSVD can recover 
%          try
%             warning('Running full SVD (maybe use less OPR modes)  (weak modes %i/%i) ', Ggather(sum(diag(S.^2)/sum(diag(S.^2)) < 2e-3)), Nmodes)
%             [U,S,V] = svd(reshape((probe),[],Np(3)) ,0);
%             U = single(U(:,1:Nmodes));
%             S = single(S(1:Nmodes,1:Nmodes));
%             V = single(V(:,1:Nmodes));
%          catch
%              keyboard
%          end            
%      end

   %Notes by YJ:
   % U is orthonormal modes. size [Np(1)*Np(2), Nmodes]
   % S is diagonal matrix of singular values. size [Nmodes,Nmodes]
   % V is conjugated orthonormal evolution matrix. size [Npos,Nmodes]
   
   % Lower dimensional representation of the reconstructed probes = USV*
   % probe_evolution = SV*
   %disp('xxxx')
   %disp(size(U))
   U = reshape(U, Np(1),Np(2),1,[]);
   %disp(size(U))
   U = apply_probe_contraints(U, mode);
   %disp('oooo')
   
   V(:,1) = mean(V(:,1)) + 0.99*(V(:,1) - mean(V(:,1)));
   V(:,2:end) = mean(V(:,2:end)) + 0.99*(V(:,2:end) - mean(V(:,2:end)));

    %% remove outliers 
    aV = abs(V);
    MAX = quantile(aV,0.99);
    V = min(aV, MAX) .* (V ./ (aV+1e-3)); 
    %disp(size(S*V'))
    probe_evolution = (S*V').'; %Note by YJ: why transpose it here?
    %disp(size(probe_evolution))

    probe = U; 
    
    
    avg = mean(abs(probe_evolution(:,1)),1); 
    
    probe = probe*avg;
    probe_evolution = probe_evolution / avg;


end

