%  GET_OPTIMAL_LSQ_STEP calculate the optimal step lenght for given update directions and chi array  
%
% [beta_probe, beta_object] = get_optimal_LSQ_step(self,chi,dO,dP,O,P, p_ind, par)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** chi       [Nx,Ny,N] array, difference between original and updated exit-wave  
% ** dO        [Nx,Ny,N] array, object update direction 
% ** dP        [Nx,Ny,N] array, probe update direction 
% ** O         [Nx,Ny,N] array,  object views 
% ** P         [Nx,Ny,1] or [Nx,Ny,N] array,  single or variable probe 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** par       structure containing parameters for the engines 
%
% returns:
% ++ beta_probe    optimal probe step 
% ++ beta_object   optimal object step 
%
% see also: engines.GPU.LSQML 


function  [beta_probe, beta_object] = get_optimal_LSQ_step(self,chi,dO,dP,O,P, p_ind, par)
    
    % find optimal step in the LSQ sense 

    import engines.GPU.GPU_wrapper.*
    import math.*
    global gpu use_gpu
    if ~isempty(gpu); wait(gpu); end 
   
    grouping = size(chi,3); 
    lambda_0 =  eps(single(1)) / prod(self.Np_p);
    
    lambda_LSQ = 0.1;  % add some small regularization to avoid instabilities

    
if use_gpu 
    % fast mex based CUDA version 
    if size(dP,3) == grouping 
        p_ind = 1:grouping;   % one update for each position 
    elseif size(dP,3) == 1 && numel(p_ind) == 1
        p_ind = ones(size(chi,3),1); % use only the one update given 
    elseif max(p_ind) <= size(dP,3) && numel(p_ind) == grouping
        
    else
        warning('Checkme, untested option')
        keyboard        
    end
    try
        [AA, Atb] = get_LSQ_step_mex(chi,dO,dP,O,P,lambda_0, uint8(p_ind)); 
   catch err 
        if any(strcmp(err.identifier, { 'MATLAB:UndefinedFunction','MATLAB:mex:ErrInvalidMEXFile'}))
            path = fullfile(replace(mfilename('fullpath'), mfilename, ''), 'private'); 
            mexcuda('-output', [path, '/get_LSQ_step_mex'], [path, '/get_LSQ_step_mex.cu'])

            [AA, Atb] = get_LSQ_step_mex(chi,dO,dP,O,P,lambda_0, uint8(p_ind)); 
        else
           rethrow(err) 
        end
    end

    AA = Ggather(AA); 
    Atb = Ggather(Atb);
    % is seems to be faster to get it first from GPU and then apply some
    % oprations because the matrices are too small 
    AA = sum(AA,4); 
    Atb = sum(Atb,4); 
    AA = AA + lambda_LSQ*diag(diag(mean(AA,3)));
    % the system of equations is so small that solving on CPU is good enough 
    [x1, x2] = solve_LSQ(AA(1,1,:), AA(2,1,:), AA(1,2,:), AA(2,2,:), Atb(1,1,:), Atb(2,1,:)); 
    LSQ_step = cat(1, x1, x2);
else
    
    
    if ~( par.share_probe || length(unique(p_ind)) == 1 )
        % in case of multiple scans !! 
        % replicate the update back to the original dP size 
        dP = dP(:,:,p_ind);
    end
    tic
    % prevent ill posed inversion, ideally it should be Garray(mean(abs(AA1)+abs(AA4))/2) but it i show  
    [AA1,AA2,AA4, Atb1,Atb2] = ...
            Gfun(@get_optimal_step_lsq, chi,dO,dP,...
            O,P, lambda_0); 

    AA1 = sum2(AA1); 
    AA2 = sum2(AA2); 
    AA4 = sum2(AA4); 
    Atb1 = sum2(Atb1); 
    Atb2 = sum2(Atb2);
   
    % it seems faster to solve it on GPU than using pagefun on GPU 
    AA1 = Ggather(AA1);AA2 = Ggather(AA2);AA4 = Ggather(AA4);Atb1 = Ggather(Atb1);Atb2 = Ggather(Atb2);
    AA3 = conj(AA2);
    
    lambda = 0.5;  % add some small regularization to avoid unstabilities
    I = lambda*[mean(AA1), mean(AA4)];

    AA = [ AA1+I(1), AA2; AA3, AA4+I(2)];
    Atb= [ Atb1; Atb2]; 
    
    [x1, x2] = solve_LSQ(AA(1,1,:), AA(2,1,:), AA(1,2,:), AA(2,2,:), Atb(1,1,:), Atb(2,1,:)); 
    LSQ_step = cat(1, x1, x2);

end



    LSQ_step = max(0, real(LSQ_step)); 
    LSQ_step = Ggather(LSQ_step);

    % prevent unwanted oscilation of step gets too high 
    beta_probe =  LSQ_step(2,1,:);
    beta_object = LSQ_step(1,1,:);
    
    beta_probe =   (par.beta_probe *par.beta_LSQ)* beta_probe;
    beta_object =  (par.beta_object*par.beta_LSQ)* beta_object; 
    
end

function [AA1,AA2,AA4, Atb1,Atb2] = ...
            get_optimal_step_lsq(chi,dO,dP,O,P, lambda)
     % fast kernel for estimation of optimal P and object steps 
    dOP = dO.*P;
    dPO = dP.*O;
    cdOP = conj(dOP);
    cdPO = conj(dPO);

    AA1 = real(dOP .* cdOP)+lambda;
    AA2 = (dOP .* cdPO);
    AA4 = real(dPO .* cdPO)+lambda;
    Atb1 = real(cdOP .* chi);
    Atb2 = real(cdPO .* chi);
end


function [x1, x2] = solve_LSQ(AA1, AA2, AA3, AA4, Atb1, Atb2)
    % GPU kernel to solve simple 2x2 system of equations  
    det = (AA1.*AA4 - AA2.*AA3);
    x1 = -conj(AA2.*Atb2-AA4.*Atb1) ./ det;
    x2 =  conj(AA1.*Atb2-AA3.*Atb1) ./ det;
end


