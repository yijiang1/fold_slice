% GET_VARIABLE_PROBE calculate variable probe (at different scan positions) for given mean/var probe and its evolution inputs 
% 
% probe = get_variable_probe(probe, probe_evolution, p_ind) 
%
% ** probe     [Nx,Ny,probe_modes,variable_modes] array,  variable probe 
% ** probe_evolution   [Npos,variable_modes] array containing evolution of the varaible modes for each position 
% ** p_ind      indices containg corresponding probe id for each processed position
%
% returns:
% ++ probe       [Nx,Ny,N] array, a different probe for each scan position 
%
% see also: engines.GPU.LSQML 


function probe = get_variable_probe(probe, probe_evolution, p_ind)
 
    Np_p = [size(probe,1),size(probe,2)];

    probe = probe(:,:,p_ind,:); 
    
    if length(p_ind) == 1  % in case that only single probe is called  
        probe = reshape(probe,prod(Np_p),[]);
        probe = reshape(probe * probe_evolution', Np_p(1), Np_p(2), []);

    else  % in case multiple scans with unshared probe 
        probe_out = 0; 
        for ii = 1:size(probe,4)
            probe_out = probe_out + probe(:,:,:,ii).*reshape(probe_evolution(:,ii),1,1,[]); 
        end
        probe = probe_out; 
    end
end
