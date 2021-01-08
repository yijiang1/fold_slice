%  APPLY_BINNING apply binning / upsampling on all relevant parameters except data and mask
%
% p = apply_binning(p, bin_data)
%
% ** p              p structure
% ** binning        if binning > 1, then data are binned , if binning < 1, data will be upsampled 
% returns:
% ++ p          p structure
% 

function p = apply_binning(p, binning)
    % apply binning / upsampling on all relevant parameters except data and mask
    if binning > 0
        assert(all(rem(p.asize, binning)==0), 'Array size cannot be divided for binning')
    end

    % Modify variables for binning
    p.ds = p.ds*binning;
    
    if check_option(p,'prop_regime', 'farfield')
        p.object_size = p.object_size + ( 1/binning-1)*p.asize ;
        for ii = 1:p.numobjs
            p.object{ii} = utils.crop_pad(p.object{ii},p.object_size(ii,:));  % crop_pad is better when if the binned reconstruction is loaded from file as an initial guess 
        end

        p.asize = p.asize/binning;


        %% always assume that no binning was applied on the provided probes
        p.probe_initial = utils.crop_pad( p.probe_initial, p.asize);
        p.probes = utils.crop_pad( p.probes, p.asize);
    else
        p.object_size = ceil(p.object_size / binning);
        for ii = 1:p.numobjs
            p.object{ii} = utils.interpolateFT(p.object{ii},p.object_size(ii,:));  % crop_pad is better when if the binned reconstruction is loaded from file as an initial guess 
        end

        p.asize = p.asize/binning;


        %% always assume that no binning was applied on the provided probes
        p.probe_initial = utils.interpolateFT( p.probe_initial, p.asize);
        p.probes = utils.interpolateFT( p.probes, p.asize);
        
        p.dx_spec = p.dx_spec * binning; 
        p.positions = p.positions / binning; % in nearfield the pixel size also changes 
        
    end

end