 function [im_out,range]=imagesc_tomo( varargin )    
    % Function [im_out]=imagesc_tomo( varargin )
    % Parameters: 
    %   data - 3D array to be shown 
    %   colormap - name of standard matlab colormap, default = bone 
    %   clim - range of the colorbar, default is [] (auto range)
    %   axis - cell of strings with axis parameters, eg. {'image', 'off'}

    import math.sp_quantile
    import plotting.imagesc3D
    

    par = inputParser;
    par.addOptional('data', [])
    par.addParameter('colormap',  'bone' )
    par.addParameter('clim',  [] , @isnumeric )
    par.addParameter('axis',  {'image', 'off'}, @iscell)
    par.addParameter('colorbar',false, @islogical)

    par.parse(varargin{:})
    r = par.Results;
    data = r.data;
    
    if ismatrix(data) 
       imagesc3D(data);
       return
    end
        
    Npix = size(data);
    cntr = ceil(Npix/2);
        
    if isa(data, 'gpuArray')
        data = gather(data); % download data from GPU 
    end
    
    if isempty(r.clim)
        range = sp_quantile(data, [1e-3, 1-1e-3], max(1, ceil(sqrt(numel(data))/1e3)));
    else
        range = r.clim;
    end
    %disp('range:')
    %disp(range)
    
    colormap(r.colormap)
           
    args = {'show_play_button', false, 'show_edit_box', false };
    
    ax(1) = subplot(2,3,1);
    imagesc3D(rot90(permute(data, [2,3,1]),1), 'init_frame', cntr(1), 'slider_position',slider_pos(ax(1)),  args{:}); 
    axis(r.axis{:});
    if isreal(data) && range(1) < range(2); caxis(range); end
    if r.colorbar; colorbar; end
    title('Centralslice xz');
    ax(2) = subplot(2,3,2);
    imagesc3D(rot90(permute(data, [1,3,2]),1), 'init_frame', cntr(2),'slider_position',slider_pos(ax(2)), args{:});
    axis(r.axis{:});
    if isreal(data) && range(1) < range(2); caxis(range); end
    if r.colorbar; colorbar; end
    title('Centralslice yz');
    ax(3) = subplot(2,3,3);
    imagesc3D(data, 'init_frame', cntr(3),'slider_position',slider_pos(ax(3)), args{:});
    axis(r.axis{:});
    if isreal(data) && range(1) < range(2); caxis(range); end
    if r.colorbar; colorbar; end
    title('Centralslice xy');
    ax(4) = subplot(2,3,4);
    imagesc3D(rot90(squeeze(sum(data,1)),1));
    axis(r.axis{:});
    title('Projection xz');
    if r.colorbar; colorbar; end
    ax(5) = subplot(2,3,5);
    imagesc3D(rot90(squeeze(sum(data,2)),1));
    axis(r.axis{:});
    title('Projection yz');
    if r.colorbar; colorbar; end
    ax(6) = subplot(2,3,6);
    imagesc3D(squeeze(sum(data,3)));
    axis(r.axis{:});
    title('Projection xy');
    if r.colorbar; colorbar; end
    
    % link together axis from the projections and the slices 
    linkaxes(ax([1,4]), 'xy')       
    linkaxes(ax([2,5]), 'xy')
    linkaxes(ax([3,6]), 'xy')
    
    if nargout > 0
       im_out = ax; 
    end

end

function slider_default = slider_pos(ax)
    pos = ax.Position;
    slider_default = [pos(1)+pos(3)/2-0.095 pos(2)-0.045 0.2 0.04];
end


