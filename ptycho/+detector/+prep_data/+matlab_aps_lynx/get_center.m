%GET_CENTER Estimate the center of the diffraction pattern

function [ p ] = get_center( p )
import utils.verbose
import io.*

detStorage = p.detectors(p.scanID).detStorage;
ctr = detStorage.ctr;
%{
if verbose>2
    
    verbose(3, 'Loading sample image.')
    % read sample image
    if p.detectors(p.scanID).params.data_stored   % if the data should be loaded from disk
        if verbose > 2
            [p] = p.detectors(p.scanID).params.get_filename(p);
            files = detStorage.files;
            if numel(files) == 0
                error('Did not find any files.')
            end
            det.params = p.detectors(p.scanID).params;
            
            sample_file = image_read(files{1}, det.params.image_read_extraargs);
            
        end
    end
    
    sz = size(sample_file.data(:,:,1));
    f = double(sample_file.data(:,:,1));
    % Look for center
    [~, cy] = max(sum(f.*detStorage.mask,2));
    [~, cx] = max(sum(f.*detStorage.mask,1));
    ctr_auto = [cy, cx];
    
    if strcmp(detStorage.check_ctr, 'auto')
        ctr = ctr_auto;
        verbose(3, sprintf('Using center: (%d, %d)', ctr(1), ctr(2)));
    elseif strcmp(detStorage.check_ctr, 'inter')
        imagesc(log(detStorage.f));
        [cx,cy] = getpts;
        close(gcf);
        ctr = round([cy, cx]);
        verbose(3, sprintf('Using center: (%d, %d) - I would have guessed it is (%d, %d)', ctr(1), ctr(2), ctr_auto(1), ctr_auto(2)));
    else
        verbose(3, sprintf('Using center: (%d, %d) - I would have guessed it is (%d, %d)', ctr(1), ctr(2), ctr_auto(1), ctr_auto(2)));
    end
else
    verbose(2, sprintf('Using center: (%d, %d)', ctr(1), ctr(2)));
end
%}
verbose(2, sprintf('Using center: (%d, %d)', ctr(1), ctr(2)));
detStorage.ctr = ctr;

if ~p.prealign_FP
    detStorage.lim_inf = ctr-p.asize/2;
    detStorage.lim_sup = ctr+p.asize/2-1;
else
    detStorage.lim_inf = ctr-p.prealign.asize/2;
    detStorage.lim_sup = ctr+p.prealign.asize/2-1;
end

verbose(2, sprintf('Selected region: (''RowFrom'', %d, ''RowTo'', %d, ''ColumnFrom'', %d, ''ColumnTo'', %d)', detStorage.lim_inf(1), detStorage.lim_sup(1), detStorage.lim_inf(2), detStorage.lim_sup(2)));
%{
if verbose>2
    if any(detStorage.lim_inf < 1) || any(detStorage.lim_sup > sz)
        error('Array size exceeds limit (according to position of center, should be < %d)', max([1-detStorage.lim_inf, detStorage.lim_sup-sz]));
    end
end
%}
end

