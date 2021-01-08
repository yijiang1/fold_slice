% SHOW_TOMOGRAM_CUTS show cuts through the reconstructed volume 
%
%     show_tomogram_cuts(tomogram, scanstomo, par, extra_string = '' )
%
% Inputs:
%  **tomogram       - reconstructed volume 
%  **scanstomo      - scan numbers, only for naming 
%  **par            - parameter structure 
%  **extra_string   - string added to the saved name , default = ''
    
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: "Data processing was carried out 
%   using the "cSAXS matlab package" developed by the CXS group,
%   Paul Scherrer Institut, Switzerland." 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided "as they are" without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



function show_tomogram_cuts(tomogram, scanstomo, par, extra_string)
import math.*

if nargin < 4
    extra_string = '';
end


if isa(tomogram, 'gpuArray')
    tomogram = gather(tomogram); 
end
    
if par.makemovie    % Open movie file
    movie_filename = fullfile(par.output_folder,['tomo_movie_', par.scale '_' par.scans_string '_' extra_string ...
        '_movie_axis_' sprintf('%01d',par.displayaxis) '.avi']);

    if exist(movie_filename,'file')
        disp(['File ' movie_filename ' exists,' ])
        userans = input('Do you want to overwrite (y/N)? ','s');
        if strcmpi(userans,'y')
            utils.verbose(0,['Saving movie to  ' movie_filename]);

        else
            utils.verbose(0,['Did not save ' movie_filename])
            return
        end
    else
        utils.verbose(0,['Saving movie to  ' movie_filename]);
    end
    writeobj = VideoWriter(movie_filename); 
    writeobj.Quality=90;
    writeobj.FrameRate=5; 
    open(writeobj); 
end

% If displayslices is empty show central slice
if isempty(par.displayslice)&&(~par.animatedslices)
    utils.verbose(1,'Displaying central slice along axis %i', par.displayaxis)
    par.displayslice = round(size(tomogram,par.displayaxis)/2);
end

par.displayslice = unique(max(1,min(size(tomogram,par.displayaxis),round(par.displayslice))));

% Determine range of tomogram
switch num2str(par.tomobaraxis)
case 'auto_per_frame'
    autobar = true;
    slices_ind = {':', ':', ':'}; 
    slices_ind{par.displayaxis} = par.displayslice; 
    par.tomobaraxis = sp_quantile(tomogram(slices_ind{:}), [1e-4, 1-1e-4],5); 
case 'auto'
    autobar = true;
    % ignore outliers 
    par.tomobaraxis = sp_quantile(tomogram, [1e-4, 1-1e-4],ceil(max(10, sqrt(numel(tomogram))/100)));
    % full range 
    %par.tomobaraxis = [min(tomogram(:), max(tomogram(:))];
otherwise
    autobar = false;
end
switch lower(par.scale)
    case 'phase'
        if autobar
            par.tomobaraxis = par.tomobaraxis/par.factor;
            tomogram = tomogram / par.factor;
        end
        strscale = 'phase';
    case 'delta'
        strscale = 'delta';
    case 'edensity'
        if autobar
            par.tomobaraxis = sort(par.tomobaraxis*par.factor_edensity);
        end
        strscale = 'electron density [e/A^3]';
    case 'amp'
        strscale = 'amplitude';
    case 'beta'
        strscale = 'beta';
    case ''
        strscale = '';
    otherwise
            error('scale should be phase, delta, amp, beta or edensity')
end

%%% Here the option for showing animation
if par.average_slices == 1
    par.animatedslices = 0; 
end
if par.animatedslices
    par.displayslice = [1:size(tomogram,par.displayaxis)];
end
   
if par.average_slices == 0
   loopdisplayslice = par.displayslice;
elseif par.average_slices == 1
   loopdisplayslice = 1;    
end

fig = plotting.smart_figure(1); 
clf()
if par.windowautopos
    screensize = get( 0, 'Screensize' );
    set(gcf,'Outerposition',[1 screensize(4)-650 640 665]);
    par.windowautopos = false;
end
rect = get(fig,'Position'); 
rect(1:2) = [0 0]; 

for showslice = loopdisplayslice
    % Determine sagital, coronal or axial slices
    slice = {':',':',':'};
    if loopdisplayslice==1
        slice{par.displayaxis} = par.displayslice;
    else
        slice{par.displayaxis} = showslice;
    end
    if showslice > size(tomogram,par.displayaxis)
        continue
    end
    sliceview = squeeze(mean(tomogram(slice{:}),par.displayaxis))';
    if par.displayaxis == 3
        sliceview = sliceview';
    end
    sectionstring = {'Coronal','Sagital',  'Axial'};
    sectionstring = sectionstring{par.displayaxis};
   
    
    switch lower(par.scale)
        case 'phase'
            sliceview = sliceview/par.factor;
        case 'delta'
            
        case 'edensity'
            sliceview = sliceview*par.factor_edensity;
        case 'amp'
        case 'beta'
        case ''
        otherwise
            error('scale should be phase, delta, amp, beta or edensity')
    end
    

    if ~par.realaxis
        imagesc(sliceview)
    else
        xaux = ([1 size(sliceview,2)]-size(sliceview,2)/2)*par.pixel_size*1e6;
        yaux = ([1 size(sliceview,1)]-size(sliceview,1)/2)*par.pixel_size*1e6;
        imagesc(xaux,yaux,sliceview)
        xlabel('microns')
        ylabel('microns')
    end
    axis xy image
    c = colormap(par.colormapchoice);
    if par.reverse_contrast
        c = flipud(c); 
        colormap(c); 
    end
    caxis(sort(par.tomobaraxis))
    h = colorbar;
    ylabel(h, strscale)
    
    if (~isempty(par.bar_length))&&par.realaxis  %% Show scale bar
        hold on
        axisaux = axis;  
        rectangle('Position', [axisaux(1)+par.bar_start_point(1)*1e6 axisaux(3)+par.bar_start_point(2)*1e6 par.bar_length*1e6 par.bar_height*1e6], ...
                    'facecolor',par.bar_color,'edgecolor','none') 
        text(axisaux(1)+par.bar_start_point(1)*1e6,...
            axisaux(3)+par.bar_start_point(2)*1e6+par.bar_height*2e6,...
            [num2str(par.bar_length*1e6) ' microns'],'Color',par.bar_color,'FontSize',12);
        hold off
    end
    if par.average_slices == 1
        title(strrep(sprintf(['Tomogram ' strscale ': ' par.scans_string, ...
        ' ' sectionstring ' section: \n Average slices ' num2str(par.displayslice(1)) ' to ' num2str(par.displayslice(end))]),'_', '\_'))
    else    
        title(['Tomogram ' strscale ': ' strrep(par.scans_string, '_', '\_') ...
            ' ' sectionstring ' section: Slice ' num2str(showslice)])
    end

    drawnow
    if par.makemovie
        currFrame = getframe(fig,rect);
        writeVideo(writeobj,currFrame);
    end
    pause(par.pausetime)
    

end


if par.makemovie == 1 
    close(writeobj); 
end

if par.writesnapshots   && ~debug()
    output_path = fullfile(par.output_folder,['tomo_cut_', par.scans_string '_' par.scale '_' extra_string '_' num2str(size(sliceview,1)) 'x' num2str(size(sliceview,2)) '_axis_' num2str(par.displayaxis)]);

    if par.average_slices == 1
        output_path = [output_path, '_average_slices_' num2str(par.displayslice(1)) '_to_' num2str(par.displayslice(end))];
    else
        output_path = [output_path, '_slice_' num2str(showslice)];
    end
    fprintf('Writting image files \n %s.png \n %s.eps\n',output_path,output_path);
    print('-f1','-dpng','-r300',[output_path,'.png']);
    print('-f1','-depsc2',[output_path,'.eps']);
end

end
