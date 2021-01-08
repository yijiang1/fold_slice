%  SAVE_MOVIE Function to save stack of images as a movie 
%
%  save_movie(img_stack, movie_filename, theta, par, varargin)
% 
%  ** img_stack    stack of input images, real or complex 
%  ** movie_filename   name of the movie to be saved 
% *optional*
%  ** theta   angles of the projections, default = []
%  ** par     structure with tomo parameters, default = [] 
%  ** varargin   other parameters, see the code for more details 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or pareters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



function save_movie(img_stack, movie_filename, theta, par, varargin)

    import math.*

    parser = inputParser;
    parser.addParameter('pixel_size',  nan , @isnumeric )
    parser.addParameter('colormap',  bone , @isnumeric )
    parser.addParameter('output_folder',  '' , @isstr )
    parser.addParameter('fps',  10 , @isnumeric )
    parser.addParameter('quality',  100 , @isnumeric ) % 100 is maximal 
    parser.addParameter('baraxis','auto') 
    parser.addParameter('windowautopos',false) 

    parser.parse(varargin{:})
    r = parser.Results;

    if nargin < 4
        par = struct(); 
    end
    
    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    movie_path = [par.output_folder, movie_filename]; 
    
    if exist(movie_path,'file')
        disp(['File ' movie_path ' exists,' ])
        userans = input(['Do you want to overwrite (Y/n)? '],'s');
        if strcmpi(userans,'n')
            disp(['Did not save ' par.saveprojfile])
            return
        end
    else
        delete(movie_path);
    end

    Nslices = size(img_stack,3);

    
    frames = 1:Nslices;
    if ~isempty(theta) && par.showsorted 
        [~, ind] = sort(theta);
        frames = ind(frames);
    end
    
    
    clf()
    if isfield(par, 'baraxis') && strcmpi(par.baraxis,'auto') && isreal(img_stack)
        % set the same range for all the frames using quantile range
        range = gather(sp_quantile(img_stack, [1e-3, 1-1e-3], 10));
    end
    
     
    
    
    
    disp(['Saving movie to  ' movie_path]);
    writeobj = VideoWriter(movie_path);
    writeobj.Quality=par.quality;
    writeobj.FrameRate=par.fps;
    
    open(writeobj); 

    % Create an animation.
    plotting.imagesc3D(img_stack(:,:,1))
    axis off image 
    colormap(par.colormap)

    set(gca,'nextplot','replacechildren');

    if  par.windowautopos
        screensize = get( groot, 'Screensize' );
        win_size = [1060 767]; 
        set(gcf,'Outerposition',[141 min(257,screensize(2)-win_size(2)) 1060 767]);
    end

    
    
    disp('Creating movie')
    for kk = 1:Nslices
       utils.progressbar(kk, Nslices)
       % use imagesc3D to image also complex valued images 
       plotting.imagesc3D(img_stack(:,:,frames(kk)))
       if isreal(img_stack) && isfield(par, 'baraxis') 
           if ~strcmpi(par.baraxis,'auto')
               caxis(par.baraxis)
           else
               caxis(range)
           end
       end
       axis xy
       % Write each frame to the file.
       currFrame = getframe;
       writeVideo(writeobj,currFrame);
    end
    
    close(writeobj); 
    
    disp('Movie finished')



end