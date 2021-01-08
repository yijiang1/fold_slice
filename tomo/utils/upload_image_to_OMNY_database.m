% UPLOAD_IMAGE_TO_OMNY_DATABASE upload images of the reconstruction to image gallery tomography database
%
% upload_image_to_OMNY_database(tomogram_delta, par)
%
% Inputs
%   **tomogram_delta      delta tomogram 
%   **par                 tomography parameter structure 

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
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



function upload_image_to_OMNY_database(tomogram_delta, par)

    % save snapshots for online viewer 
    plotting.smart_figure(1)
    par.displayslice = [];

    par.displayaxis = 3;      
    tomo.show_tomogram_cuts(tomogram_delta, par.scanstomo, par)
    path = sprintf('%s_xy_tomo.png',par.online_tomo_path); 
    print('-f1','-djpeg','-r300',path); 
    system(sprintf('convert -trim %s %s', path, path));
    system(sprintf('cp %s %s', path, sprintf('%s/preview_xy_tomo.png',par.output_folder)));

    
    par.displayaxis = 1;           
    tomo.show_tomogram_cuts(tomogram_delta, par.scanstomo, par)
    path = sprintf('%s_xz_tomo.png',par.online_tomo_path); 
    print('-f1','-djpeg','-r300',path); 
    system(sprintf('convert -trim %s %s', path, path));
    system(sprintf('cp %s %s', path, sprintf('%s/preview_xz_tomo.png',par.output_folder)));

    par.displayaxis = 2;           
    tomo.show_tomogram_cuts(tomogram_delta, par.scanstomo, par)
    path = sprintf('%s_yz_tomo.png',par.online_tomo_path);
    print('-f1','-djpeg','-r300',path);    
    system(sprintf('convert -trim %s %s', path, path));
    system(sprintf('cp %s %s', path, sprintf('%s/preview_yz_tomo.png',par.output_folder)));

    
    path = fullfile(par.output_folder,'Database_xy_tomo.png');
    print('-f1','-dpng','-r300',path);
    system(sprintf('convert -trim "%s" "%s"', path, path));

    utils.verbose(0, 'Uploading images to OMNY database')
    
    unix_cmd = sprintf('/work/sls/spec/local/XOMNY/bin/upload/upload_tomography_slice.sh %d %s',par.tomo_id, path);
    utils.verbose(0,'Uploading to database %s',path)
    fprintf('%s\n',unix_cmd)
    unix(unix_cmd);
    delete('upload.php')  % remove some leftover from the upload_tomography_slice
    
end 
