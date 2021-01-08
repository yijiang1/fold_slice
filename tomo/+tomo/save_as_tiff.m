% SAVE_AS_TIFF simple function to save data as tiff images 
%
%  save_as_tiff(rec, param)
%
% Inputs:
%   **rec           reconstructed volume 
%   **p             parameter  structure 
%   **extra_string  extra string added to the saved name for example name of the reconstruction method 
% Parameters 
%   params.scans_string = 'name'
%   params.save_as_stack = false; 
%   params.tiff_compression = 'none';
%   params.tiff_folder_name = 'folder'; 

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


function save_as_tiff(rec, p, extra_string)

    import utils.*
    
    tiff_folder_name = fullfile(p.output_folder,p.tiff_subfolder_name);
    if ~exist(tiff_folder_name, 'dir')
        mkdir(tiff_folder_name); 
    elseif ~isfield(p, 'force_overwrite') || p.force_overwrite == false
        display(['Folder exists: ' tiff_folder_name])
        userans = input(['Do you want to overwrite TIFFs in this folder (Y/n)? '],'s');
        if ~strcmpi(userans,'n')
            disp('Overwritting');
        else
            error('Writting of TIFFs aborted')
        end
    end

    cutoff =  [min(rec(:)),max(rec(:))];

    rec = (rec-cutoff(1))/(cutoff(2)-cutoff(1)); 
    rec_uint16 =  uint16( (2^16-1) *rec);
    
    verbose(0,['Saving to:', tiff_folder_name '/' p.name_prefix '_' p.scans_string '_' ...
                extra_string '.tif'])
            
    Nlayers = size(rec,3);
    for j=1:Nlayers
        progressbar(j, Nlayers)
        if p.save_as_stack
            image_filename_with_path = [tiff_folder_name '/' p.name_prefix '_' p.scans_string '_' ...
                p.filter_type '_freqscl_' sprintf('%0.2f',p.freq_scale) '.tif'];
            if j == 1
                imwrite(rec_uint16,image_filename_with_path,'tiff',...
                    'Compression',p.tiff_compression);
            else
                imwrite(rec_uint16,image_filename_with_path,'tiff',...
                    'Compression',p.tiff_compression,'WriteMode','append');
            end
        else
            imwrite(rec_uint16(:,:,j),[tiff_folder_name '/' p.name_prefix '_' p.scans_string '_' extra_string '_' sprintf('%04d',j) '.tif'], ...
                'tiff','Compression',p.tiff_compression);
        end
    end


    fid=fopen([tiff_folder_name '_cutoffs.txt'],'w');
    fprintf(fid, '# low_cutoff = %e\n', cutoff(1));
    fprintf(fid, '# high_cutoff = %e\n', cutoff(2));
    fprintf(fid, '# factor = %e\n', p.factor);
    fprintf(fid, '# pixel size = %e\n', p.pixel_size);
    fprintf(fid, '# factor_edensity = %e\n', p.factor_edensity);
    fprintf(fid, '# Conversion formula\n');
    fprintf(fid, '# im_delta_from_tiff = im_tiff*(high_cutoff-low_cutoff)/(2^16-1) + low_cutoff;\n');
    fprintf(fid, '# im_edensity_from_tiff = im_delta_from_tiff*factor_edensity;\n');
    fclose(fid);
end