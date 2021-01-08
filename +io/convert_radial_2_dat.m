% convert_radial_2_dat converts all radial integration mat files in 
%    readpahtmask into dat files, pauses 10 minutes and repeats
% 
% Inputs:
%    **readpathmask     A cell containing the input file string masks
%    **outpathmask      A cel containint the corresponding output
%                       directories
%
% Example:
% readpathmask{1} = '~/Data10/analysis/radial_integration/*.mat';
% outpathmask{1}  = '~/Data10/analysis/radial_integration_dat/';
% readpathmask{2} = '~/Data10/analysis/radial_integration_waxs/*.mat';
% outpathmask{2}  = '~/Data10/analysis/radial_integration_waxs_dat/';
% convert_radial_2_dat(readpathmask, outpathmask)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2019 by Paul Scherrer Institute (http://www.psi.ch)    |
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


% clear

% readpathmask{1} = '~/Data10/analysis/radial_integration/*.mat';
% outpathmask{1}  = '~/Data10/analysis/radial_integration_dat/';
% 
% readpathmask{2} = '~/Data10/analysis/radial_integration_waxs/*.mat';
% outpathmask{2}  = '~/Data10/analysis/radial_integration_waxs_dat/';

function convert_radial_2_dat(readpathmask, outpathmask)

while 1==1
for ii = 1:numel(readpathmask)
    if ~exist(outpathmask{ii},'dir')
        mkdir(outpathmask{ii})
    end
    files = dir(readpathmask{ii});
    for jj = 1:numel(files)
        currentradial = fullfile(files(jj).folder,files(jj).name);
        [auxpath, auxname, auxext] = fileparts(currentradial);
        outputradial  = fullfile(outpathmask{ii},[auxname '.dat']);
        s = load(currentradial);
        
        save_data = [s.q.', ...
            reshape(s.I_all     , [size(s.I_all,1) size(s.I_all,2)*size(s.I_all,3) ]) , ...
            reshape(s.I_std     , [size(s.I_std,1) size(s.I_std,2)*size(s.I_std,3) ]) , ...
            reshape(s.norm_sum  , [size(s.norm_sum,1) size(s.norm_sum,2)*size(s.norm_sum,3) ]) , ...
            ];
        
        fprintf('Saving %s\n',outputradial);
        save( outputradial , 'save_data', '-ascii','-double');
    end
    clear files
end
fprintf('Pausing 10 minutes\n')
pause(60*10)
end
