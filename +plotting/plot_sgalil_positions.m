import io.image_read

% this script  is used to visualize  the sgalil positions
clear
scan_init=11103;
number_of_scans = 10;

for jj=1:number_of_scans
this_scan=scan_init+jj-1;
filename = sprintf('~/Data10/sgalil/S%05d.dat',this_scan);

data = beamline.read_position_file(filename);
x = data.Avg_x;
y = data.Avg_y;

pos_arrayx(:,jj) = x;
pos_arrayy(:,jj) = y;

end
plot(pos_arrayx,pos_arrayy,'-bo')
grid on
axis equal
title(sprintf('S%05d - S%05d',scan_init,scan_init+number_of_scans))


pos_arrayy_aux = pos_arrayy;
pos_arrayy_aux(:,2:2:end) = flipud(pos_arrayy_aux(:,2:2:end));
offset = mean(mean(abs(diff(pos_arrayy_aux(1:6,1:6),1,2))));

fprintf('**************************************************************************** \n');
fprintf('The mean step size is %.05f microns \n', mean(abs(diff(y(:, 1))))*1000);
fprintf('The standard deviation is %.05f microns \n', std(abs(diff(y(:, 1))))*1000);
fprintf('The offset between lines (absolute value) is is %.05f microns \n', offset*1e3);


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
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%