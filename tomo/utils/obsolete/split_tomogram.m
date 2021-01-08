% FUNCTION [ sub_thetaT sub_ind_T] = split_tomogram( indices, theta, num_sub_tomograms )
% Divides up the tomogram into sub tomograms according to the matrix
% indices.
% A combined tomogram is in each row. The single tomograms contained in each 
% row are given in the columns.
% For example, four sub tomograms out of 8 with a pairing of 1 2 , 3 4 etc
% would be given by indices = [[1 2];[3 4];[5 6];[7 8]]
% All of the first 4 sub tomograms will be returned by indices = [[1] [2] [3] [4]]
% The output are cells. 
 
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


function [ sub_thetaT sub_ind_T] = split_tomogram( indices, theta, num_sub_tomograms )

[sortedTheta, ind] = sort(theta);

%% no_subtomos = 16
    n = log2(num_sub_tomograms);
    if mod(n,1) > 0.000001
        disp(['Error - impossible number of subtomograms chosen'])
        subtomo_order = NaN;
        return
    end
    subtomos = zeros(1,num_sub_tomograms);
    k= 1;
    jj=1;
    k_prev = k;
    for i=1:n
        
        k = k_prev/2;
        while k <= 1
            jj = jj+1;
            subtomos(jj)=k;
            k=k+k_prev;
        end
        k_prev = k_prev/2;  
    end
    [dump,subtomo_order]= sort(subtomos);


%%

subtomo_order = [1 5 3 7 2 6 4 8];


for i = 1: num_sub_tomograms;
   sub_theta_u{i} = sortedTheta(subtomo_order(i):num_sub_tomograms:end);
   sub_ind_u{i} = ind(subtomo_order(i):num_sub_tomograms:end);
end

 for i = 1 : size(indices,1);
         Combined_Angle = [sub_theta_u{indices(i,1:end)}];
         Combined_Index = [sub_ind_u{indices(i,1:end)}];
         Combined = sortrows(transpose([  Combined_Index ;  Combined_Angle ]),2);
         sub_ind_T{i} = Combined(:,1);
         if std(diff(Combined(:,2)))/mean(diff(Combined(:,2))) > 0.03
            disp(['Combined sub tomogram ', num2str(i), ' has too large deviations in the angles. The sub tomograms are probably not equally spaced']) 
         end  
       %  std(diff(Combined(:,2)))/mean(diff(Combined(:,2)))
         sub_thetaT{i} = Combined(:,2);    
 end
end