% GET_INTEGRATION_MATRIX Generate sparse integration matrix that sums up
% values according to the provided integration mask 
% 
%
%  int_matrix = get_integration_matrix(mask)
%
%  Inputs: 
%   **mask - 2D integer array, 0 = ignored regions, 1:max(mask) are different sectors that will be summed separatelly 
%  Outputs: 
%   ++int_matrix - 2D sparse matrix 
% 
%
%%%%%%%%%%%%%%%%%%%%% HOW TO USE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% % create some data 
% img = single(imread('cameraman.tif'));
% img = repmat(img, 1,1,10);  % just add there 3rd dimension 
% Np = size(img);
%
% %% define parameters of the integration matrix  
% Nrad = 20;
% Nsec = 30;
% center_pos = Np/2-30;
% 
% 
% %% generate radial and sector masks, needs to be modified if center != Np/2
% [mask, radial_mask, sector_mask] = get_radial_integration_mask(Np, Nrad, Nsec, center_pos);
% 
% %% check the generated sector mask 
% figure(1)
% imagesc(mask); axis off image 
% title('Radial & Angular sectors')
% colormap(hsv)
% drawnow 
% 
% % generate the integration 2D sparse matrix 
% T = get_integration_matrix(mask);
% 
% 
% %% perform sparse matrix based integration
% tic
% img_sum = single(reshape((T*reshape(double(img), prod(Np(1:2)), [])), Nrad,Nsec, []));
% toc
%
% %% perform matlab based integration for comparison 
% tic
% img_sum_0 = zeros(Nrad,Nsec, size(img,3)); 
% for nz = 1:size(img,3)
%     im = img(:,:,nz);
%     for i = 1:Nrad
%         m = radial_mask == i;
%         for j = 1:Nsec
%             img_sum_0(i,j,nz) = sum(im( m & sector_mask == j ));
%         end
%     end
% end
% toc
% 
% 
% % show the first frame to check that the methods are identical 
% figure
% subplot(1,2,1)
% imagesc(img_sum_0(:,:,1)); axis image 
% title('Matlab')
% subplot(1,2,2)
% imagesc(img_sum(:,:,1)); axis image 
% title('Sparse matrix')




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
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%  
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.



function int_matrix = get_integration_matrix(mask)
    N = max(mask(:));
    Np = size(mask);
    int_matrix = zeros(prod(Np),2);
    ind_start = 1; 
    for id = 1 : max(mask(:))
        [i,j] = find(mask == id); 
        ind_end = ind_start + length(i)-1; 
        int_matrix(ind_start:ind_end,:) = [id*ones(length(i),1),i+(j-1)*Np(1)];
        ind_start = ind_end + 1; 
    end
    % convert to sparse matrix 
    int_matrix = sparse(int_matrix(1:ind_end,1),int_matrix(1:ind_end,2),ones(ind_end,1), N, prod(Np)); 

end

