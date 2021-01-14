% GET_RADIAL_INTEGRATION_MASK create 2D integer array that serves as a
% template for get_integration_matrix for radial integration 
% 
%
% [radial_integration_mask, radial_mask, sector_mask] = get_radial_integration_mask(Np, Nrad, Nsec, center_pos)
%
%  Inputs: 
%   **Np   - size of the integrated frames 
%   **Nrad - number of radial rings 
%   **Nsec - number of angular sectors 
%   **center_pos - position of center in pixels, e.g. Np/2 for well centered dataset 
%  Outputs: 
%   ++radial_integration_mask -  2D integer array integration mask 
%   ++radial_mask   -  2D integer array integration mask of only radial rings 
%   ++sector_mask   -  2D integer array integration mask of only angular sectors 
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



function [radial_integration_mask, radial_mask, sector_mask] = get_radial_integration_mask(Np, Nrad, Nsec, center_pos)
    % generate 2D integration masks - radial + sectors 
    
    offset = center_pos - Np/2; 
    xgrid = (-floor(Np(2)/2)+1:floor(Np(2)/2))+offset(2);
    ygrid = (-floor(Np(1)/2)+1:floor(Np(1)/2))+offset(1);
    
    [X,Y] = meshgrid(xgrid, ygrid);
    R = sqrt(X.^2 + Y.^2);
    Phi = atan2(X,Y);

    % calculate array corresponding to rings 
    r_all = linspace(0, max(Np(1:2))/2, Nrad+1);
    radial_mask = zeros(Np(1:2));
    for i = 1:Nrad
        radial_mask(R >= r_all(i) & R < r_all(i+1)) = i;
    end

    % calculate array corresponding to sectors 
    sec_all = linspace(-pi,pi,Nsec+1);
    sector_mask = zeros(Np(1:2));
    for i = 1:Nsec
        sector_mask(Phi >= sec_all(i) & Phi < sec_all(i+1)) = i;
    end
    
    
    % generate joined integration mask 
    radial_integration_mask = double(radial_mask + (sector_mask-1) .* Nrad);
    radial_integration_mask(radial_mask ==0 | sector_mask == 0) = 0;

end