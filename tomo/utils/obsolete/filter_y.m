%   tomo_filtered = filter_y(tomo,freq_scale)
%   Receives a tomogram and applies filtering only along the third index of
%   a 3D matrix. 
% Inputs:
%   tomo        Input tomogram matrix
%   freq_scale  Frequency cutoff
% Manuel Guizar Feb 14, 2016


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


function tomo_filtered = filter_y(tomo_filtered,freq_scale)


%%% Filter a tomogram along the last index %%%
% freq_scale = 0.1;
% 
% %%% Create example %%%
% N = 200;
% tomo_filtered = zeros([N N N]);
% for ii = 1:N
%     tomo_filtered(ii,:,:) = phantom(N);
% end
% tomo_filtered = repmat(tomo_filtered,[1 1 N]);
%%%%%%%%%%%%%%%%%%%%%%

%%% Create filter %%%

Nfilt = size(tomo_filtered,3);

filt = zeros([Nfilt 1]);
d = freq_scale;

w = [-Nfilt/2:Nfilt/2-1]+mod(Nfilt/2,2);
w = 2*pi*w/Nfilt;

filt = (1+cos(w/d)) / 2;
filt(abs(w)/d>pi) = 0;
filt = ifftshift(filt);

% figure(1);
% plot(filt);
%%%%%%
%

filt3D = repmat(reshape(filt,[1,1,Nfilt]),[size(tomo_filtered,1) size(tomo_filtered,2) 1]);

tomo_filtered = real(ifft( fft(tomo_filtered,[],3).*filt3D ,[],3));



% %% Image %%%
% figure(1)
% imagesc(squeeze(tomo_filtered(2,:,:)));
% % imagesc(squeeze(filt3D(1,:,:)));
% colorbar

