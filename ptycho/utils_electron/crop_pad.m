% CROP_PAD adjusts the size by zero padding or cropping
% Inputs: 
%   **img                input image
%   **outsize            size of final image
% *optional:*
%   **fill               value to fill padded regions 
% returns: 
%   ++imout              cropped image 


function [ imout ] = crop_pad( img, outsize, fill)

if nargin < 1
    fprintf('CROP_PAD: adjusts the size by zero padding or cropping\n');
    fprintf('crop_pad(img, outsize)\n');
    return
end

Nin = size(img);

if isempty(outsize) || all(outsize(1:2) == Nin(1:2))
    imout = img;   % if outsize == [], return the same image without changes 
    return
end

Nout = outsize(1:2);

if nargin < 3
    fill = 0;
end



center = floor(Nin(1:2)/2)+1;

imout = zeros([Nout,Nin(3:end)],'like',img);

if fill ~= 0 
    imout = imout + fill; 
end

centerout = floor(Nout/2)+1;

cenout_cen = centerout - center;
imout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2)),:,:) ...
    = img(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)),:,:);

if ~isreal(img)
    imout = complex(imout);
end

end

