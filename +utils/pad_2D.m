% pad_2D pads a 2D array (in both direction). It's based on Matlab's padarray function
% Inputs: 
%   **img                input image
%   **outsize            size of final image
% *optional:*
%   **padval               value to fill padded regions 
% returns: 
%   ++imout              cropped image 
% Written by YJ

function [ imout ] = pad_2D( img, outsize, padval)

Nin = size(img);
Nout = outsize(1:2);

if Nout(1)<Nin(1) || Nout(2)<Nin(2)
    disp(size(img))
    error('Output size is smaller than input image!')
end

if nargin < 3
    padval = 0;
end

pad_pre = [0,0];
pad_post = [0,0];
%calculate how much to pad 
if mod(Nin(1),2)==0 %if input image size is even
    pad_post(1) = ceil((Nout(1)-Nin(1))/2);
    pad_pre(1) = floor((Nout(1)-Nin(1))/2);
else %odd
    pad_post(1) = floor((Nout(1)-Nin(1))/2);
    pad_pre(1) = ceil((Nout(1)-Nin(1))/2);
end

if mod(Nin(2),2)==0 %if input image size is even
    pad_post(2) = ceil((Nout(2)-Nin(2))/2);
    pad_pre(2) = floor((Nout(2)-Nin(2))/2);
else %odd
    pad_post(2) = floor((Nout(2)-Nin(2))/2);
    pad_pre(2) = ceil((Nout(2)-Nin(2))/2);
end

%imout = padarray(img, [(Nout(1)-Nin(1))/2, (Nout(2)-Nin(2))/2],padval);


imout = padarray(img, pad_pre,padval,'pre');
imout = padarray(imout, pad_post,padval,'post');

if ~isreal(img)
    imout = complex(imout);
end

end