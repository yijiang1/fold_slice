function [ ] = save_tiff_image_stack(input, par)
% Save 3D arry as a (8-bit) tiff stack 
% Written by YJ
% ** input: 3D real array
% ** par: parameters

input = double(input);

if ~isfield(par,'axis'); par.axis = 3; end
if ~isfield(par,'range'); par.range = [0,0]; end
if ~isfield(par,'norm_single_slice'); par.norm_single_slice = false; end
if ~isfield(par,'name'); par.name = 'image_stack.tiff'; end

if any(par.range)
    range = par.range;
else
    range = [min(input(:)),max(input(:))];
end

switch par.axis
    case 3
        if par.norm_single_slice
            imwrite(mat2gray(input(:,:,1)),par.name,'tiff')
        else
            imwrite(mat2gray(input(:,:,1),range),par.name,'tiff')
        end
        for i=2:size(input,3)
            if par.norm_single_slice
                imwrite(mat2gray(input(:,:,i)),par.name,'tiff','WriteMode','append')
            else
                imwrite(mat2gray(input(:,:,i),range),par.name,'tiff','WriteMode','append')
            end
        end
    case 2
        if par.norm_single_slice
            imwrite(mat2gray(squeeze(input(:,1,:))),par.name,'tiff')
        else
            imwrite(mat2gray(squeeze(input(:,1,:)),range),par.name,'tiff')
        end
        for i=2:size(input,2)
            if par.norm_single_slice
                imwrite(mat2gray(squeeze(input(:,i,:))),par.name,'tiff','WriteMode','append')
            else
                imwrite(mat2gray(squeeze(input(:,i,:)),range),par.name,'tiff','WriteMode','append')
            end
        end
    case 1
        if par.norm_single_slice
            imwrite(mat2gray(squeeze(input(1,:,:))),par.name,'tiff')
        else
            imwrite(mat2gray(squeeze(input(1,:,:)),range),par.name,'tiff')
        end
        for i=2:size(input,1)
            if par.norm_single_slice
                imwrite(mat2gray(squeeze(input(i,:,:))),par.name,'tiff','WriteMode','append')
            else
                imwrite(mat2gray(squeeze(input(i,:,:)),range),par.name,'tiff','WriteMode','append')
            end
        end
    otherwise
        disp('Wrong axis! par.axis should be 1, 2, or 3.')
end

end

