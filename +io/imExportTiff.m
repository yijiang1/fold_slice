function [ ] = imExportTiff( matrix_input,name,section,varargin)
%Wirte a image stack as a 'Tiff' file
%imExportTiff( matrix,name,section )
%   matrix: 3D matrix(y,x,z)
%   name: output file name
%   section: choose corss sections
%           'XY' matrix(:,:,i)
%           'YZ' matrix(:,i,:)
%           'ZX' matrix(i,:,:)
%       

matrix = matrix_input;
matrix = matrix-min(min(min(matrix)));
s = size(matrix);
m = double(max(max(max(matrix))));

if isempty(varargin)  %gray 
    switch section
        case 'XY'
            imwrite(double(matrix(:,:,1))./m, name,'tiff')
            for i=2:s(3)
                imwrite(double(matrix(:,:,i))./m, name,'tiff', 'WriteMode','append')
            end
        case 'ZY'
            imwrite(double(squeeze(matrix(:,1,:))./m), name,'tiff')
            for i=2:s(2)
                imwrite(double(squeeze(matrix(:,i,:))./m), name,'tiff', 'WriteMode','append')
            end
        case 'ZYs'
            imwrite(mat2gray(double(squeeze(matrix(:,1,:)))), name,'tiff')
            for i=2:s(2)
                imwrite(mat2gray(double(squeeze(matrix(:,i,:)))), name,'tiff', 'WriteMode','append')
            end
        case 'ZX'
            imwrite(double(squeeze(matrix(1,:,:)))./m, name,'tiff')
            for i=2:s(1)
                imwrite(double(squeeze(matrix(i,:,:)))./m, name,'tiff', 'WriteMode','append')
            end
        case 'XYs'
            imwrite(mat2gray(matrix(:,:,1)), name,'tiff')
            for i=2:s(3)
                imwrite(mat2gray(matrix(:,:,i)), name,'tiff', 'WriteMode','append')
            end
        case 'ZXs'
            imwrite(squeeze(matrix(1,:,:))./m, name,'tiff')
            for i=2:s(1)
                imwrite(mat2gray(squeeze(matrix(i,:,:))), name,'tiff', 'WriteMode','append')
            end
    end
else
    a = varargin{1}; %color
    switch section
        case 'XY'
            imwrite(double(matrix(:,:,1))*a*64/m,jet, name,'tiff')
            for i=2:s(3)
                imwrite(double(matrix(:,:,i))*a*64/m,jet, name,'tiff', 'WriteMode','append')
            end
        case 'ZY'
            imwrite(squeeze(matrix(:,1,:))*a*64/m,jet, name,'tiff')
            for i=2:s(2)
            imwrite(squeeze(matrix(:,i,:))*a*64/m,jet, name,'tiff', 'WriteMode','append')
            end

        case 'ZX'
            imwrite(squeeze(matrix(1,:,:))*a*64/m,jet, name,'tiff')
            for i=2:s(1)
            imwrite(squeeze(matrix(i,:,:))*a*64/m,jet, name,'tiff', 'WriteMode','append')
            end
        case 'XYs'
            imwrite(mat2gray(matrix(:,:,1))*a*64,parula, name,'tiff')
            for i=2:s(3)
                imwrite(mat2gray(matrix(:,:,i))*a*64,parula, name,'tiff', 'WriteMode','append')
            end    
    end
end

end

