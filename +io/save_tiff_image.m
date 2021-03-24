function [] = save_tiff_image(I, path)
% Wrapper function for saving 16-bit .tiff image
% Written by YJ
% ** I: 2D image
% ** path: file path

tagstruct.ImageLength     = size(I,1);
tagstruct.ImageWidth      = size(I,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t = Tiff(path,'w');
t.setTag(tagstruct)
t.write(uint16(mat2gray(I)*2^16));
t.close();

end

