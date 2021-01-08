% [object] = generate_virtual_object(dataset, img_path, Npix, angle)
%   create a virtual thickness map that can be used to generate test object
%   for ptycho 
%   Inputs:
%       dataset - number or string with name 
%       img_path - path to folder with source images 
%       Npix - size of the output ptrojections 
%       angle - rotation of the 3D object 
% Outputs: 
%       object -  normalized thickness map 

%   Publications most relevant to the Difference-Map implementation
%       + P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, 
%       "High-Resolution Scanning X-ray Diffraction Microscopy," Science 321, 379-382 (2008)
%       + P. Thibault, M. Dierolf, O. Bunk, A. Menzel, F. Pfeiffer,
%       "Probe retrieval in ptychographic coherent diffractive imaging,"
%       Ultramicroscopy 109, 338–343 (2009)
%
%   Publications most relevant to the Maximum Likelihood refinement
%       + M. Guizar-Sicairos and J. R. Fienup, "Phase retrieval with transverse
%       translation diversity: a nonlinear optimization approach," Opt. Express 16, 7264-7278 (2008)
%       + P. Thibault and M. Guizar-Sicairos, "Maximum-likelihood refinement for
%       coherent diffractive imaging," New J. Phys. 14, 063004 (2012).

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [object] = generate_virtual_object(p, dataset, img_path, Npix, angle)
    %% preload some artificial object 
    switch dataset
        case 1
            %% MONA LISA
            object = imread(fullfile(img_path, 'ML512.jpg'));
            object = double(object(:,:,1));
            object = double(object)./double(max(object(:)));


        case 2 
            %% autumn
            object = imread('autumn.tif');
            object = double(object(:,:,1));
            object = double(object)./double(max(object(:)));



        case 3 
            %% camera man
            object = imread('autumn.tif');
            object = double(object(:,:,1));
            object = double(object)./double(max(object(:)));
            

        case 4 
            %% test pattern
            object = imread('testpat1.png');
            object = double(object(:,:,1));
            object = double(object)./double(max(object(:)));
            
        case 5
            %% USAF target 
            object = 255-imread(fullfile(img_path, 'USAF-1951.png'));
            object = single(object(:,:,1) < 128);
        case 6
            %% MANDRIL 
            object = imread(fullfile(img_path, 'mandrill.png'));
            object = double(object)./double(max(object(:)));
            
        case 7
            %% CHIP
            object = imread(fullfile(img_path, 'chip_phantom.png'));
            object = single(object(:,:,1)) / 255; 
        
        case 8
            %% Snellen chart 
            object = imread(fullfile(img_path, 'Snellen_chart.png'));
            object = single(object(1:2:end,1:2:end,1) > 50); 
            object = utils.crop_pad(object, size(object)+10,1); 
            
       case '3D_phantom'
            %% create a phantom that looks like a porous material and henerate projections from it, 
            import utils.*
            import math.*
            

%             Nx = ceil(Npix(1)/128)*128;
%             Nz = ceil(Npix(2)/128)*128;
%             Ny = Nx;
            Nx = 512; 
            Ny = 512; 
            Nz = 512; 
            
            if check_option(p, 'use_gpu')
                dtype = gpuArray.ones(1,'single'); 
                verbose(0,'Processing phantom on GPU')
            else
                dtype = ones(1,'single'); 
            end
            
            phantom_filename = fullfile(img_path, sprintf('glass_data_%ipx.mat', Nx));
            try
                d = load(phantom_filename); 
                volData = single(d.volData)/255;  % convert from uint8 to singles 
                [Nx, Ny, Nlayers] = size(volData);
            catch
                verbose(0, 'Cached phantom not availible in %s', phantom_filename)
                disp('Creating phantom')
                % create a phantom that looks like a porous material
                rng default 
                volData = randn([Nx, Ny, Nz],'like', dtype);
                utils.progressbar(1,6);
                volData = imgaussfilt3_fft(volData, 2); 
                utils.progressbar(2, 6);
                volData = volData / max(abs(volData(:)));
                utils.progressbar(3,6);
                volData = single(volData <= 0);
                utils.progressbar(5, 6);

%                 volData = volData*0.2 + 0.8; % add some background 
                % store as uint8 to save space 
                volData = uint8(volData/max(volData(:)) * 255); 
                
                if isa(volData, 'gpuArray')
                    volData = gather(volData);
                end
                try; save('-v6', phantom_filename, 'volData'); end
                volData = single(volData)/255;
                utils.progressbar(6, 6);
            end
            
            if ~isempty(angle)
                verbose(0,'Phantom volume is rotated by %6.3g degrees', angle)
                if check_option(p, 'use_gpu')
                    volData = gpuArray(volData); 
                end
                % apply circular mask 
                xgrid = -Nx/2+1 : Nx/2;
                ygrid = -Ny/2+1 : Ny/2;
                [X,Y]  = meshgrid(ygrid, xgrid);
                air_gap = 200; % pixels 
                volData = volData .* imgaussfilt(single(X.^2+Y.^2 < (Nx/2 - air_gap/2)^2), 3);
                
                % rotate the volume using FFT  (volData is real valued at this point)
                volData = utils.imrotate_ax_fft(volData, angle, 3);
                volData = max(0,min(volData,1));
            end
            
            object = 1-permute(volData, [3,1,2]); 

        otherwise
            error('Sorry, sample %d is missing. Feel free to add it!', dataset)
                        
            
    end
    
    if size(object,3) == 1
        % in case of a simple 2D object replicate the pattern to avoid void
        % space around 
        object = repmat(object, 5,5); 
    end

    %% adjust the object to the requested size and surround by empty space 
    object = utils.crop_pad(object, Npix, mean(mean(mean(object(:,[1,end],:))))); 

    object = reshape(object,size(object,1), size(object,2),1,size(object,3));
end
