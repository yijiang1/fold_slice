%ALIGNED_FSC_TEMPLATE
% Script to align images and compute FSC
%
% References relevant to this code:
% For using this FSC code with ptychography: J. Vila-Comamala, et al., "Characterization of high-resolution diffractive X-ray optics by ptychographic coherent diffractive imaging," Opt. Express 19, 21333-21344 (2011).
% For subpixel alignment: M. Guizar-Sicairos, et al., "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156 (2008).
% For matching of phase ramp by approximate least squared error:  M. Guizar-Sicairos, et al., "Phase tomography from x-ray coherent diffractive imaging projections," Opt. Express 19, 21345-21357 (2011).
%

addpath ../base/
% clear;
% close all;
params = struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openwithGUI = 1;    % Use GUI for choosing files, otherwise specify parameters below
scan1 = [89];            % Scan number of first image
scan2 = [90];            % Scan number of second image
sample_name = '';    % File prefix
suffix = 'test_1_recons.h5';         % File suffix
analysis_folder = '../../analysis/'; % /mnt/das-gpfs/work/p12345/analysis/  % /sls/X12SA/Data10/e12345/analysis/
filenamewithpath1 = ['image1.tif'];     % Give the full filename and path - Overrides the parameters above; Can be in .mat or any format supported by 'imread'
filenamewithpath2 = ['image2.tif'];     % Give the full filename and path - Overrides the parameters above; Can be in .mat or any format supported by 'imread'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Alignment parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.verbose_level = 3;   % adjust output level
params.plotting = 2;        % (3) show everything, (2) show aligned images + FSC, (1) show FSC, (0) none
params.remove_ramp = 1;     % Try to remove ramp from whole image before initial alignment
params.image_prop = 'phasor'; % = 'complex' or = 'phasor' (phase with unit amplitude) or = 'phase'  (Note: phase should not be used if there is phase wrapping)
params.crop = 'manual'; 
            % '' for using the default half size of the probe
            % 'manual' for using GUI to select region. This will display the range, e.g. {600:800, 600:800}     
            % {600:800, 600:800} for custom vertical and horizontal cropping, respectively
params.flipped_images = 0; % If images are taken with a horizontal flip, e.g. 0 & 180 for tomography
params.GUIguess = 0;       % To click for an initial alignment guess, ignores the values below
params.guessx = [];       % Some initial guess for x alignment
params.guessy = [];
%%%%%%%%%%%%%%%%%%%%%%
%%% FSC parameters %%%
%%%%%%%%%%%%%%%%%%%%%%
params.taper = 20;             % Pixels of image tapering (smoothing at edges) - Increase until the FSC does not change anymore
params.SNRt = 0.5;      % SNRt = 0.2071 for 1/2 bit threshold for resolution of the average of the 2 images
                        % SNRt = 0.5    for 1   bit threshold for resolution of each individual image
params.thickring = 10;  % Thickness of Fourier domain ring for FSC in pixels
params.freq_thr = 0.05;  % (default 0.05) To ignore the crossings before freq_thr for determining resolution 


%%%%%%%%%%%%
%%% misc %%%
%%%%%%%%%%%%
params.prop_obj = false;     % propagation distance at the sample plane; leave empty to use the value from the reconstruction p structure; set to "false" for no propagation
params.apod = [];            % if true, applies an apodization before propagating by params.prop_obj, the apodization border is around the valid reconstruction region; leave empty to use the value from the reconstruction p structure
params.lambda = [];           % wavelength; needed for propagating the object; leave empty to use the value from the reconstruction p structure
params.pixel_size = [];       % pixel size at the object plane; leave empty to use the value from the reconstruction p structure

%%%%%%%%%%%%%%%%%%%%
%%% FP parameter %%%
%%%%%%%%%%%%%%%%%%%%

%%% the following parameters are ignored, unless p.fourier_ptycho==true %%%
params.filter_FFT = true;           % apply a circular mask to the reconstructed spectrum (needs p.plot_maskdim)
params.crop_factor = 0.9;           % crop final image by the given factor
params.crop_asize = [800 800];      % crop object before applying the FFT
params.z_lens = 49.456e-3;          % FZP focal distance 






%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do not modify below %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
caller = dbstack;
if length(caller)==1
    
    addpath('utils')
    scanfolder1 = utils.compile_x12sa_dirname(scan1(1)); % Looks for the file in this folder, I leave a variable so that the folder can be overriden
    scanfolder2 = utils.compile_x12sa_dirname(scan2(1));
    %%% Opening file %%%
    if openwithGUI
        disp('Using GUI open mode')
        if exist([analysis_folder scanfolder1],'dir')
            uipath1 = [analysis_folder scanfolder1];
        else
            uipath1 = [];
        end
        if exist([analysis_folder scanfolder2],'dir')
            uipath2 = [analysis_folder scanfolder2];
        else
            uipath2 = [];
        end
        filetypes = {'*.h5;*.mat','Reconstruction files (*.h5,*.mat)'; '*.*',  'All Files (*.*)'};
        [filename, pathname] = uigetfile(filetypes,'Open first reconstruction', uipath1);
        file{1} = fullfile(pathname,filename);
        
        [filename, pathname] = uigetfile(filetypes,'Open second reconstruction', uipath2);
        file{2} = fullfile(pathname,filename);
    else
        % Checking recons 1 %
        if ~isempty(filenamewithpath1)
            file{1} = filenamewithpath1;
        else
            file{1} = fullfile(analysis_folder,scanfolder1,[sample_name '*' suffix]);
            D = dir(file{1});
            if numel(D) == 0
                error(['I did not find any file: ' file{1}])
            elseif numel(D) > 1
                warning(['I found many files with the mask: ' file{1}]);
                warning(['I selected ' D(1).name]);
            end
            file{1} = fullfile(analysis_folder,scanfolder1,D(1).name);
        end
        
        % Checking recons 2 %
        if ~isempty(filenamewithpath2)
            file{2} = filenamewithpath2;
        else
            file{2} = fullfile(analysis_folder,scanfolder2,[sample_name '*' suffix]);
            D = dir(file{2});
            if numel(D) == 0
                error(['I did not find any file: ' file{2}])
            elseif numel(D) > 1
                warning(['I found many files with the mask: ' file{2}]);
                warning(['I selected ' D(1).name]);
            end
            file{2} = fullfile(analysis_folder,scanfolder2,D(1).name);
        end
    end
    
    % Making a JPEG of FSC
    [~,filename] = fileparts(file{1});
    params.out_fn = sprintf('%sonline/ptycho/%s_FSC.jpg', analysis_folder, filename);
    
    [resolution] = aligned_FSC(file{1}, file{2}, params);
end




% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
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
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
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
