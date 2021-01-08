% function convert(scan, raw_data_path)
% convert .raw files to .hdf5
%  Inputs
%     scan              scan number
%     raw_data_path     path to the raw data (optional)

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


function convert2hdf5(scan, raw_data_path)

    convertor_path = '/sls/X12SA/data/x12saop/EigerPackage/slsDetectorsPackage/bin/hdf5MakerOMNY'; 
    load_dir = utils.compile_x12sa_dirname(scan);
    if exist('raw_data_path','var')&&exist(fullfile(raw_data_path,load_dir),'dir')
        load_dir = fullfile(raw_data_path,load_dir);
    elseif exist(['~/Data10/eiger_4/'],'dir')
        load_dir = ['~/Data10/eiger_4/' load_dir]; 
    elseif exist([raw_data_path,'eiger_4/'])
        load_dir = [raw_data_path,'/eiger_4/' load_dir]; 
    elseif exist([raw_data_path,'/eigeromny/'])
        load_dir = [raw_data_path,'/eigeromny/' load_dir]; 
    end
    
    if ~exist(load_dir, 'dir')
        warning('Raw data path %s not found', load_dir)
        return
    end

    testDir = [load_dir, '/deleteMe'];
    % test for write permissions by creating a folder and then deleting it
    isWritable = mkdir(testDir);
    % check if directory creation was successful
    if isWritable == 1
        rmdir(fullfile(testDir));
    end

    list_h5 = dir([load_dir, '/run_*.h5']); 
        
    
    file_sizes = [list_h5.bytes]; 
    if any(file_sizes < 1e6)  % find files < 1MB
        warning('H5 files in scan %i seem damaged, generate again ... ', scan)
        list_raw = dir([load_dir, '/run_d0_f0000000*.raw']);
        if isempty(list_raw)
           error('RAW data are missing, data cannot be converted')
        else
           delete(sprintf('%s/*.h5',load_dir))
        end
        list_h5 = dir([load_dir, '/run_*.h5']); 
    end
    
    if isempty(list_h5)    
        if ~isWritable
            warning('Conversion failed because folder %s is not writable', load_dir)
            return
        end
        
        list_raw = dir(fullfile(load_dir, 'run_d0_f0000000*.raw')); 

        Nscans = length(list_raw); 
        
        for ii = 1:Nscans
            ind_scans(ii) = str2num(list_raw(ii).name(16:17)); 
        end
        
        for ii = 1:Nscans
%            [~,out] = system([cbf_maker_path ' ' fullfile(load_dir,sprintf('run_d0_f00000000%02i000_%i.raw',ind_scans(ii) ,scan))]);     
            systemcall = [convertor_path ' ' fullfile(list_raw(1).folder,list_raw(1).name)];
            fprintf('%s\n',systemcall);
            [~,out] = system(systemcall);
            out
        end
        
        list_h5 = dir([load_dir, '/*.h5']);
        if isempty(list_h5)
           error(sprintf('After conversion did not find any h5 in %s\n',load_dir))
           return
        end
        if numel(list_h5)>1
           error(sprintf('After conversion I found more than one h5 in %s\n',load_dir))
           return
        end
            
        h5fileinfo = h5info(fullfile(list_h5.folder,list_h5.name));
            
        nframes_converted = h5fileinfo.Groups.Groups(1).Datasets.Dataspace.Size(3);
        out = splitlines(out); 
        nframes_expected = str2num(out{end-2}(14:end)); 
        fprintf('Frames expected: %i, frames converted %i \n', nframes_expected, nframes_converted)
        if nframes_converted == nframes_expected
            fprintf('Scan %i succefully converted to H5\n', scan); 
             delete(sprintf('%s/*.raw',load_dir))
        else
            error('Scan %i WAS NOT CONVERTED to H5\n', scan)    
             delete(sprintf('%s/*.h5',load_dir))
        end
    end
end
