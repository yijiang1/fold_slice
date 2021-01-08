%GET_FILENAMES_CSAXS compile filenames of raw data files
% receives 

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

function [p] = get_filenames_cSAXS(p)
import utils.*

% get detector paramters
det = p.detectors(p.scanID).params;
read_path = p.raw_data_path_full{p.scanID};
detStorage = p.detectors(p.scanID).detStorage;

if isfield(det, 'filename_pattern')
    if iscell(det.filename_pattern)
        %% if filename patterns exist, use them to restrict the file search
        det.filename_pattern_full = [p.detector.data_prefix];
        fill = [];
        
        for ii=1:size(det.filename_pattern,2)
            fill = [fill det.filename_pattern{ii}.str det.filename_pattern{ii}.del];
            det.filename_pattern_full = [det.filename_pattern_full det.filename_pattern{ii}.str det.filename_pattern{ii}.del];
            del{ii} = det.filename_pattern{ii}.del;
            
            switch det.filename_pattern{ii}.content
                case 'burst'
                    burst = ii;
                case 'scan'
                    scan = ii;
                case 'pos'
                    pos = ii;
            end
            
        end
        
        det.filename_pattern_full = [det.filename_pattern_full det.file_extension];
        
        if exist('burst', 'var')
            det.filename_pattern_burst = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            parts{burst} = '*';
            det.filename_pattern_burst = [det.filename_pattern_burst strjoin(parts, del) det.file_extension];
        end
        
        if exist('pos', 'var')
            det.filename_pattern_pos = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            parts{pos} = '*';
            det.filename_pattern_pos = [det.filename_pattern_pos strjoin(parts, del) det.file_extension];
        end
        
        if exist('scan', 'var')
            det.filename_pattern_scan = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            if exist('burst', 'var')
                parts{burst} = '*';
            end
            if exist('pos', 'var')
                parts{pos} = '*';
            end
            det.filename_pattern_scan = [det.filename_pattern_scan strjoin(parts, del) det.file_extension];
        end
        
        input_vars = {};
        if isfield(det, 'filename_pattern_pos')
            k = 1;
            for jj=1:length(det.filename_pattern)
                switch det.filename_pattern{jj}.content
                    case 'pos'
                        continue
                    case 'burst'
                        input_vars{k} = det.filename_pattern{jj}.start;
                    case 'scan'
                        input_vars{k} = p.scan_number(p.scanID);
                end
                k = k+1;
            end
            filename_pattern_pos = fullfile(read_path, det.filename_pattern_pos);
            filename_pos = sprintf(filename_pattern_pos, input_vars{:});
            [~, pos_files] = find_files(filename_pos);
            numpos = size(pos_files,2);
            
            % apply natural sorting order i.e. sort 1,2,3,10,200 and not 1 10 100 2 20 200
            % important if the file makes are not defined as S%05i but rather S%i
            [~,idx] = natsort({pos_files.name});
            pos_files = pos_files(idx);
            
        end
        
        if isfield(det, 'filename_pattern_burst')
            k = 1;
            for jj=1:length(det.filename_pattern)
                switch det.filename_pattern{jj}.content
                    case 'pos'
                        input_vars{k} = det.filename_pattern{jj}.start;
                    case 'burst'
                        continue
                    case 'scan'
                        input_vars{k} = p.scan_number(p.scanID);
                end
                k = k+1;
                
            end
            filename_pattern_burst = fullfile(read_path, det.filename_pattern_burst);
            filename_burst = sprintf(filename_pattern_burst, input_vars{:});
            [~, burst_files] = find_files(filename_burst);
            numburst = size(burst_files,2);
        else
            numburst = 1;
        end
        
        
        
        if numburst > 1
            % if burst frames exist, we need to make sure that the file order is correct
            file_args = '[';
            for ii=1:length(det.filename_pattern)
                switch ii
                    case burst
                        file_args = [file_args ' det.filename_pattern{ii}.start + burstID-1'];
                    case pos
                        file_args = [file_args ' det.filename_pattern{ii}.start + posID-1'];
                    case scan
                        file_args = [file_args ' p.scan_number(p.scanID)'];
                end
            end
            file_args = [file_args ']'];
            for posID=1:numpos
                for burstID=1:numburst
                    files(burstID+(posID-1)*numburst).name = sprintf(det.filename_pattern_full, eval(file_args));
                end
            end
            datadir = read_path;
            
        else
            % if there are no burst frames, use the pos files
            datadir = read_path;
            files = pos_files;
        end
        
        if numel(files)==0
            error('Could not find any files using the filename pattern %s. \n ', fullfile(read_path, strrep(det.filename_pattern_full, '%', '%%')))
        end
    else
        % use wildcards
        files = find_files(fullfile(read_path, [det.filename_pattern det.file_extension]));
        if numel(files)==0
            error('Could not find any files using the filename pattern %s. \n ', fullfile(read_path, strrep(det.filename_pattern_full, '%', '%%')))
        end
    end
else
    % if no filename pattern was specified, just load everything containing the specified file extension
    [datadir, files] = find_files([read_path '*.' det.file_extension]);
    
    if numel(files)==0
        error('Could not find any files using the filename pattern %s.\n ', fullfile(read_path, ['*.' det.file_extension]))
    end
end


detStorage.files = [];

for ii=1:length(files)
    detStorage.files{ii} = fullfile(datadir, files(ii).name);
end

for ii=1:length(det.image_read_extraargs)
    if strcmpi(det.image_read_extraargs{ii}, 'H5Location')
        detStorage.h5_group{1} = det.image_read_extraargs{ii+1};
        break;
    end
end

end

