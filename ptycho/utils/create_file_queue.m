% CREATE_FILE_QUEUE create directories and .dat files for a file queue
% only needed if OMNY/flOMNI is not available
% ** start          start scan number
% ** stop           stop scan number
% ** step           reconstruction bundle size
% 
% *optional*        given as name/value pair
% ** dirpath        path to the .dat files directory; default: '../reconstruction'
% ** prop           structure with parameters that need to be saved to the .dat file
% ** prop_name      structure name; default: 'p'
% ** spec           path to spec file; needed for check2detpos
% ** check2detpos   check detector positions to avoid repeated scans
% ** det_motor      detector motor; default 'dettrx'
% ** verbose        set verbose level; default 0
%
% EXAMPLE:
%   s.lockfile = false;
%   s.energy = 6.2015;
%   s.check_nextscan_started = 1;
%   
%   create_file_queue(646, 650, 2, 'prop', s)
%
%
%   create_file_queue(287, 2553, 2, 'prop', s, 'spec', '/das/work/p16/p16812/data/pilatus/e16403/', 'check2detpos', true, 'det_motor', 'dettrx')

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

function create_file_queue(start, stop, step, varargin)
import utils.*

par = inputParser;
par.addParameter('dirpath', '../reconstruction', @ischar) 
par.addParameter('prop', [], @isstruct)
par.addParameter('prop_name', 'p', @ischar)
par.addParameter('spec', [], @ischar)
par.addParameter('check2detpos', false, @logical)
par.addParameter('det_motor', 'dettrx', @ischar)
par.addParameter('verbose', 0, @isnumeric)

par.parse(varargin{:})

var = par.Results;

% reduce verbose level
vbl = utils.verbose;
utils.verbose(var.verbose);

% prepare directories
if ~exist(fullfile(var.dirpath), 'dir')
    mkdir(fullfile(var.dirpath))
end
if ~exist(fullfile(var.dirpath, 'in_progress'), 'dir')
    mkdir(fullfile(var.dirpath, 'in_progress'))
    mkdir(fullfile(var.dirpath, 'done'))
    mkdir(fullfile(var.dirpath, 'failed'))
end

% filename pattern
fname_pattern = repmat('scan%05d_',1,step);
fname_pattern = [fname_pattern(1:end-1) '.dat'];
scan_number_pattern = repmat('%d ',1,step);
scan_number_pattern = scan_number_pattern(1:end-1);

% get value for p.detector.check_2_detpos
if var.check2detpos
    
    S = io.spec_read(var.spec, 'ScanNr', start);
    detpos1 = S.(var.det_motor);
    for ii=start+1:stop
        S = io.spec_read(var.spec, 'ScanNr', ii);
        if S.(var.det_motor)~=detpos1
            detpos2=S.(var.det_motor);
            break
        end
    end
    var.prop.detector.check_2_detpos = abs(detpos1-detpos2)/2+detpos2;
end


% if needed, prepare structure
if ~isempty(var.prop)
    prop_fnames = utils.struc2cell(var.prop);
else
    prop_fnames = [];
end

utils.verbose(0, 'Creating file queue.');

ii = start;
if utils.verbose < 2
    utils.progressbar(1, round((stop-start)/step))
end
while ii<=stop
    % check if detector positions are repeated
    if var.check2detpos
        detpos = [];
        kk= 1;
        for jj=ii:ii+step-1
            S = io.spec_read(var.spec, 'ScanNr', jj);
            detpos(kk) = S.(var.det_motor);
            kk = kk+1;
        end
        if all(detpos == detpos(1))
            utils.verbose(2,'Found repeated detector positions.')
            ii = ii+1;
            continue
        end
    end
    
    % write file            
    fname = sprintf(fname_pattern, ii:ii+step-1);
    fid = fopen(fullfile(var.dirpath, fname), 'w');
    fprintf(fid, ['p.scan_number \t ' scan_number_pattern '\n'], ii:ii+step-1);
    if ~isempty(prop_fnames)
        for jj=1:numel(prop_fnames)
            cval = eval(['var.prop.' prop_fnames{jj}]);
            if islogical(cval) || isnumeric(cval)
                %eval(['var.prop.' prop_fnames{jj} '= double(var.prop.' prop_fnames{jj} ')']);
                prop_val = num2str(cval);
            else
                prop_val = cval;
            end
            prop_val = replace(prop_val, '%', '%%');   % avoid interpretation of special characters in strings 
            fprintf(fid, [var.prop_name '.' prop_fnames{jj} '\t ' prop_val '\n']);
        end
    end
    
    fclose(fid);    
    
    ii = ii+step;
    if utils.verbose < 2
        utils.progressbar(ii-start, round((stop-start)/step))
    end
end

utils.verbose(0, '\nDone.');

% revert changes to verbose level
utils.verbose(vbl);

end
