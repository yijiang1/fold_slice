% show recons.m
% Warning: Currently working only for square pixels
% close all   % Recommended if improfile will be used (there is a bug with cursor positioning otherwise)
% Mayor changes, basically rewritten, made on Oct 19, 2015 in order to accomodate waiting for
clear
import utils.*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show recons parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_path='~/Data10/';
addpath '~/Data10/matlab'
addpath ~/Data10/matlab/ptycho/
colorbarphase = [-1 1]*pi; % Give the range, or 'auto'
saveplots = 1;  % Saves JPEGs of reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters to find the file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Option 1 %%%    % Path + name
filefixname = [];   % Provide a full path and filename. Leave empty to use the next option.

    %%% Option 2 %%%    % Range of scan numbers
scans      = fliplr([17000:25000]);    % Specify a range of scan numbers, the code will try to be smart and find the reconstructions. Leave empty to use the next option.
    prefix = '';    % Define a prefix to choose one reconstruction if there are many in the folder. Leave empty to just grab the first one.        
    suffix = '_recons';    % Alternatively you can define a suffix.

    %%% Option 3 %%%    Specify an OMNY/flOMNI dat file path. The code will
                        % keep looking in this folder, it finds the file and moves it to the
                        % second folder. If the name of second folder is specified, if left
                        % empty it will not move it
queue_path     = ['~/Data10/specES1/recontruct/done/'];
queue_path_out = ['~/Data10/specES1/recontruct/done_shown/'];   % If you leave this empty it will not move the files, but then this is pretty useless, eh?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

continuewithloop = true;
verbose(0); % Change to   = 2   in order to have more output on the 
varargs{1} = 'PhaseColorBarAxis';
varargs{2} = colorbarphase;
if saveplots
    varargs{3} = 'ImageSaveFolder';
    varargs{4} = fullfile(base_path,'analysis/online/ptycho/show_recons/');
end
scans_to_do = scans;

scans_plotted = 0;

while(~isempty(scans_to_do))
%     while(continuewithloop)
        if ~isempty(filefixname)    % Just use the fixed name
            file = filefixname;
            continuewithloop = false;
            plotting.ptycho_show_recons(file,varargs);

        else
            if isempty(scans)   % Use OMNY reconstruct dat file
                %%% Find file with task to plot
                lookforafile = true;
                verbose(1,['queue_path is active, looking for files in the queue in ' queue_path]);
                while(lookforafile)
                    files_recons = dir([queue_path 'scan*']);
                    if ~isempty(files_recons)
                        file_dat = fullfile(queue_path,files_recons(1).name);
                        verbose(1,['Found file in queue ' file_dat]);
                        lookforafile = false;
                    else
                        verbose(1,sprintf('Did not find files in the queue: %s, pausing 10 sec',queue_path));
                        pause(10)
                    end
                end
                %%%
                p_out = parse_queue_file(file_dat);
                scanstoplot = p_out.scan_number;
                % Now move the file
                if ~exist(queue_path_out,'dir')
                    warning(sprintf('Creating folder %s',queue_path_out))
                    mkdir(queue_path_out);
                end
                verbose(1,sprintf('Moving %s to %s',file_dat,queue_path_out))
                movefile(file_dat,queue_path_out)
            else    % Use the scan numbers provided
                continuewithloop = false;
                scanstoplot = scans;
            end
            
            % Ok, now it knows which scans to plot
            num_orig_args = numel(varargs);
            for scannum = scans_to_do%scanstoplot
                file = find_ptycho_filename(base_path,scannum,prefix,suffix); % Compile name
%                 waitingforrecons = false;
%                 while(waitingforrecons)
                    if iscell(file)
                        warning('More than one file found, using the 1st one. Consider specific filename.');
                        file = file{1};
                    end
                    if exist(file,'file')
                        verbose(1,sprintf('Found %s',file));
                        verbose(1,'Waiting 3 sec to make sure the file is written')
                        pause(3)
                        varargs{num_orig_args+1} = 'ScanNumber';
                        varargs{num_orig_args+2} = scannum;
                        try
                            JavaObj = java.lang.Runtime.getRuntime;
                            fprintf('Free memory: %d\n', JavaObj.freeMemory/1e6)
                            plotting.ptycho_show_recons(file,varargs);
                            scans_plotted = scans_plotted+1;
                            save('debug11.mat','scans_plotted','scans_to_do');
                        catch
                            fprintf('failed to open file\n')
                        end
                        waitingforrecons = false;
                        scans_to_do(find(scans_to_do==scannum)) = [];
                    else
                        verbose(1,sprintf('Reconstruction %s not found, pausing 10 sec',file));
                        %                     pause(10)
                    end
%                 end
            end
        end
%     end
end

return

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


