% Detect an interruption of the scan in the Log file
% From spec compile  post_scan.mac
% = 1, if the scan is ready for processing
% = 0, if the scan was interrupted and is not yet repeated, or if the scan
%      is being repeated now

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

function out = check_scan_interruption(specLogFile,scanno)
    
    specLogFile = find_specLogFile(specLogFile);
    
    cmd = sprintf('cat -n %s | grep "S %d " | cut -f 1 -d ''#''', specLogFile,scanno);
    [~,sysout] = system(cmd);
    arroutbegin = regexp(sysout,'[:\n]','split');
    cmd = sprintf('cat -n %s | grep "X %d" | cut -f 1 -d ''#''', specLogFile,scanno);
    [~,sysout] = system(cmd);
    cmd = sprintf('cat -n %s | grep "X %d" | cut -f 2 -d ''#''', specLogFile,scanno);
    [~,sysout_2] = system(cmd);
    arroutend = regexp(sysout,'[:\n]','split');
    arroutend_text = regexp(sysout_2,'[:\n]','split');
    numbegin = [];
    numend = [];
    for ii = 1:numel(arroutbegin)-1
        numbegin(ii) = str2num(arroutbegin{ii});
    end
    count = 1;
    for ii = 1:numel(arroutend)-1
        if strcmp(arroutend_text{ii}(3:end),num2str(scanno))
            numend(count) = str2num(arroutend{ii});
            count = count+1;
        end
    end
    if numel(numbegin)~=numel(numend)
        % The experiment is being repeated now
        out = 0;
        return
    else
        cmd = sprintf('head -n %d %s | tail -n %d|grep chk',numend(end) ,specLogFile,numend(end)-numbegin(end));
        [~,sysout] = system(cmd);
        errorout = regexp(sysout,'[:\n]','split');
        if isempty(errorout{1})
            out = 1;
        else
            % There was an error in this scan
            out = 0;
        end
        return
    end
    
    
%     if any(strcmp(arrout,sprintf('#X %d',scanno)))
%         out = true;
%     else
%         out = false;
%     end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function specLogFile = find_specLogFile(specLogFile)
  while (exist(specLogFile,'file') ~= 2)
% if the variable specLogFile is not a complete path to a file
%   try to guess where a spec data file can be found, by
%   - look for directories called 'spec' or 'dat-files'
%   - look for files called '*.dat'
%   - take the newest one
    compare_str = specLogFile;
    fname = dir(specLogFile);
    if (exist(specLogFile,'dir'))
      if (specLogFile(end) ~= '/')
        specLogFile = strcat(specLogFile,'/');
      end
      
      for ii=1:numel(fname)
        if (regexp(fname(ii).name,'.log$'))
          specLogFile = strcat(specLogFile,'*.log');
          fname = [];
          break;
        end
      end
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,'log-files'))
          specLogFile = strcat(specLogFile,fname(ii).name);
          fname = [];
          break;
        end
      end
      for ii=1:numel(fname)
        if (strcmp(fname(ii).name,'specES1'))
          specLogFile = strcat(specLogFile,fname(ii).name);
          break;
        end
        if (strcmp(fname(ii).name,'spec'))
          specLogFile = strcat(specLogFile,fname(ii).name);
          break;
        end
      end
    else
      if (numel(fname)>0)
        [~,ii] = max(cell2mat({fname.datenum}));
        specLogFile = regexprep(specLogFile,'\*\.log$',fname(ii).name);
      else
        error('''%s'' cannot be found.', specLogFile);
        break
      end
    end
    if (strcmp(specLogFile,compare_str))
      break
    end
  end
end