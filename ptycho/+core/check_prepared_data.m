%CHECK_PREPARED_DATA compares prepared data file with current
%reconstruction parameters. If force_update == true, data preparation will
%be forced (see core.ptycho_prepare_scans)
% ** p                  p structure
% 
% returns:
% ++ force_update       true if data preparation has be done enforced
% 
% see also: core.ptycho_prepare_scans

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

function [ force_update ] = check_prepared_data( p )


force_update = false;

file = fullfile(p.prepare_data_path, p.prepare_data_filename);

if ~exist(file, 'file')
    error('Prepared data file %s does not exist!', file);
end

h = h5info(file);

if ~isempty(find(contains({h.Attributes.Name}, 'format'),1)) && h.Attributes(find(contains({h.Attributes.Name}, 'format'),1)).Value==2
    h5format = 'LibDetXR';
else
    h5format = 'Matlab';
end

switch h5format
    case 'LibDetXR'
        
        % check for sharing, asize and numpts
        idx = contains({h.Groups.Name}, '/measurement');
        fn = length(h.Groups(idx).Groups);
        probe_ID = zeros(1,fn);
        object_ID = zeros(1,fn);
        data_size = zeros(fn,3);
        for ii=1:fn
            group = h.Groups(idx).Groups; 
            idx_sub = ismember({group.Name}, ['/measurement/n' num2str(ii-1)]);
            attrib = group(idx_sub).Attributes;
            idx_pr = ismember({attrib.Name}, 'probe');
            probe_ID(ii) = attrib(idx_pr).Value + 1;
            idx_ob = ismember({attrib.Name}, 'object');
            object_ID(ii) = attrib(idx_ob).Value + 1;
            data_size(ii,:) = group(idx_sub).Datasets(find(contains({group(idx_sub).Datasets.Name}, 'data'),1)).Dataspace.Size;
            
        end
        
        if any(data_size(1,[2,1])~=p.asize/2^p.detector.binning)
            force_update = true;
        end
        
        if  (size(data_size,1) ~= length(p.numpts)  || any(data_size(:,3)'~=p.numpts)) && ~p.fourier_ptycho
            force_update = true;
        end
        
        if length(unique(probe_ID))~=length(unique(p.share_probe_ID)) || any(probe_ID~=p.share_probe_ID)
            force_update = true;
        end
        
        if length(unique(object_ID))~=length(unique(p.share_object_ID)) || any(object_ID~=p.share_object_ID)
            force_update = true;
        end        
        
        
    case 'Matlab'
        
        idx = contains({h.Groups.Name}, '/measurements');
        fn = length(h.Groups(idx).Groups);
        probe_ID = zeros(1,fn);
        object_ID = zeros(1,fn);
        scan_ID = zeros(1,fn);
        for ii=1:fn
            idx_pr = contains({h.Groups(idx).Groups(ii).Attributes.Name}, 'probe');
            idx_det = contains({h.Groups(idx).Groups(ii).Attributes.Name}, 'detector');
            idx_ob = contains({h.Groups(idx).Groups(ii).Attributes.Name}, 'object');
            scan_ID(ii) = h.Groups(idx).Groups(ii).Attributes(idx_det).Value + 1;
            probe_ID(ii) = h.Groups(idx).Groups(ii).Attributes(idx_pr).Value + 1;
            object_ID(ii) = h.Groups(idx).Groups(ii).Attributes(idx_ob).Value + 1;
            
        end  
        numscans = unique(scan_ID);
        numpts = zeros(1,numel(numscans));
        for ii=1:numel(numscans)
            numpts(ii) = sum(scan_ID==numscans(ii));
            if ~all(object_ID(scan_ID==numscans(ii)) == p.share_object_ID(ii)) || ~all(probe_ID(scan_ID==numscans(ii)) == p.share_probe_ID(ii))
                force_update = true;
            end
        end
        
        
        if (numel(numpts) ~= numel(p.numpts)) || any(numpts~=p.numpts)
            force_update = true;
        end
            
        % asize
        asize = hdf5_load(file, '/probes');
        asize = asize(1:2);
        if any(asize ~= p.asize)
            force_update = true;
        end
        
    otherwise
        error('Unknown prepared data format')
end
        


end

