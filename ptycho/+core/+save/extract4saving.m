%EXTRACT4SAVING extracts datasets and parameters from p and creates HDF5
%structure
% ** p              p structure
% ** append         boolean; true if data will be appended to an existing file
%
% returns:
% ++ s              structure for save2hdf5
%
% see also: io.HDF.save2hdf5
%

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

function [ s ] = extract4saving(p, append)
import utils.update_param
import math.double2int

s = [];



p = rmfield_safe(p, 'probe');
p = rmfield_safe(p, 'positions_temp');
p = rmfield_safe(p, 'scanidxs');
p = rmfield_safe(p, 'share_pos');

p.positions = transpose(p.positions);
p.positions_real = transpose(p.positions_real);
p.positions_orig = transpose(p.positions_orig);

%%%%%%%%%%%%%%%%%%%
%%% measurement %%%
%%%%%%%%%%%%%%%%%%%
%% external link to data file
% should be saved with a relative path
% s.measurement.data = ['ext:' p.prepare_data_path p.prepare_data_filename ':/'];

%% Meta data 
if ~isempty(p.meta)
    s.measurement.meta_all = p.meta;
end
p = rmfield_safe(p, 'meta');


%% Detector settings
% s.measurement.detector
p = rmfield_safe(p, p.detector.name);

%% fmask and fmag
p = rmfield_safe(p, 'fmask');
p = rmfield_safe(p, 'fmag');



%%%%%%%%%%%%%%%%%%%%%%
%%% reconstruction %%%
%%%%%%%%%%%%%%%%%%%%%%

%% engines
em_indx = 0;
for ii=1:length(p.engines)
    tmp = p.engines{ii};

%     % add object_final
%     if ~isempty(tmp.object_final)
%         s.reconstruction.engines{ii}.object_final = tmp.object_final;
%     end
%     % add probes_final
%     if ~isempty(tmp.probes_final)
%         s.reconstruction.engines{ii}.probes_final = tmp.probes_final;
%     end
    if isfield(tmp, 'error_metric_final')
    for jj=1:length(tmp.error_metric_final)
        % add error_metric_final
        if iscell(tmp.error_metric_final)
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).iteration = tmp.error_metric_final{jj}.iteration;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).value = tmp.error_metric_final{jj}.value;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).method = tmp.error_metric_final{jj}.method;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).err_metric = tmp.error_metric_final{jj}.err_metric;
            s.reconstruction.p.engines{ii}.error_metric_final.Attributes.MATLAB_class = 'cell';
        else
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).iteration = tmp.error_metric_final.iteration;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).value = tmp.error_metric_final.value;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).method = tmp.error_metric_final.method;
            s.reconstruction.p.engines{ii}.error_metric_final.(['em_' num2str(jj-1)]).err_metric = tmp.error_metric_final.err_metric;
            s.reconstruction.p.engines{ii}.error_metric_final.Attributes.MATLAB_class = 'cell';
        end
        % add error_metric
%         s.reconstruction.p.error_metric.(['em_' num2str(em_indx)]) = ['int_soft:/reconstruction/p/engines/' fn{ii} '/error_metric_final/em_' num2str(jj-1)];
%         em_indx = em_indx + 1;
    end
    s.reconstruction.p.engines{ii} = update_param(s.reconstruction.p.engines{ii}, double2int(tmp), 'force_update', 0);

    end
    
    tmp = rmfield_safe(tmp, 'object_final');
    tmp = rmfield_safe(tmp, 'probes_final');
    tmp = rmfield_safe(tmp, 'error_metric_final');
    tmp = rmfield_safe(tmp, 'fdb');
    
end

p = rmfield_safe(p, 'engines');
p = rmfield_safe(p, 'error_metric');
% p = rmfield_safe(p, 'err');
% p = rmfield_safe(p, 'rfact');

%% probe (dataset)
if ~append
    for ii=1:p.numprobs
        s.reconstruction.p.probes.(['probe_' num2str(ii-1)]) = permute(squeeze(p.probes(:,:,ii,:)), [2 1 3]);
    end
end

p = rmfield_safe(p, 'probes');




%% object (dataset)
if ~append
    for ii=1:p.numobjs
        s.reconstruction.p.objects.(['object_' num2str(ii-1)]) = permute(p.object{ii}, [2 1 3 4]);
    end
end
p = rmfield_safe(p, 'object');



%% probe mask
if isfield(p, 'probe_mask')
    s.reconstruction.p.probe_mask = p.probe_mask;
    
    p = rmfield_safe(p, 'probe_mask');
end

%% ctr
s.reconstruction.p.ctr = uint32(transpose(p.ctr));

p = rmfield_safe(p, 'ctr');

%% everything else

s.reconstruction.p = update_param(s.reconstruction.p, double2int(p), 'force_update', 0);




end

function p = rmfield_safe(p, val)
if isfield(p, val)
   p = rmfield(p, val);
end    
end
    
