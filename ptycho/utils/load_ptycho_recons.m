%LOAD_PTYCHO_RECONS Load data from cxs/h5 or mat file and return it as 
% structure, single dataset or directly into the workspace.
% An additional argument can be passed to select subsections of the data.
% Loading single datasets is only supported for at least 2 output
% arguments.
%
%   file...     path to cxs/h5 or mat file
%
%   *optional*
%   section...  'full', 'probe', 'object', 'recon' or 'p' to select 
%               subsections of the data; default: 'full'
%
%   EXAMPLES:
%     %% recommended usage %%
%       % load into a structure
%       S = load_ptycho_recons('./recon.h5');
%
%       % load a subset
%       S = load_ptycho_recons('./recon.h5', 'probe');
%
%       % load into single datasets
%       [object, probe, p] = load_ptycho_recons('./recon.h5');
%
%     %% not recommended, only works in 'base' workspace %%
%       % load directly into workspace
%       load_ptycho_recons('./recon.h5');
%
%
%   full = object, probe (current scan) and p
%   recon = object and probe (current scan)
%   probe = probe (current scan)
%   object = object (current scan)
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

function varargout  = load_ptycho_recons( filename_with_path, varargin )

import io.HDF.hdf5_load

varargout = {};

if ~ischar(filename_with_path)
   error('First argument has to be string') 
end

filename_with_path = utils.abspath(filename_with_path); 

if ~exist(filename_with_path, 'file')
    error('Could not find reconstruction file %s', filename_with_path)
end

if nargin > 1
    switch varargin{1}
        case {'pr'; 'probe'; 'probes'}
            section = 'probe';
        case {'ob'; 'obj'; 'objects'}
            section = 'object';
        otherwise
            section = varargin{1};
    end
else
    section = 'full';
end

if ~nargout
    output = 0;
elseif nargout >=2
    output = 2;
else
    output = 1;
end

    function assign_struct(val, val_name)
        switch output
            case 1
                varargout{1}.(val_name) = val;
            case 2
                varargout{end+1} = val;
            otherwise 
                assignin('base', val_name, val);
        end
    end

    function assign_val(struc)
        switch output
            case 1
                varargout{1} = struc;
                
            case 2
                if isfield(struc, 'object')
                    varargout{end+1} = struc.object;
                end
                if isfield(struc, 'probe')
                    varargout{end+1} = struc.probe;
                end
                if isfield(struc, 'p')
                    varargout{end+1} = struc.p;
                end
                    
            otherwise
                fn = fieldnames(struc);
                for ii=1:length(fn)
                    assignin('base', fn{ii}, struc.(fn{ii}))
                end
        end
    end


% check if it is a .mat file or a .cxs file
[~, ~, ext] = fileparts(filename_with_path);
switch ext
    case '.mat'
        switch section
            case 'recon'
                S = load(filename_with_path, 'object', 'probe');
                assign_val(S);

            case 'full'
                S = load(filename_with_path);
                assign_val(S);

            case 'object'
                S = load(filename_with_path, 'object');
                assign_val(S);

            case 'probe'
                S = load(filename_with_path, 'probe');
                size(S)
                assign_val(S);

            case 'p'
                S = load(filename_with_path, 'p');
                assign_val(S);

            otherwise
                error('Unknown data section %s', section);
        end
        
    case {'.cxs','.h5'}
        if io.HDF.hdf5_dset_exists(filename_with_path, 'object', '/reconstruction', true)
            h5_path = '/reconstruction';
        else
            h5_path = '';
        end
        
        % reconstruction
        switch section
            case 'recon'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
            case 'full'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
                % load p
                p = convert2p(hdf5_load(filename_with_path, '/reconstruction/p', '-c'));
                if io.HDF.hdf5_dset_exists(filename_with_path, 'meta_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/meta_all', '-c');
                elseif io.HDF.hdf5_dset_exists(filename_with_path, 'spec_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/spec_all', '-c');
                end
                assign_struct(p, 'p');                
            case 'object'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
            case 'probe'
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
            case 'p'
                % load p
                p = convert2p(hdf5_load(filename_with_path, '/reconstruction/p', '-c'));
                if io.HDF.hdf5_dset_exists(filename_with_path, 'meta_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/meta_all', '-c');
                elseif io.HDF.hdf5_dset_exists(filename_with_path, 'spec_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/spec_all', '-c');
                end
                assign_struct(p, 'p');
            otherwise
                error('Unknown data section %s', section);
        end
                
        
        
    otherwise
        error('Unknown ptycho datatype %s.', ext)
end



end

function tmp = load_data_cell(h)

    fn = fieldnames(h);
    num_end = str2double(subsref(strsplit(fn{1}, '_'), struct('type', '{}', 'subs',{{length(strsplit(fn{1},'_'))}})));
    if length(fn)==2 && (strcmpi(fn{1}, 'i') || strcmpi(fn{1}, 'r'))
        tmp = permute(h.r + 1i*h.i, [2,1,3,4]);
    elseif isnumeric(num_end) && ~isnan(num_end)
        for ii=1:length(fn)
            if isstruct(h.(fn{ii}))
                tmp{ii} = load_data_cell(h.(fn{ii}));
            else
                if isnumeric(h.(fn{ii}))
                    tmp{ii} = double(h.(fn{ii}));
                else
                    tmp{ii} = h.(fn{ii});
                end
            end
        end
%         tmp = h;
    else
        for ii=1:length(fn)
            if isstruct(h.(fn{ii}))
                tmp.(fn{ii}) = load_data_cell(h.(fn{ii}));
            else
                if isnumeric(h.(fn{ii}))
                    tmp.(fn{ii}) = double(h.(fn{ii}));
                else
                    tmp.(fn{ii}) = h.(fn{ii});
                end
            end
        end
    end

end

function tmp = convert2p(h)
    
    fn = fieldnames(h);
    for ii=1:length(fn)
        if isstruct(h.(fn{ii}))
            h.(fn{ii}) = load_data_cell(h.(fn{ii}));
        elseif isnumeric(h.(fn{ii}))
            h.(fn{ii}) = double(h.(fn{ii}));                
        else
            continue;
        end
    end
    tmp = h;
    
    % object
    for ii=1:length(h.objects)
        tmp.object{ii} = permute(h.objects{ii}, [2 1 3 4]);
    end
    tmp = rmfield(tmp, 'objects');
    
    % probes
    pr = tmp.probes;
    tmp.probes = [];
    for ii=1:length(pr)
        tmp.probes(:,:,ii,:) = permute(pr{ii}, [2 1 3 4]);
    end
    
    % positions
    tmp.positions = transpose(tmp.positions);
    tmp.positions_real = transpose(tmp.positions_real);
    tmp.positions_orig = transpose(tmp.positions_orig);
    
    % ctr
    tmp.ctr = transpose(tmp.ctr);
            
end

