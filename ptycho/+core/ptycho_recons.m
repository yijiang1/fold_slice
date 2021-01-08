% [OUT, STATUS] = PTYCHO_RECONS(P, PREPARE_ONLY = false)
% Runs the reconstruction using parameters in the structure p.
% Returns a structure OUT containing all necessary information.
% 
% ** p                  p structure
% 
% *optional*    
% ** prepare_only       stop after the data preparation; default: false
%
% returns:
% out                   updated p structure
%
% see also: template_ptycho

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

function [out,status] = ptycho_recons(p, prepare_only)

if ~exist('+math/argmax.m','file')
    if exist(p.cSAXS_matlab_path, 'dir')  %% avoid repeated loading by testing availibility of argmax.m file
        addpath(p.cSAXS_matlab_path);
    elseif ~isempty(p.cSAXS_matlab_path)
        warning('Nonexistent cSAXS_matlab_path: "%s"', p.cSAXS_matlab_path)
    end
end
import utils.*

if nargin < 2
    prepare_only = false;
end

caller = dbstack;
p.caller = caller(end).name;

if ~isfield(p, 'queue')
    p.queue = struct();
end
if ~isfield(p, 'io')
    p.io = struct();
end

if ~isfield(p.io, 'SMS_sleep')
    p.io.SMS_sleep = 1800;
end

if ~isfield(p.io, 'phone_number')
    p.io.phone_number = [];
end

if ~isfield(p.io, 'send_failed_scans_SMS')
    p.io.send_failed_scans_SMS = false;
end
if ~isfield(p.io, 'send_finished_recon_SMS')
    p.io.send_finished_recon_SMS = false;
end
if ~isfield(p.io, 'send_crashed_recon_SMS')
    p.io.send_crashed_recon_SMS = false;
end

utils.verbose(struct('prefix', {'init'}))
p.   run_name = ''; 

p.getReport.completed = false;

p = core.ptycho_prepare_paths(p, true);

if isfield(p.queue,'path')&&~isempty(p.queue.path) && ~isfield(p.queue, 'name')
    verbose(0,'Missing setting of p.queue.name, using default p.queue.name=''filelist'' ')
    p.queue.name = 'filelist';
end

if ~isfield(p.queue, 'isreplica')
    p.queue.isreplica = false;
end
if ~isfield(p.queue, 'remote_recons')
    p.queue.remote_recons = false;
end

p.recon_success = false; %% added by YJ

% the existence of this file will cancel the calculation, 
% CTRL-C replacement, checked each iteration
p.io.break_check_name = '/tmp/break_ptycho';
% Function for clean exit, currently used to exit matlab to clean memory
% from MEX upon Ctrl-c
% c = onCleanup(@()ptycho_exit);

% Check screen size
try
    p.plot.scrsz = get(0,'ScreenSize');
catch
    p.plot.scrsz = [1 1 2560 1024];
end

finishup = utils.onCleanup(@(x) ptycho_exit(x), p);

    function ptycho_call()
        import utils.*

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check for file queue  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        [p, status] = scans.get_queue(p, false);
        if ~status
            return
        end
        
        p.   run_name = [p.prefix core.generate_scan_name(p) '_' num2str(p.asize(1)/2^p.detector.binning) 'x' num2str(p.asize(2)/2^p.detector.binning) '_b' num2str(p.detector.binning) '_' p.suffix];  % If empty: automatically generated

        finishup.update(p);
        
        if ~p.queue.remote_recons || p.queue.isreplica
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % initialize ptycho and prepare object and probe
            [p, status] = core.initialize_ptycho(p); %p.positions are created here. unit: pxiel
            
            if ~status || prepare_only
                finishup.update(p);
                out = {p};
                return;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            verbose(struct('prefix', {'ptycho'}))
            verbose(0, ['Reconstructing ' repmat('S%05d ', 1, numel(p.scan_number))], p.scan_number)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%  MAIN PTYCHOGRAPHY CODE %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % get all engines
            verbose(struct('prefix', {'ptycho'}))
            for ieng=1:length(p.engines)
                etic = tic();
                p.current_engine_id = ieng;
                
                % engine call
                verbose(1, 'Calling engine %s', p.engines{ieng}.name)
                verbose(struct('prefix', {p.engines{ieng}.name}))
                
    
                [p, fdb] = core.run_engine(p,ieng);
                if fdb.status.status ~= 0
                    error('Engine %s returned with exit status %d from %s [%d].\n', p.engines{ieng}.name, fdb.status.status, fdb.status.ln(1).name, fdb.status.ln(1).line);
                end
                
                fdb = [];
                
                if p.ortho_probes && size(p.probes,4)>1
                    % orthogonalize probes
                    p.probes = core.probe_modes_ortho(p.probes);
                end
                
                % save reconstructed object, probe and feedback in the p structure of
                % the currently used engine
                p.engines{ieng}.object_final = p.object;
                p.engines{ieng}.probes_final = p.probes;
                p.engines{ieng}.error_metric_final = p.error_metric;
                
                % store images of current engine
                if p.save.store_images_intermediate
                    p = core.save.save_results(p, 0);
                end
               
                if ieng~=length(p.engines)  && p.use_display
                    % intermediate results, not yet final plotting 
                    core.analysis.plot_results(p);
                end
                
                etoc = toc(etic);
                verbose(struct('prefix', {'ptycho'}))
                verbose(1, 'Elapsed time for engine %s: %0.1f s', p.engines{ieng}.name, etoc)
                
            end
            

            verbose(struct('prefix', {'saving'}))
            try
                p = core.save.save_results(p, 1);
            catch ME
                if p.verbose_level > 3
                    keyboard
                else
                    disp('#######  Failed to save data. ##########');
                    rethrow(ME)
                end
            end
                
        else
            % remote reconstruction
%             p.queue.isreplica = false;
            p = core.export4remote(p);
            finishup.update(p);
            [p, status] = core.remote_status(p);
            if ~status
                core.remote_cleanup(p, status);
                finishup.update(p);
                out = {p};
                return;
            end
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if p.queue.remote_recons
            status = core.remote_cleanup(p, true);
        end
        
        [p, status] = scans.get_queue(p, true);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        p.getReport.completed = true;
        
        finishup.update(p);
    end

if p.verbose_level > 3
    fprintf('\n###############################################################\n')
    fprintf('########################## DEBUG MODE #########################\n')
    fprintf('###############################################################\n')
    ptycho_call();
else
    try
        ptycho_call();
        p.recon_success = true; %% added by YJ

    catch ME
        p.getReport.ME = ME;
        p.getReport.crashed = true;
        finishup.update(p);
        status = false;
        p.recon_success = false; %% added by YJ

    end
end



out = p;
return

end


