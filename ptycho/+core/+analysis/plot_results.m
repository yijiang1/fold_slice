%%% PLOT_RESULTS plotting routine for ptychographic reconstructions
% ** p              p structure from a ptychographic reconstruction
% 
% *optional*
% ** use_display    show plots (default: true)
% ** store_images    write images to disk (default: false)
% ** final          show all error metrics for a final plot (default: false)
% ** save_path      change default save_path (p.save_path) for saving jpgs
%
% EXAMPLES:
%    core.analysis.plot_results(p);
%    core.analysis.plot_results(p, 'store_images', true);
%    core.analysis.plot_results(p, 'use_display', false, 'store_images', true);
%
%
% see also: plotting.ptycho_show_recons
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

function plot_results(p, varargin)
import plotting.*
import utils.*
import math.sp_quantile


% parse p.plot inputs 
check_input = @(x) islogical(x) || isnumeric(x);
parse_p = inputParser;
parse_p.KeepUnmatched = true;

parse_p.addParameter('fourier_ptycho', false, check_input)
parse_p.addParameter('object_spectrum', [], check_input)
parse_p.addParameter('filt', true, check_input)
parse_p.addParameter('plot_layers', true, check_input)
parse_p.addParameter('plot_layers_stack', true, check_input)
parse_p.addParameter('remove_phase_ramp', false, check_input)
parse_p.addParameter('prop_obj',0, check_input)
parse_p.addParameter('plot_conj', false, check_input)
parse_p.addParameter('obj_apod', false, check_input)

parse_p.parse(p.plot);
p.plot = utils.update_param(p.plot, parse_p.Results);

% parse p.save 
parse_p = inputParser;
parse_p.KeepUnmatched = true;
parse_p.addParameter('store_images_format', 'png', @(x)ismember(x, {'png', 'jpg'}))
parse_p.addParameter('store_images_dpi', 150, @math.isint)
parse_p.parse(p.save);
p.save = utils.update_param(p.save, parse_p.Results);


if isempty(varargin) || ischar(varargin{1})
    par = inputParser;
    par.addParameter('use_display', true, check_input)
    par.addParameter('store_images', false, check_input)
    par.addParameter('final', false, check_input)
    par.addParameter('save_path',[], @ischar)
    
    par.parse(varargin{:})
    vars = par.Results;
end


if p.fourier_ptycho
    p.plot.fov_box = false;
end

if isempty(p.plot.object_spectrum)
    p.plot.object_spectrum = (utils.verbose>=3);
end

if ~isfield(p.plot, 'log_scale')
    p.plot.log_scale = [false false];
elseif isscalar(p.plot.log_scale)
    p.plot.log_scale = repmat(p.plot.log_scale,1,2);
end



% Subplot geometry
% p.subplwin = [floor(sqrt(p.numscans)) ceil(p.numscans/floor(sqrt(p.numscans)))];
numwinobj = p.numobjs*p.object_modes;
p.plot.subplwinobj = [floor(sqrt(numwinobj)) ceil(numwinobj/floor(sqrt(numwinobj)))];
if p.numprobs == 1 || p.probe_modes == 1
    % distribute as efficiently as possible 
    numwinprob = p.numprobs*p.probe_modes;
    p.plot.subplwinprob = [floor(sqrt(numwinprob)) ceil(numwinprob/floor(sqrt(numwinprob)))];
else
    % show scans in columns and probe modes in rows 
    p.plot.subplwinprob = [p.probe_modes, p.numprobs];
end
    
    
if check_option(p.plot, 'subplwinobj_dir', 'vertical') || (p.plot.show_layers && size(p.object{1},4) > 1)
    % prefer to stack the object verticaly , useful for multilayer object plotting 
    p.plot.subplwinobj = sort(p.plot.subplwinobj, 'descend');
end


%%%%%%%%%%%%%%%%%%%%
%%%%%% OBJECTS %%%%%
%%%%%%%%%%%%%%%%%%%%


[fig1, fig2] = core.analysis.plot_objects(p, vars.use_display); 

 
%%%%%%%%%%%%%%%%%%%%
%%%%%% PROBES %%%%%%
%%%%%%%%%%%%%%%%%%%%



fig3 = core.analysis.plot_probes(p, vars.use_display);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ERROR METRIC %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


try

     fig4 = core.analysis.plot_error_metric(p, vars.final, vars.use_display); 
        
catch ME
    warning('Failed to plot error metrics.')
    disp([ME.getReport]);
    fig4 = plotting.smart_figure(4);
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% POSITIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

if p.plot.positions
     core.analysis.plot_positions(p); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PROBES @ DETECTOR %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.plot.probe_spectrum
     fig5 = core.analysis.plot_probes_at_detector(p, vars.use_display); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OBJECT SPECTRUM %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.plot.object_spectrum
    fig6 = core.analysis.plot_object_spectrum(p, vars.use_display); 
end


drawnow;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write in figures folder %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you configured to dump it also saves the figures in
% analysis/online/ptycho
if vars.store_images

    
    if isempty(vars.save_path)
        % split the save_path and find the last occurrence of 'analysis'
        tmpPath = strsplit(p.save_path{1}, '/');
        analysisPos = find(strcmpi(tmpPath, 'analysis'));
        save_path_online = strjoin(tmpPath(1:max(1,analysisPos(end)-1)),'/');
        
        if ~isfield(p, 'datasetID')
            p.datasetID = 0;
        end
        
        
        % load sample name if provided in .dat files and append it to the
        % save name suffix 
        if isfield(p, 'samplename')
            suffix = sprintf('dset_%s_%05d', p.samplename, p.datasetID);
        else
            suffix = sprintf('dset_%05d', p.datasetID);
        end
        subdir = fullfile(save_path_online,'analysis/online/ptycho/', suffix);
        gallery = fullfile(save_path_online,'analysis/online/ptycho/gallery/');
        
    else
        subdir = vars.save_path;
        gallery = fullfile(subdir,'/gallery/');
    end
    
    
    if ~exist(gallery,'dir')
        mkdir(gallery);
    end
    if ~exist(subdir,'dir')
        mkdir(subdir);
    end
    
    utils.verbose(0, 'Saving images to %s', subdir)
    
    % ignore prefix in the run name -> make sorting by name equivalent to
    % sorting by scan number -> easier preview and browsing through image
    % gallery 
    image_name = p.run_name(1+length(p.prefix):end); 
    
    width = 6*p.plot.subplwinobj(2);
    height = 4*p.plot.subplwinobj(1); 

    if any(p.save.store_images_ids == 1)
        fig1.PaperPosition = [3 3 width height];  % adjust size of the resulting image 
        save_figs(fig1, '%s_amplitude.%s', image_name, subdir, gallery, p.save);
    end
    if any(p.save.store_images_ids == 2)
        fig2.PaperPosition = [3 3 width height];  % adjust size of the resulting image 
        save_figs(fig2, '%s_phase.%s', image_name, subdir, gallery, p.save);
    end
    if any(p.save.store_images_ids == 6)
        fig6.PaperPosition = [3 3 width height];  % adjust size of the resulting image 
        save_figs(fig6, '%s_object_spectrum.%s', image_name, subdir, gallery, p.save);
    end  

    width = 4*p.plot.subplwinprob(2);
    height = 4*p.plot.subplwinprob(1); 
    if any(p.save.store_images_ids == 3)
        fig3.PaperPosition = [3 3 width height];  % adjust size of the resulting image 
        save_figs(fig3, '%s_probe.%s', image_name, subdir, gallery, p.save);
    end
    if any(p.save.store_images_ids == 4)
        save_figs(fig4, '%s_err.%s', image_name, subdir, gallery, p.save);
    end
    if any(p.save.store_images_ids == 5)
        fig5.PaperPosition = [3 3 width height];  % adjust size of the resulting image 
        save_figs(fig5, '%s_probe_spectrum.%s', image_name, subdir, gallery, p.save);
    end    

end

end

function save_figs(fig_handle, fname,run_name, subdir, gallery, params)
    try

        fname = sprintf(fname,run_name, params.store_images_format);
        utils.verbose(3, 'saving %s',fullfile(subdir,fname));

        switch  params.store_images_format
            case 'png' , printer = '-dpng'; 
            case 'jpg' , printer = '-djpeg'; 
            otherwise, error('Unsupported image extension')
        end
        print(fig_handle, printer,['-r', num2str(params.store_images_dpi)],fullfile(subdir,fname));
        % trim borders around the images 
        system(sprintf('convert -trim  %s %s', fullfile(subdir,fname),  fullfile(subdir,fname)));
        % make a symbolic link to a gallery folder 
        system(sprintf('ln -sf  %s %s', fullfile(subdir,fname), fullfile(gallery, fname)));
    catch err 
        warning('Saving plot handle fig%i failed: %s',fig_handle.Number, err.message)
    end
end

