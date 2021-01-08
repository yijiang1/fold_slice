%INITIALIZE_PTYCHO 
% everything that needs to be done before triggering the reconstruction. This includes 
% inter alia initial checks, intial guess preparations and loading the data. 
%
% ** p          p structure
%
% returns:
% ++ p          p structure
% ++ status     status flag
%
%
% see also: core.ptycho_recons
%
% Based on cSAXS code, modified by Yi Jiang

function [ p, status ] = initialize_ptycho( p )

import utils.*
import io.*

%%% read meta data %%%
if ~isfield(p, 'src_metadata')
    verbose(0,' p.src_metadata is not set, using default p.src_metadata = ''spec''')
    p.   src_metadata = 'spec';    % load meta data from file; currently only 'spec' is supported;
end

% check store_images flag
if ~isfield(p.save, 'store_images')
    p.save.store_images = true;
end
if ~p.save.store_images
    close all
end

% prepare container for meta data
assert( isnumeric(p.scan_number), 'p.scan_number has to contain an integer number')
p.numscans = length(p.scan_number); % Number of scans
p.meta = cell(1,length(p.scan_number));

p = scans.read_metadata(p);

for ii = 1:p.numscans
    p.   scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
end

% write procID
write_procID(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checks and defaults %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = core.initial_checks(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% prepare paths, note that it was already initialized in ptycho_recons %%%

p = core.ptycho_prepare_paths(p);


%%% load detector settings %%%
p = detector.load_detector(p);


if isfield(p, 'ds') && ~isempty(p.ds)
    warning(['Defining ds in the template is not supported anymore and ' ...
    'will not change the pixel size. Please make sure '...
    'that it is set correctly in +detector/+%s/%s.m and remove ds from your template.'], p.detector, p.detector)
end
for ii=1:length(p.detectors)
    assert(p.detectors(1).params.pixel_size==p.detectors(ii).params.pixel_size, 'Different detector pixel sizes are not supported at the moment.')
end
p.ds = p.detectors(1).params.pixel_size;


if check_option(p, 'prop_regime', 'nearfield')
    % nearfield ptychography 
    assert(check_option(p,'focus_to_sample_distance'), 'Undefined p.focus_to_sample_distance that is required for nearfield ptychography')
    p.nearfield_magnification = (p.z-p.focus_to_sample_distance)/p.focus_to_sample_distance; 
    verbose(1, 'Propagation in nearfield regime, magnification = %g', p.nearfield_magnification)
    p.dx_spec = [p.ds,p.ds] / p.nearfield_magnification;
    p.z = p.z / p.nearfield_magnification;
else
    % standard farfield ptychography
    % modified by YJ for electron pty
    if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
        p.dx_spec = 1./p.asize/p.dk; %angstrom
    else
        p.dx_spec = p.lambda*p.z ./ (p.asize*p.ds);                   % resolution in the specimen plane
    end
    p.dx_spec = p.dx_spec ./ cosd(p.sample_rotation_angles(1:2));      % account for a tilted sample ptychography      
end


%%% prepare positions %%%
p = scans.read_positions(p);

%%% find which positions belongs to each object 
p = core.find_shared_IDs(p); 


%%% prepare positions
% Prepare positions, note the output is already in probe positions which
% are different from object (scan) positions by a minus sign
p = core.ptycho_adjust_positions(p);
% p.positions_orig = p.positions;
% p.numpts_orig = p.numpts;
p.numpos = sum(p.numpts);


p.asize_nobin = p.asize;

%%%%%%%%%%%%%%%%%%%%%%
%%% prepare scans %%%%
%%%%%%%%%%%%%%%%%%%%%%


% make sure that all scans have a ctr
numctr = size(p.ctr,1);
if p.numscans > numctr
    for ii=numctr+1:p.numscans
        p.ctr(end+1,:) = p.ctr(numctr,:);
    end
elseif p.numscans < numctr
    p.ctr(p.numscans+1:end,:) = [];
end

%%%%  load data, mask and generate initial estimate of the probe
if p.prepare.auto_prepare_data
    [p, status]=core.ptycho_prepare_scans(p);
else
    if ~isa(p.prepare.prepare_data_function, 'function_handle')
        error(['Expected function handle as p.prepare.prepare_data_function. '... 
            'Please update p.prepare.auto_prepare_data or set p.prepare.auto_prepare_data=true.'])
    else
        [p, status] = p.prepare.prepare_data_function(p);
    end
end
if ~status
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot prepared data %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.plot.prepared_data && p.use_display
    core.analysis.plot_raw_data(p)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot initial guess %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define combined strings for figure title
p.plot.obtitlestring = '';
p.plot.prtitlestring = '';
p.plot.errtitlestring = '';
if p.share_object
    p.plot.obtitlestring =  [core.generate_scan_name(p) ' '];
end
if p.share_probe
    p.plot.prtitlestring = [core.generate_scan_name(p) ' '];
end
p.plot.errtitlestring = [core.generate_scan_name(p) ' '];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot initial guess %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.use_display
    p.plot.extratitlestring = sprintf(' (%dx%d) - Initial guess', p.asize(2), p.asize(1));
    core.analysis.plot_results(p, 'use_display', p.use_display);
end
p.plot.extratitlestring = sprintf(' (%dx%d)', p.asize(2), p.asize(1));

if ~isfield(p.plot, 'windowautopos')
    p.plot.windowautopos = false; % So resizing after first time display is respected
end

verbose(1, 'Finished data preparation and initialization.')

end




