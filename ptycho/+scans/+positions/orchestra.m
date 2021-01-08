%OMNY Load positions from Orchestra scan file
function [ p ] = orchestra( p )
import beamline.*
import utils.*

if isempty(p.positions_file)
    error('OMNY positions file is not specified. Please check p.positions_file in your template.')
end

if ~isfield(p,'angular_correction_setup') || isempty(p.angular_correction_setup)
    error('p.angular_correction_setup is not specified. Please check p.angular_correction_setup in your template.')
end

if isfield(p,'omny_interferometer')
    error(' p.omny_interferometer is not supported, use p.angular_correction_setup')
end

if ~isfield(p.detector,'burst_frames')||isempty(p.detector.burst_frames)
    p.detector.burst_frames = 1;
end

switch lower(p.angular_correction_setup)
    case 'omny'
        p.   orchestra.laser_height=-10.0e-3;                                 % Height of horizontal laser beam on the sphere compared to pin tip (only for p.fromspec='opos_angle', 13.5e-3 for OMNI (not fully tested, better with opos than opos_angle), -10.0e-3 for OMNY)
        p.   orchestra.mirrdis=-9.0e-3;                                       % Distance mirror-pin tip (only for p.fromspec='opos_angle', 22.0e-3 for OMNI (not fully tested, better with opos than opos_angle), -9.0e-3 for OMNY)
        p.   orchestra.beam_separation=7.5e-3;                                % Distance mirror-pin tip (only for p.fromspec='opos_angle', 13.0e-3 for OMNI (not fully tested, better with opos than opos_angle), 7.5e-3 for OMNY)
        apply_correction = true; 
    case 'flomni'
        p.   orchestra.laser_height=-13.5e-3;                                  % Height of horizontal laser beam on the sphere compared to pin tip (only for p.fromspec='opos_angle', 13.5e-3 for OMNI (not fully tested, better with opos than opos_angle), -10.0e-3 for OMNY)
        p.   orchestra.mirrdis=-17.4e-3;                                       % Distance mirror-pin tip (only for p.fromspec='opos_angle', 22.0e-3 for OMNI (not fully tested, better with opos than opos_angle), -9.0e-3 for OMNY)
        p.   orchestra.beam_separation=-16e-3;                                 % Distance mirror-pin tip (only for p.fromspec='opos_angle', 13.0e-3 for OMNI (not fully tested, better with opos than opos_angle), 7.5e-3 for OMNY)
        apply_correction = true; 
    case {'lamni', 'none'}
        apply_correction = false;
    otherwise
        error('Wrong  p.angular_correction_setup, choose from ''omny'', ''flomni'',''lamni'',''none'' ')
end



for ii = 1:length(p.scan_number)
    p.scan.is_cont = true;  % So that burst data is prepared normally rather than integrated
    if ~exist(sprintf(p.positions_file,p.scan_number(ii)), 'file' )
        error('Missing OMNY specs file  %s', sprintf(p.positions_file,p.scan_number(ii)))
    end
    out_orch = read_omny_pos(sprintf(p.positions_file,p.scan_number(ii)));
    if ~isfield(out_orch,'Average_y_st_fzp') ||  ~isfield(out_orch,'Average_x_st_fzp')
        out_orch.Average_y_st_fzp = out_orch.Average_y;
        out_orch.Average_x_st_fzp = out_orch.Average_x;
    end
    if ~isfield(out_orch, 'Average_rotz_st')
        apply_correction =false;
    end
    if ~apply_correction
        if isfield(out_orch, 'Average_y_st_fzp')
            positions_real = [out_orch.Average_y_st_fzp*1e-6 out_orch.Average_x_st_fzp*1e-6];
        else  % outdated position format
            positions_real = [out_orch.Average_y*1e-6 out_orch.Average_x*1e-6];
        end
    else
        deltax = p.orchestra.laser_height*out_orch.Average_rotz_st*1e-6/p.orchestra.beam_separation; % p.orchestra.beam_separation: separation between two laser beams for angular measurement
        % p.orchestra.laser_height: height of horizontal laser beam on the sphere compared to pin tip
        deltay = p.orchestra.mirrdis*out_orch.Average_rotz_st*1e-6/p.orchestra.beam_separation; % p.orchestra.beam_separation: separation between two laser beams for angular measurement
        % p.orchestra.mirrdis dist mirror-pin tip
        posx = out_orch.Average_x_st_fzp*1e-6 - deltax;
        posy = out_orch.Average_y_st_fzp*1e-6 - deltay;
        positions_real = [posy posx];
    end
    
    p.numpts(ii) = size(positions_real,1)*p.detector.burst_frames;
    
    positions_tmp = zeros(p.numpts(ii), 2);
    positions_tmp(:,1) = reshape(repmat(positions_real(:,1)',[p.detector.burst_frames 1]),[],1); 
    positions_tmp(:,2) = reshape(repmat(positions_real(:,2)',[p.detector.burst_frames 1]),[],1); 
    
    p.positions_real = [p.positions_real ; positions_tmp];
    %size(p.positions_real)
    
end

end

