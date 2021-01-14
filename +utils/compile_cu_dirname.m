%COMPILE_APS_DIRNAME returns the default APS directory tree for a 
% given scan number
% 
% EXAMPLE:
%   scan_dir = utils.compile_x12sa_dirname(10);
%       -> scan_dir = 'S00000-00999/S00010/'
%
% written by Yi Jiang, based on PSI's code 

function scan_dir = compile_cu_dirname(scan_no)

scan_dir = sprintf('S%05d-%05d/S%05d/',floor(scan_no/1000)*1000, ...
        floor(scan_no/1000)*1000 + 999, ...
        scan_no);
    
end

