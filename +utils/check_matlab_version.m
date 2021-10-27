%CHECK_MATLAB_VERSION check matlab version to make sure it is compatible
% ver...    compatible version number in STRING
%
% EXAMPLE:
%    check_matlab_version('9.2')
%
% MATLAB 9.0 - 2016a
% MATLAB 9.1 - 2016b
% MATLAB 9.2 - 2017a
% MATLAB 9.3 - 2017b

% modified by YJ for newer versions

function check_matlab_version( ver )

current_version = version;
ver_str = strsplit(current_version, '.');
ver_input_str = strsplit(ver, '.');

ver_input_num = [str2double(ver_input_str{1}), str2double(ver_input_str{2})];
ver_num = [str2double(ver_str{1}), str2double(ver_str{2})];

if ver_num(1) == ver_input_num(1)
    if ver_num(2) < ver_input_num(2)
        warning('You are using Maltab version %d.%02d but the code was designed and tested with %d.%02d.',...
            ver_num(1),ver_num(2),ver_input_num(1),ver_input_num(2));
    end
elseif ver_num(1) < ver_input_num(1)
    warning('You are using Maltab version %d.%02d but the code was designed and tested with %d.%02d.',...
        ver_num(1),ver_num(2),ver_input_num(1),ver_input_num(2));
end

end

