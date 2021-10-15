function [name] = get_host_name()
%Return host name
%   Written by YJ for I/O
if isunix() 
    name = getenv('HOSTNAME');
else 
    name = getenv('hostname'); 
end

end

