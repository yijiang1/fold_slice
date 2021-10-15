function [name] = get_user_name()
%Return account username
%   Written by YJ for I/O
if isunix() 
    name = getenv('USER');
else 
    name = getenv('username'); 
end

end

