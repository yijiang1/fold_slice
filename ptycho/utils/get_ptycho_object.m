function [object] = get_ptycho_object(reconDir, varargin)
%Get object function from a ptychographic reconstruction
%   Input: reconDir-dir to the recon file (.h5)
%   Output: object-complex object function
if nargin==1
    objectPath = '/reconstruction/p/objects/object_0';
else
    objectPath = varargin{1};
end
h = h5read(reconDir,objectPath);
object = h.r + 1i*h.i;

end

