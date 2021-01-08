%   [pos] = read_scan_positions_spec(specdatafile,scans,motors)
%   Returns SPEC scan positions.
%   positions = read_scan_positions(specdatafile,scan,motors) returns the 
%   positions of a SPEC file corresponding to two given motor names
%
%   Input parameters: 
%   specdatafile: string with the file name of the SPEC data file 
%                 (e.g.'~/Data10/specES1/dat-files/specES1_started_2017_07_11_1633.dat' )
%   scans: scan number(s) (e.g. [2:77,88])
%   motors: 1x2 cell array with scanning motor names (e.g. {'samx','samy'})
%
%   The output is a structure array with length(scans) number of elements,
%   each element containing in the field 'data' an Nx2 array with N
%   being the number of points in the scan. The two columns contain motor
%   positions for motors samx and samy, respectively, if motors={'samx','samy'}
%
%   Example:  
%   pos=read_scan_positions_spec('~/Data10/specES1/dat-files/specES1_started_2017_07_11_1633.dat',[58:59],{'samx','samy'}); 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function [pos] = read_scan_positions_spec(specdatafile,scans,motors)

Sp = io.spec_read(specdatafile,'ScanNr',scans);

if numel(Sp) == 1
    S{1} = Sp;
else
    S = Sp;
end
clear Sp

field='data';
value=cell(1,length(S));
pos=struct(field,value);

for ii=1:length(S)
    
    C = strsplit(S{ii}.S,' ');
    
    %Is the scan 1D or 2D
    if numel(C)<10
        scanis1D = true;
    else
        scanis1D = false;
    end
    
    fastmotorname = C{4};
        
    % Which is the fast axis?
    if strcmpi(fastmotorname,motors{1})      % Fast axis is x
        fast_axis_index = 1;
        slow_axis_index = 2;
    elseif strcmpi(fastmotorname,motors{2})  % Fast axis is y                                  % Fast axis is y
        fast_axis_index = 2;
        slow_axis_index = 1;
    else
        error(sprintf('The scan motor %s does not match any of the inputs %s %s',fastmotorname,motors{1},motors{2}))
    end
    
    % Is this a continuous scan? Identified by having only one position in
    % the fast axis
    if numel(getfield(S{ii},fastmotorname))==1
        iscont = true;
    else
        iscont = false;
    end
    
num=str2num(C{7})+1;    

    if scanis1D
        positions=zeros(num,2);
        if iscont
            step=(str2num(C{6})-str2num(C{5}))/str2num(C{7});
            positions(:,fast_axis_index) = [str2num(C{5}):step:str2num(C{6})];
        else
            positions(:,fast_axis_index) = getfield(S{ii}, motors{fast_axis_index});
        end
        positions(:,slow_axis_index) = zeros(num,1) + getfield(S{ii}, motors{slow_axis_index});
    else
        num2=str2num(C{11})+1;
        num=num*num2;
        slowmotorname = C{8};
        positions=zeros(num,2);
        if ~strcmpi(slowmotorname,motors{slow_axis_index}) % Check slow axis name
            error(sprintf('The scan motor %s does not match the input %s',slowmotorname,motors{slow_axis_index}))
        end
        positions(:,fast_axis_index) = getfield(S{ii}, motors{fast_axis_index});
        positions(:,slow_axis_index) = getfield(S{ii}, motors{slow_axis_index});
    end
    pos(ii) = struct(field,positions);
end

return;
    