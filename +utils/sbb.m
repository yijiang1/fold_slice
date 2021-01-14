%SBB returns station board for PSI West
% 
% syntax: 
%   as station board:
%       sbb (from, number_of_connections)
%   as connection finder:
%       sbb (from, to, number_of_connections)
%   
%   03/2017

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

function sbb( varargin )
stationboard = true;
if nargin<1
  dep = 'VilligenPSIWest';
else
  dep = varargin{1};
end

if nargin>=2
    if ischar(varargin{2})
        stationboard=false;
        dest = varargin{2};
    end
end
if stationboard
    if nargin <2
        limit = 5;
    else
        limit = varargin{2};
    end
end


server = 'http://transport.opendata.ch/v1';

if stationboard
    req = sprintf('/stationboard?station=%s&limit=%d',dep,limit);
    
    data = webread([server req]);
    
    for i=1:limit
        c_time_stamp = strsplit(data.stationboard(i).stop.departure,'T');
        c_day = c_time_stamp(1);
        c_day = c_day{1};
        c_time = c_time_stamp(2);
        c_time = c_time{1};
        c_time = c_time(1:8);
        stat_dep = data.station.name;
        fprintf('From: %s\n', stat_dep);
        fprintf('To: %s\n', data.stationboard(i).to);
        fprintf('Departure time: %s, %s\n\n', c_time, c_day);
    end
else
    limit = 3;
    req = sprintf('/connections?from=%s&to=%s&limit=%d',dep,dest,limit);
    data = webread([server req]);
    for i=1:limit
        c_data = data.connections(i);
        c_time_stamp = strsplit(c_data.from.departure,'T');
        c_day = c_time_stamp(1);
        c_day = c_day{1};
        c_time = c_time_stamp(2);
        c_time = c_time{1};
        c_time = c_time(1:8);
        
        c_dur = strsplit(c_data.duration, 'd');
        fprintf('From: %s\n', c_data.from.station.name);
        fprintf('To: %s\n', c_data.to.station.name);
        fprintf('Departure time: %s, %s\n', c_time, c_day);
        
        if str2double(c_dur{1})==0
            fprintf('Duration: %s\n\n', c_dur{2})
        else
            fprintf('Duration: %s days %s\n\n', c_dur{1}, c_dur{2});
        end
    end
end
    


end

