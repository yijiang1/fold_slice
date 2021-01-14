% [varargout] = find_capillary(varargin)

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

function [varargout] = find_capillary(varargin)
import io.spec_read

    vararg = varargin;
    vararg{end+1} = 'UnhandledParError';
    vararg{end+1} = 0;
    [S] = spec_read(vararg{1},vararg{2:end});
    vararg = vararg(2:end);
    for jj = 1:2:length(vararg)
        name = vararg{jj};
        value = vararg{jj+1};
        switch name
            case 'Counter'
            counter = S.(value);
        end
    end
    
    arrout = regexp(S.S,' +','split');
    motor = S.(arrout{4});
    threshold = .1 * max(counter);
    
    motor = motor(counter>threshold);
    counter = counter(counter>threshold);
    

    threshold = .9 * max(counter);
    i_i = find(counter>threshold,1,'first');
    i_f = find(counter>threshold,1,'last');
    
    p = polyfit(motor(counter>threshold),counter(counter>threshold),1);
    dy = polyval(p,motor)-counter;
    COM = sum(motor(i_i:i_f).*dy(i_i:i_f))/sum(dy(i_i:i_f));
    
    do_plot = 0;
    if (do_plot)
        figure(1)
        plot(motor,dy)
        hold on
        area(motor(i_i:i_f),dy(i_i:i_f))
        plot([1 1]*COM,[0 max(dy(i_i:i_f))],'r','LineWidth',2)
        hold off
    end
    
    varargout{1} = COM;
    
    

end
