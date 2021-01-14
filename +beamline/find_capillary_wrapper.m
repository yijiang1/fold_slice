% this script contains the necessary loop for 'find_capillary.m' to be called
%    as function of 'spec'
% written by  (last change: 2011-06-16)
%    in case of bugs, problems, and suggestions for improvements, please contact
%    CXS group
%
% note that EPICS communication works only on local machines at the beamline, i.e.,
%    NOT on the compute nodes
% run this (or related scripts that use EPICS for communication), for instance, on
%    x12sa-cons-1

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

lastscan = 0;

while(1)
    scall = sprintf('caget ''X12SA-ES1-DOUBLE-00''');
    [err,io] = system(scall);
    arrout = regexp(io,' +','split');
    scannr = str2double(arrout{2});
    if (scannr > lastscan)
        try
            COM = +beamline.find_capillary('..','ScanNr',scannr,'Counter','diode')
        catch   
            fprintf('Failed find capillary, pausing 5 sec and retrying\n')
            pause(5)
            COM = +beamline.find_capillary('..','ScanNr',scannr,'Counter','diode')
        end
    else
        pause(1);
    end
    scall = sprintf('caputq X12SA-ES1-DOUBLE-02 %f',COM);
    [err,io] = system(scall);
    scall = sprintf('caputq X12SA-ES1-DOUBLE-01 %d',scannr);
    [err,io] = system(scall);
    lastscan = scannr;
end
