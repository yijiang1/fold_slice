% function [varargout] = banana_scan(basepath, scan, Nx, Ny)
%
% This function only works for mesh scans done to see the undulator banana
% with spec. This is an example pf such scans:
% dmesh idgap -0.05 0.05 50 sl1cv -0.4 0.4 8 0.2
% 50 steps (51 points) along the fast axis
% 8 steps (9 points) along the slow axis
%
%
% Input parameters:
%
% basepath: main path where you are working
% scan: scan number of spec mesh scan
% Nx: number of points in fast axis
% Ny: number of steps in slow axis
%
% please reports bugs, problems, suggestions for improvements to:
% CXS group
%
% the case of duplicate scannumbers in the spec file, and how to
%   address them remains to be implemented

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




function [varargout] = banana_plot(basepath, scan, Nx, Ny)

addpath([basepath 'matlab/'])
import io.*

S=io.spec_read(basepath,'ScanNr',scan);
m=reshape(S.bpm4i,Nx,Ny);
idgap=S.idgap(1:Nx);
m_sl1cv=reshape(S.sl1cv,Nx,Ny);
sl1cv=m_sl1cv(1,:);

imagesc(idgap,sl1cv,m')
title(sprintf('scan %05d',scan));
xlabel('idgap')
ylabel('sl1cv')

end

