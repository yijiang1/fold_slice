import plotting.plot_radial_integ
import io.spec_read

scan_nr=345;
backgroundscan=343;

base_path='~/Data10/';
eaccount='e15581';
spec_data = spec_read(base_path,'ScanNr',scan_nr);
transmission_data=mean(spec_data.diode);
datafile=sprintf('%sanalysis/radial_integration/%s_1_%05d_00000_00000_integ.mat',base_path,eaccount,scan_nr);
backgroundfile=sprintf('%sanalysis/radial_integration/%s_1_%05d_00000_00000_integ.mat',base_path,eaccount,backgroundscan);
spec_data_bgr = spec_read(base_path,'ScanNr',backgroundscan);
transmission_bgr=mean(spec_data_bgr.diode);
correction=transmission_data/transmission_bgr;

plot_radial_integ(...
datafile,...
'FigNo',101,...             number of the figure for plotting the integrated intensities, default is 100
'NewFig',1,...              open a new figure for each file, default is 0
'ClearFig',0,...            clear the figure before plotting, default is 1
'XLog',1,...                logarithmic scaling of the x-axis, default is 0
'YLog',1,...                logarithmic scaling of the y-axis, default is 1
'PlotQ',1,...               plot as a function of momentum transfer q rather than pixel no., default is 0
'PlotAngle',0,...           plot as a function of the azimuthal angle rather than q or the radius, default is 0
'RadiusRange',[],...        for azimuthal plots the intensity over this radius range is averaged, default is [] for all radii
'QMulPow',[],...        multiply intensity with q to the power of this value, default is [ ] for no multiplication
'Inverse_nm',1,...          plot q in inverse nm rather than inverse Angstroem, default is 0
'PixelSize_mm',0.172,...    pixel size for q-calculation, default is 0.172 mm
'DetDist_mm',7281.9,...     sample to detector distance for q-calculation, default is 7200.000 mm
'E_keV',11.2,...           x-ray energy for q-calculation, default is 11.200 mm
'SegAvg',1,...              average over angular segments rather than plotting them with different line colours, default is 1
'SegRange',[],...           segment range to plot, default is [] for all segments
'LegendMulSeg',1,...        show a legend in case of multiple segments being plotted, default is 1
'PointAvg',1,...            plot the average of all intensity curves in the file, which typically means the average of a scan line, default is 1
'PointRange',[],...         %point range to plot, default is [] for all points in a file
'BgrFilename',backgroundfile,... %         background to subtract from each intensity profile, must have the same dimensions the data have
'BgrScale',correction);        %        scaling factor to apply to the backgroubnd data, default is 1.000e+00
%'Axis',<[ x_from x_to y_from y_to ]>   fixed scale for the plot
%'SleepTime',<seconds>                wait time after each plot, default is 0.000
%'XLog',<0-no, 1-yes>                 logarithmic scaling of the x-axis, default is 0
%'YLog',<0-no, 1-yes>                 logarithmic scaling of the y-axis, default is 1
%'FilenameIntegMasks',<filename>      Matlab file containing the integration masks, needed for normalization in case of averaging over radii, default is '~/Data10/analysis/data/pilatus_integration_masks.mat'
%'BgrFilename',<'filename'>           background to subtract from each intensity profile, must have the same dimensions the data have
%'BgrScale',<value>                   scaling factor to apply to the backgroubnd data, default is 1.000e+00
%'BgrPoint',<integer>                 point to use from the file BgrFilename, default is 1, use [] to subtract 1:1


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