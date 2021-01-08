% IMAGESC_HSV for plotting complex valued arrays , similar to imagesc3D but with more options 
%  imagesc_hsv(varargin)
% 
% ** varargin       see the code 

%     
% Academic License Agreement
% 
% Source Code
% 
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
% 
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
% 
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 


function imagesc_hsv(varargin)

    import utils.*
    import math.*

    par = inputParser;
    par.addOptional('data', [])
    par.addParameter('scale',  nan , @isnumeric )
    par.addParameter('clim',  [] , @isnumeric )
    par.addParameter('inverse',  false , @islogical )  % use white background 
    par.addParameter('show_ROI',  false , @islogical )  % show only intersting area
    par.addParameter('points',  [] , @isnumeric )  % plot dots 
    par.addParameter('enhance_contrast',  false , @islogical )  % plot dots 
    par.addParameter('axis',  [] , @isnumeric )  % plot dots 
    par.addParameter('stabilize_phase',  true , @islogical )  % plot dots 
    par.addParameter('show',  true , @islogical )  % plot dots 

    par.parse(varargin{:})
    r = par.Results;
    data = r.data;
    clim = r.clim;

    if all(data(:) == 0)
        warning('Empty data to plot')
        return 
    end
          
    
    
    [W,H] = size(data);
    
    if ~isempty(r.axis)
        X = linspace(r.axis(1),r.axis(2),W)*1e6;
        Y = linspace(r.axis(3),r.axis(4),H)*1e6;
    else
        if ~isnan(r.scale)
            scale = ones(2,1).*r.scale(:);
            X = [-W/2:W/2-1]* scale(1)*1e6;
            Y = [-H/2:H/2-1]* scale(2)*1e6;
        else
            X = 1:W; Y = 1:H;
        end
    end
    if r.show_ROI
        asum = abs(sum(data,3));
        try           
            T1 =  (graythresh_new((sum(asum,1))));
            T2 =  (graythresh_new((sum(asum,2))));
            asum(:,sum(asum,1) < T1) = 0;
            asum(sum(asum,2) < T2,:) = 0;
            [ROI] = get_ROI(asum > 0.01*quantile(asum(:), 0.99), 0);
            data = data(ROI{:});
            X = X(ROI{1});
            Y = Y(ROI{2});
        catch
            warning('ROI estimation failed')
        end
    end
    [W,H] = size(data);

    if ~isempty(clim)
        ind_min = abs(data) < clim(1); 
        ind_max = abs(data) > clim(2);
        data(ind_min) =  data(ind_min) ./ abs(data(ind_min)) * clim(1);
        data(ind_max) =  data(ind_max) ./ abs(data(ind_max)) * clim(2);
    end
    
    adata = abs(data);
  
    
    
    alpha  = 1e-3;
    tmp= sort(adata(:));
    MAX = tmp(ceil(end*(1-alpha)));
    ind = adata > MAX;
    data(ind) = MAX * data(ind) ./ abs(data(ind));
    if r.enhance_contrast
      data = data ./ sqrt(alpha+abs(data));
      clim = sqrt(clim);
    end
    if r.stabilize_phase
        data = stabilize_phase(data, abs(data), abs(data), 'remove_ramp', false);
    end
        
    adata = abs(data);

    if isempty(clim)
        range = sp_quantile(adata(:), [1e-2, 1-1e-2],10);
    else
        range =  clim;
    end

    adata = (adata  - range(1) ) ./ ( range(2) - range(1) );
    ang_data = angle(data); 
      
    if r.enhance_contrast && r.stabilize_phase
        ang_range = max(abs(sp_quantile(ang_data(:), [1e-2, 1-1e-2],10)));
        ang_range = max(1e-3, ang_range);
        ang_data = 2*pi*ang_data  ./ (2* ang_range);
    end
       
    
    if r.inverse
        hue =  mod(ang_data+1.5*pi, 2*pi)/(2*pi);
        hsv_data = [ hue(:) , adata(:), ones(W*H,1) ];
    else 
        hue =  mod(ang_data+2.5*pi, 2*pi)/(2*pi);
        hsv_data = [ hue(:) , ones(W*H,1), adata(:) ];
    end
    hsv_data = min(max(0, hsv_data),1);


    rgb_data = hsv2rgb(hsv_data);

    rgb_data = reshape(rgb_data, W,H,3);      
    rgb_data = min(1,rgb_data);

    
    if r.show
        hh = imagesc(Y,X, rgb_data );
        axis image 
    end


    if r.show
        % Get the parent Axes of the image
        axis image 

        if ~isempty(r.points)  && ~any(isnan(r.scale))
            hold on
                points = r.scale.*1e6.*r.points;
                plot( points(:,1),points(:,2), '.w')
            hold off
        end
    end
   
    
end



