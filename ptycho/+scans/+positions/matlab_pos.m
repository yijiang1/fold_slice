%MATLAB_POS calculate scan parameters based on the values set in the
%template

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
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
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


function [ p ] = matlab_pos( p )

for ii = 1:p.numscans
    positions_real = zeros(0,2); 
    switch p.scan.type
        case 'raster'
            %for iy=0:p.scan.ny
            %    for ix=0:p.scan.nx
            for iy=1:p.scan.ny %modified by YJ. seems odd to begin with 0...
                for ix=1:p.scan.nx
                    xy = [iy * p.scan.step_size_y, ix *  p.scan.step_size_x] + ...
                            randn(1,2).*p.scan.step_randn_offset.*[ p.scan.step_size_y,  p.scan.step_size_x];
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
            
        case 'round'
            dr = (p.scan.radius_out - p.scan.radius_in)/ p.scan.nr;
            for ir=1:p.scan.nr+1
                rr = p.scan.radius_in + ir*dr;
                dth = 2*pi / (p.scan.nth*ir);
                for ith=0:p.scan.nth*ir-1
                    th = ith*dth;
                    xy = rr * [sin(th), cos(th)];
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
        case 'round_roi'
            rmax = sqrt((p.scan.lx/2)^2 + (p.scan.ly/2)^2);
            nr = 1 + floor(rmax/p.scan.dr);
            for ir=1:nr+1
                rr = ir*p.scan.dr;
                dth = 2*pi / (p.scan.nth*ir);
                for ith=0:p.scan.nth*ir-1
                    th = ith*dth;
                    xy = rr * [sin(th), cos(th)];
                    if( abs(xy(1)) >= p.scan.ly/2 || (abs(xy(2)) > p.scan.lx/2) )
                        continue
                    end
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
            
        case 'fermat'
            % this should be changed to have the same variable
            % conventions as in its spec implementation
            phi=2*pi*((1+sqrt(5))/2.) + p.scan.b*pi;
            start = 1;
            if ~isempty(p.scan.lx)
                for ir=start:p.scan.n_max
                    r=p.scan.step*0.57*sqrt(ir);
                    if abs(r*sin(ir*phi))> p.scan.ly/2
                        continue
                    end
                    if abs(r*cos(ir*phi))> p.scan.lx/2
                        continue
                    end
                    xy  = [r*sin(ir*phi)+p.scan.cenxy(1) r*cos(ir*phi)+p.scan.cenxy(2)];
                    positions_real(end+1,:) = xy;
                end
            else
                for ir=start:p.scan.n_max
                    r=p.scan.step*0.57*sqrt(ir);
                    xy  = [r*sin(ir*phi)+p.scan.cenxy(1) r*cos(ir*phi)+p.scan.cenxy(2)];
                    positions_real(end+1,:) = xy;
                end
            end
            
            
            
        case 'custom'
            fn_splt = strsplit(p.scan.custom_positions_source,'.');
            if length(fn_splt)>1
                % file already has an extension
                ext = fn_splt(end);
                if strcmp(ext, 'm')
                    [~, positions_real, ~] = p.scan.custom_positions_source(p);
                elseif strcmp(ext, 'mat')
                    posi = load(p.scan.custom_positions_source, 'pos');
                    positions_real = posi.pos;
                    clear posi;
                else
                    error('File extenstion %s is not supported.', ext)
                end
                
            else
                % file does not have an extension
                if exist([p.scan.custom_positions_source '.m'], 'file')
                    [~, positions_real, ~] = p.scan.custom_positions_source(p);
                elseif exist([p.scan.custom_positions_source '.mat'], 'file')
                    posi = load(p.scan.custom_positions_source, 'pos');
                    positions_real = posi.pos;
                    clear posi;
                else
                    error('Could not find function or data file %s', p.scan.custom_positions_source);
                end
            end
            
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real];
end
    
end

