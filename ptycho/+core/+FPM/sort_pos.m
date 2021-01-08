%SORT_POS sort positions
% par = sort_pos(par, varargin)
%
% sort_pos is a helper function of FP_prealign
% it minimizes the path length, similar to the travelling salesman problem,
% albeit optimized for the peculiarities of a Fourier ptychographic scan
%
% ** par            FP_prealign structure
%
% *optional*
% ** sort_type      'raster', 'round' or 'raster_lim'
%
% returns:
% ++ par            updated FP_prealign structure
%
%
% see also: core.FPM.FP_prealign

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

function par = sort_pos(par, varargin)
pos = par.pos;
sel_axis = par.axis;

if nargin > 1
    sort_type = varargin{1};
else
    sort_type = 'raster';
end

if strcmp(sort_type, 'raster')

    % first find highest and lowest measurements
    pos_min = min(pos(:,1));
    pos_max = max(pos(:,1));
    offset = 2e-6;
    stepsz = 10e-6;
    
    % now select rows of with width of stepsz
    
    stpindx = 0;
    offset = offset + stepsz;
    row = {};
    while true
        lb = pos_min-offset + stpindx*stepsz;
        tb = pos_min-offset + (stpindx+1)*stepsz;
        if pos_max-tb<stepsz
            tb = tb + stepsz;
        end
        if lb >= pos_max
            break;
        end
        row{stpindx+1} = [];
        for i=1:size(pos,1)
            
            if (pos(i,sel_axis)>=lb) && (pos(i,sel_axis)<tb)
                row{stpindx+1} = [row{stpindx+1} i];
                
            end
        end
        
        stpindx = stpindx + 1;
        
    end
    
    % get pos in correct order for alignment
    ascend=true;
    alid = [];
    if sel_axis==1
        axsort = 2;
    else
        axsort=1;
    end
    for i=1:size(row,2)
        if isempty(row{i})
            continue
        else
            if ascend
                direction = 'ascend';
            else
                direction = 'descend';
            end
            [~,I] = sort(pos(row{i}, axsort), direction);
            row_sel = row{i};
            alid = [alid row_sel(I)];
            ascend = ~ascend;
        end
        
    end
elseif strcmp(sort_type, 'round')
    clear alid;
    alid{1} = [];
    if sel_axis==1
        ascend = true;
    else
        ascend = false;
    end
    dr = 10e-6;
    radii = sqrt(par.pos(:,1).^2 + par.pos(:,2).^2);
    
    rad_max = max(radii);
    rad_min = min(radii);
    
    stepindx = 0;
    shell = {};
    last_run = false;
    % find positions in shell
    while true
        lb = rad_min + stepindx*dr;
        tb = rad_min + (stepindx+1)*dr;
        if rad_max-tb<dr
            tb = tb + dr;
            last_run = true;
        end

        shell{stepindx+1} = [];
        for ii=1:size(par.pos,1)
            if radii(ii)>=lb && radii(ii)<tb
                shell{stepindx+1} = [shell{stepindx+1} ii];
            end
        end
        stepindx = stepindx + 1;
        if last_run
            break;
        end
    end
    
    % get pos in correct order for alignment
    for shindx=1:size(shell,2)
        if isempty(shell{shindx})
            fprintf('Warning: Empty shell in path optimization!')
            continue
        else
            if ascend
                direction = 'ascend';
            else
                direction = 'descend';
            end
        [~,I] = sort(atan2(par.pos(shell{shindx},1),par.pos(shell{shindx},2)), direction);
        shell_sel = shell{shindx};
        alid{1} = [alid{1} shell_sel(I)];
        ascend = ~ascend;
        
        end
    end
    
    par.align_section = [];
    par.align_section{1} = 1:par.sz(3)-1;
    par.align_section{1} = par.align_section{1}';
    
elseif strcmp(sort_type, 'raster_lim')
    

    clear alid;
    alid{1} = [];
    
    % split positions into 4 subsections
    par.align_section{1} = find(par.pos(:,mod(sel_axis,2)+1)<-par.rad_filt_min);
    par.align_section{2} = find(par.pos(:,mod(sel_axis,2)+1)>par.rad_filt_min);
    par.align_section{3} = find(par.pos(:,mod(sel_axis+1,2)+1)<-par.rad_filt_min);
    par.align_section{4} = find(par.pos(:,mod(sel_axis+1,2)+1)>par.rad_filt_min);
    
    offset = 2e-6;
    stepsz = 10e-6;
    
    
    
    % now select rows of with width of stepsz
    for jj=1:length(par.align_section)
        
        % first find highest and lowest measurements
        pos_min = min(pos(par.align_section{jj},sel_axis));
        pos_max = max(pos(par.align_section{jj},sel_axis));
        stpindx = 0;
        offset = offset + stepsz;
        row = {};
        last_run = false;
        while true
            lb = pos_min-offset + stpindx*stepsz;
            tb = pos_min-offset + (stpindx+1)*stepsz;
            if pos_max-tb<stepsz
                tb = tb + stepsz;
                last_run = true;
            end
            row{stpindx+1} = [];
            for ii=par.align_section{jj}'
                
                if (pos(ii,sel_axis)>=lb) && (pos(ii,sel_axis)<tb)
                    row{stpindx+1} = [row{stpindx+1} ii];
                    
                end
            end
            
            stpindx = stpindx + 1;
            
            if last_run
                break;
            end
            
        end
        
        % get pos in correct order for alignment
        ascend=true;
        alid{jj} = [];
        if sel_axis==1
            axsort = 2;
        else
            axsort=1;
        end
        for ii=1:size(row,2)
            if isempty(row{ii})
                continue
            else
                if ascend
                    direction = 'ascend';
                else
                    direction = 'descend';
                end
                [~,I] = sort(pos(row{ii}, axsort), direction);
                row_sel = row{ii};
                alid{jj} = [alid{jj} row_sel(I)];
                ascend = ~ascend;
            end
            
        end
    end
    
end

if par.plot_alignment
    fig1 = plotting.smart_figure(1);
    clf;
    hold on
    plot(pos(:,1), pos(:,2))
    title('Alignment')
    colors = jet(length(par.align_section));
    for jj=1:length(par.align_section)
        for ii=1:size(alid{jj},2)-1
            set(groot,'CurrentFigure',fig1);
            plot(pos(alid{jj}(ii),1), pos(alid{jj}(ii),2), 'Color', colors(jj,:), 'Marker', 'x')
            pause(0.01)
        end
    end
    hold off
end
par.alid = [];
par.alid = alid;

end
