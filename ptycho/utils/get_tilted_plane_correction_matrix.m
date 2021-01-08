%GET_TILTED_PLANE_CORRECTION_MATRIX 
%   create a sparse matrix that applies deformation on measured farfield to correct for effects of
%   sample tilt 
% 
%   T =  get_tilted_plane_correction_matrix(Npix, detectorDistance,detectorPixel,Chi,Theta,Psi)
%
% Inputs: 
%   **Npix          - size of the dataset 
%   **detectorDistance - sample to detector distance 
%   **detectorPixel - size of detector pixel 
%   **Theta         - rotation around X axis 
%   **Chi           - rotation around Y axis 
%   **Psi           - rotation around beam axis 

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


function T =  get_tilted_plane_correction_matrix(Npix, detectorDistance,detectorPixel,Theta,Chi,Psi)

  
    Npix = Npix(1); %% hardcoded assumption of square image 
    
    [qxx,qyy] = tilted_ewald_plane(Npix, Npix,detectorDistance,detectorPixel, Theta,Chi, Psi);
    
    % transform coordinates to regular grid 1:Npix
    qxx = qxx' / cosd(Chi) +Npix/2+1; 
    qyy = qyy' / cosd(Theta) +Npix/2+1; 
 
 
    M1 = [00,0;0,1;1,0;1,1];
    W = zeros( Npix, Npix, 4, 'single');
   
    pos = cat(3, qxx, qyy); 
    pos_ = floor(pos);
    dX = 1 - (pos - pos_);
    clear pos
    dX = dX(:,:,[2,1]);    
    W(:,:,1) = dX(:,:,1) .* dX(:,:,2);
    W(:,:,2) = (1-dX(:,:,1)) .* dX(:,:,2);
    W(:,:,3) = dX(:,:,1) .* (1-dX(:,:,2));
    W(:,:,4) = (1-dX(:,:,1)) .* (1-dX(:,:,2));
    clear dX 
      
    
    ind_wrong = any(pos_<=0,3) | pos_(:,:,1) > Npix | pos_(:,:,2) > Npix ;
    W = W .* ~ind_wrong;
    
    clear ind_wrong   
    for j = 1:2
        P{j} = zeros( Npix, Npix,4, 'single');
        for i = 1:4
            P{j}(:,:,i) = pos_(:,:,j) + M1(i, j);
        end
    end
    clear pos_
    [Y,X] = meshgrid(1:Npix, 1:Npix);
    ind_x = single((X-1)*Npix+Y);
    ind_x = repmat(reshape(ind_x, [Npix, Npix]), [1,1,4]);
    ind_y = (P{1}-1)*Npix+P{2};
    out = P{1} > Npix | P{2} > Npix | P{1} < 1 | P{2} < 1; 
    S = [ind_x(:), ind_y(:), W(:)];
    S(out(:),:) = []; 
    
    T = sparse(double(S(:,1)), double(S(:,2)), double(S(:,3)),  Npix*Npix, Npix*Npix );

    % renormalize the deformation to preserve flux !!
    T = T ./ max(eps, sqrt(sum(T'*T,1)));


end

function [qyy,qxx] = tilted_ewald_plane(dimx, dimy,detectorDistance,detectorPixel, Theta,Chi,Psi)

xgrid = linspace(-dimx/2,dimx/2, dimx);
ygrid = linspace(-dimy/2,dimy/2, dimy);

[y, x] = meshgrid(xgrid,ygrid);


x = x * detectorPixel;
y = y * detectorPixel;

R = sqrt(x.^2 + y.^2 + detectorDistance^2); 
qx = x ./ R; 
qy = y ./ R; 
qz = (detectorDistance ./ R -1); 


RotMat = utils.get_rotation_matrix_3D(Chi, Theta, -Psi); % make rotation consistent with imrotate in matlab


rotCoords = RotMat*[qx(:)';qy(:)';qz(:)'];% rotate Ewald sphere 
qxx = rotCoords(1,:);
qyy = rotCoords(2,:);
qxx =reshape(qxx, dimx, dimy); 
qyy =reshape(qyy, dimx, dimy); 
% the qzz coordinates are lost 

%calculate Cartesian FFT of model to use for interpolation
% [Y, X] = meshgrid(xgrid,ygrid);


%convert coordinates to correct units
% X = cosd(Theta)*X;
% Y = cosd(Chi)*Y;
scale = detectorPixel / detectorDistance; 
qxx = qxx / scale; 
qyy = qyy / scale; 

% imagesc(X)

% if method == 1
%     img = interp2(Y,X,double(img),qyy,qxx);
% else
%     img = griddata(qyy,qxx,double(img),Y,X);
% end

end

