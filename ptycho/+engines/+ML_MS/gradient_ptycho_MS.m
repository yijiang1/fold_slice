% This function is intended to be only a distributor for either only 'func', or 'func' and 
% 'grad' calculation.

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

function [func, grad, p] = gradient_ptycho_MS(xopt,p,fmag,fmag2,fmask,...
    numobjs,object_size,numprobs,numscans,scanindexrange,initialerror,...
    fnorm,probe_mask,errtitlestring,plot_mask,plot_ind,creg,smooth_gradient)
import utils.fract_hanning
import utils.fract_hanning_pad
import utils.verbose

 
% if p.object_modes ~= 1 || p.probe_modes ~= 1 || numprobs ~= 1 || numobjs ~= 2 || ...
%         p.share_object || p.N_layer ~= 2 || ...
%         p.ms_opt_flags_local(1) ~= 1 || p.ms_opt_flags_local(2) ~= 1 || ~strcmpi(p.opt_errmetric, 'l1')
%     error(['Only double layer single-modes object/single-probe calculations are currently supported. ', ...
%         'E.g. N_layer = 2, numobjs = 2, numprobs = 1, p.object_modes = 1, p.probe_modes = 1, ', ...
%         'p.share_object = false, p.share_probe = true, p.ms_opt_flags_local(1) = true, ', ...
%         'p.ms_opt_flags_local(2) = true, p.opt_errmetric = l1'])
% end

if nargout == 1
    % calculate only the error metric 'func'
    func = feval(@ptycho_MS_err,xopt,p,fmag,fmag2,fmask,...
        numobjs,object_size,numprobs,numscans,scanindexrange,initialerror,...
        fnorm,probe_mask,errtitlestring,plot_mask,plot_ind,creg,smooth_gradient);
    grad = [];
    
else
    % calculate the error metric 'func' *and* gradients 'grad'
    [func, grad] = feval(@ptycho_MS_grad,xopt,p,fmag,fmag2,fmask,...
        numobjs,object_size,numprobs,numscans,scanindexrange,initialerror,...
        fnorm,probe_mask,errtitlestring,plot_mask,plot_ind,creg,smooth_gradient);
end
end

function func = ptycho_MS_err(xopt,p,fmag,fmag2,fmask,...
    numobjs,object_size,numprobs,numscans,scanindexrange,initialerror,...
    fnorm,probe_mask,errtitlestring,plot_mask,plot_ind,creg,smooth_gradient)
import utils.fract_hanning
import utils.fract_hanning_pad
import utils.verbose
N_layer = p.N_layer;
if p.ms_opt_flags_local(3)
    delta_z = xopt(end-N_layer+2:end);
else
    delta_z = p.delta_z;
end

% Precalculate 'prop' and 'derivative_prop' values
k = 2*pi/p.lambda;
prop_der = 1i * k * sqrt(1-(p.lambda*p.Fx).^2-(p.lambda*p.Fy).^2);

prop = cell(1, N_layer-1);
for n = 1:N_layer-1
    prop{n} = exp(delta_z(n)*prop_der);
end
prop = reshape([real(prop{1}(:)).'; imag(prop{1}(:)).'], [], 1);

asize = p.asize;
positions = p.positions;

fmask = logical(fmask); % the mex functions does not check data type!
func = engines.ML_MS.calc_ms_err(xopt, numobjs, N_layer, object_size(1,1), object_size(1,2), ...
    uint64(positions(1:end, 1)), uint64(positions(1:end, 2)), asize(1), asize(2), ...
    prop, 1/fnorm, delta_z, fmask, fmag);

func = func + initialerror;

%verbose(2, ['(Before regul.) Error = ' num2str(func, 20)])

%%% Object regularization %%%
% Normalized regularization to avoid the reduction of object amplitudes
% with the setting of intensity invariant
if (creg > 0) && p.ms_opt_flags_local(1)
    last_ind = 0;
    for obnum = 1:numobjs
        o_size = object_size(obnum, :);
        o_numel = prod(o_size);
        
        for obmode = 1:p.object_modes
            for n = 1:N_layer %% unused??
                ob = reshape(xopt(last_ind+1:2:last_ind+2*o_numel), o_size) + ...
                    1i*reshape(xopt(last_ind+2:2:last_ind+2*o_numel), o_size);
                last_ind = last_ind+2*o_numel;
                
                R = sum(sum(abs(diff(ob, 1, 1)).^2)) + sum(sum(abs(diff(ob, 1, 2)).^2)) - ...
                    sum(abs(diff(ob(:,end), 1, 1)).^2) - sum(abs(diff(ob(end,:), 1, 2)).^2);
                norm_r = sum(sum(abs(ob).^2));
                func = func + creg*R/norm_r;
            end
        end
    end
end

verbose(2, ['delta_z = ' num2str(delta_z*1e6,7) ' um, Error = ' num2str(func)]);


end

function [func, grad] = ptycho_MS_grad(xopt,p,fmag,fmag2,fmask,...
    numobjs,object_size,numprobs,numscans,scanindexrange,initialerror,...
    fnorm,probe_mask,errtitlestring,plot_mask,plot_ind,creg,smooth_gradient)
import utils.verbose
import utils.fract_hanning
import utils.fract_hanning_pad

optimize_object_layer = p.ms_opt_flags_local(1);
optimize_probes = p.ms_opt_flags_local(2);
optimize_delta_z = p.ms_opt_flags_local(3);

N_layer = p.N_layer;
if optimize_delta_z
    delta_z = xopt(end-N_layer+2:end);
else
    delta_z = p.delta_z;
end

verbose(2, 'Computing gradient');

% Precalculate 'prop' and 'derivative_prop' values
k = 2*pi/p.lambda;
prop_der = 1i * k * sqrt(1-(p.lambda*p.Fx).^2-(p.lambda*p.Fy).^2);

prop = cell(1, N_layer-1);
for n = 1:N_layer-1
    prop{n} = exp(delta_z(n)*prop_der);
end

prop = reshape([real(prop{1}(:)).'; imag(prop{1}(:)).'], [], 1);
prop_der = reshape([real(prop_der(:)).'; imag(prop_der(:)).'], [], 1);

asize = p.asize;
positions = p.positions;

fmask = logical(fmask); % the mex functions does not check data type!
assert(isa(xopt, 'double'), 'Inputs has to be double')
[func, grado, gradp, gradz] = engines.ML_MS.calc_ms_grad(xopt,  numobjs, N_layer, object_size(1,1), object_size(1,2), ...
    uint64(positions(:, 1)), uint64(positions(:, 2)), asize(1), asize(2), ...
    prop, prop_der, 1/fnorm, delta_z, fmask, fmag, double(optimize_delta_z));

verbose(2, ['(Before regul.) Error = ' num2str(func)])

%%% Sieves preconditioning %%%
if (any(smooth_gradient(:)) ~= 0) && optimize_object_layer
    for obnum = 1:numobjs
        for obmode = 1:p.object_modes
            for n = 1:N_layer
                grado{obnum,obmode,n} = conv2(grado{obnum,obmode,n}, smooth_gradient, 'same');
            end
        end
    end
end

%%% Object regularization %%%
% Normalized regularization to avoid the reduction of object amplitudes
% with the setting of intensity invariant
if (creg > 0) && optimize_object_layer
    last_ind = 0;
    for obnum = 1:numobjs
        o_size = object_size(obnum, :);
        o_numel = prod(o_size);
        
        for obmode = 1:p.object_modes
            for n = 1:N_layer
                ob = reshape(xopt(last_ind+1:2:last_ind+2*o_numel), o_size) + ...
                    1i*reshape(xopt(last_ind+2:2:last_ind+2*o_numel), o_size);
                last_ind = last_ind+2*o_numel;
                
                diff1 = diff(ob, 1, 1);
                diff2 = diff(ob, 1, 2);
                R = sum(sum(abs(diff1).^2)) + sum(sum(abs(diff2).^2)) - ...
                    sum(abs(diff(ob(:,end), 1, 1)).^2) - sum(abs(diff(ob(end,:), 1, 2)).^2);
                norm_r = sum(sum(abs(ob).^2));
                func = func + creg*R/norm_r;
                
                diff1 = padarray(diff(diff1, 1, 1), [1, 0]);
                diff2 = padarray(diff(diff2, 1, 2), [0, 1]);
                diff1(:,1) = 0; diff1(:,end) = 0;
                diff2(1,:) = 0; diff2(end,:) = 0;
                ob(1,:) = 0; ob(end,:) = 0; ob(:,1) = 0; ob(:,end) = 0;
                 
                grado{obnum,obmode,n} = grado{obnum,obmode,n} + 2*creg*(R/norm_r*ob-diff1-diff2);
                
                % -- Update only some roi by setting part of the gradient to 0
                if ~isempty(p.ms_grado_roi) 
                    apod = 100;
                    % mask_apod = fftshift(fract_hanning_pad(ob_dims(1), ob_dims(1)-asize+apod, ob_dims(1)-asize));
                    temp_grado = grado{obnum,obmode,n}(:,p.ms_grado_roi(1):p.ms_grado_roi(2));
                    mask_apod = fftshift(fract_hanning_pad(size(temp_grado,2), size(temp_grado,2)-apod*2+apod, size(temp_grado,2)-apod*2));
                    
                    temp_grado_apod = temp_grado .* mask_apod(round(size(temp_grado,2)/2),:);

                    grado{obnum,obmode,n} = zeros(size(grado{obnum,obmode,n}));
                    grado{obnum,obmode,n}(:,p.ms_grado_roi(1):p.ms_grado_roi(2)) = temp_grado_apod;

                end
            end
        end
    end
end

if ~optimize_delta_z
    core.errorplot(func);
end
iteration = length(core.errorplot([]));
verbose(2, 'Scan %s ; ML_iteration # %d of %d',errtitlestring, iteration, p.ms_opt_iter);

if optimize_delta_z
    verbose(2, ['Starting linesearch, delta_z currently ' num2str(delta_z*1e6,5) ' um, Error = ' num2str(func)]);
else
    verbose(2, ['Starting linesearch, delta_z fixed at ' num2str(delta_z*1e6,5) ' um, Error = ' num2str(func)]);
end

% if ~isempty(p.ms_grado_roi) 
%     if 0 %iteration==1 || mod(iteration,10)==0
%         save(sprintf('%smatlab/ptycho/temp_grado_it%d.mat', p.base_path, iteration),'temp_grado');
%         save(sprintf('%smatlab/ptycho/temp_grado_apod_it%d.mat', p.base_path, iteration),'temp_grado_apod');
%     end
% end

%%% Probe support constratint %%%
if p.use_probe_support && optimize_probes
    for prnum = 1:numprobs
        for prmode = 1:p.probe_modes
            gradp(:,:,prnum,prmode) = probe_mask.*gradp(:,:,prnum,prmode);
        end
    end
end

%%% Scaling preconditioning %%%
avobint = 0;
if p.scale_gradient && optimize_probes
    for ii = 1:numscans
        for n = 1:N_layer
            if p.share_probe
                
                grado_sum = 0;
                for obmode = 1:p.object_modes
                    grado_sum = grado_sum + sum( abs(grado{obnum,obmode,n}(:)).^2 );
                end
                
                avobint = avobint + grado_sum;
                if ii == numscans
                    avobint = avobint/numscans;
                    gradp = sqrt( avobint/sum( abs(gradp(:)).^2 ) )*gradp;
                end
                
            else
                grado_sum = 0;
                for obmode = 1:p.object_modes
                    grado_sum = grado_sum + sum( abs(grado{obnum,obmode,n}(:)).^2 );
                end
                
                gradp(:,:,ii,:) = sqrt(grado_sum / sum(sum(sum( abs(gradp(:,:,ii,:)).^2 )))) * gradp(:,:,ii,:);
            end
        end
    end
end

%%% Arranging gradients vector %%%
if ~any(p.ms_opt_flags)
    error('At least one element of flags must be 1');
end

grad = []; % Optimization vector

%figure(1000); clf
if optimize_object_layer
    for obnum = 1:numobjs
        for obmode = 1:p.object_modes
            for n = 1:N_layer
                fprintf('sum abs grado{%d,obmode,n%d} = %s \n', obnum, n, num2str(sum(abs(grado{obnum,obmode,n}(:)))) );
                grad = [grad; reshape([real(grado{obnum,obmode,n}(:)).'; imag(grado{obnum,obmode,n}(:)).'], [], 1)];
                
                %subplot(N_layer,1,n); 
                %hold on; plot(abs(grado{obnum,obmode,n}(:))); 
                %grid on; title(sprintf('Slice %d',n)); 
                %axis tight
                %drawnow 
            end
        end
    end
end

if optimize_probes
    fprintf('sum abs gradp = %s \n', num2str(sum(abs(gradp(:)))));
    grad = [grad; reshape([real(gradp(:)).'; imag(gradp(:)).'], [], 1)];
    
end

if optimize_delta_z
    fprintf('sum abs gradz = %s \n', num2str(sum(abs(gradz(:)))));
    grad = [grad; gradz(:)];
end
end

