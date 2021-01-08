%LOAD_DATA prepare filenames and load data

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

function [ p ] = load_data( p )
    import math.*
    import utils.*
    
    detStorage = p.detectors(p.scanID).detStorage;


    %% VIRTUAL LOADING FUNCTION THAT GENERATES ARTIFICIAL DATA 
    
    
    %% adjust positions , apply affine corrections from template and distorsion from p.simulation.affine_matrix
    if isempty(p.affine_matrix)
       p.affine_matrix = [1,0;0,1]; 
    end

    if ~isfield(p.simulation, 'affine_matrix') || isempty(p.simulation.affine_matrix)
       p.simulation.affine_matrix = [1,0;0,1]; 
    end

    %% calculate the positions 
    positions_0 = p.positions; % store the p.positions and return the values at the end of this function 
    tmp = p;
    tmp.affine_matrix = -p.simulation.affine_matrix * inv(p.affine_matrix); % first remove the already applied affine matrix and then apply affine matrix from simulation
    
    for ii = 1:p.numscans
        % add there a small random global offset for the positions 
        rng(p.scan_number(ii))  % reset randomization to guarantee repeatability 
        offset = 0.3; % times average step 
        avg_step = sqrt(prod(max(tmp.positions_real(p.scanidxs{ii},:)) - min(tmp.positions_real(p.scanidxs{ii},:))) / p.numpts(ii)); 
        tmp.positions_real(p.scanidxs{ii},:) = tmp.positions_real(p.scanidxs{ii},:) + randn(1,2) * avg_step * offset; 
    end

    % standard farfield ptychography
    if check_option(p.simulation,  'z') && strcmpi(p.prop_regime, 'farfield')
        tmp.dx_spec = tmp.lambda*tmp.simulation.z ./ (p.asize*p.ds);                   % resolution in the specimen plane
        tmp.dx_spec = tmp.dx_spec ./ cosd(tmp.sample_rotation_angles(1:2));      % account for a tilted sample ptychography      
    end
    
    tmp.share_object_ID = ones(p.numscans,1); 
    tmp = core.ptycho_adjust_positions( tmp ); 

    
    if p.simulation.position_uncertainty > 0 
        tmp.positions = tmp.positions + randn(sum(p.numpts),2)*p.simulation.position_uncertainty*mean(p.dx_spec);
        utils.verbose(3, 'Included position errors, std=%3.2gnm , %3.2gpx', p.simulation.position_uncertainty*1e9, p.simulation.position_uncertainty / mean(p.dx_spec) ); 
    end
    
    % Extra offset of positions given in simulation -> subtract padding
    % provided in p structure and instead add padding from simulation structure 

    p.simulation.positions = tmp.positions;
    p.simulation.positions_real = tmp.positions_real;
    
    %% create object
    if ~check_option(p.simulation, 'positions_pad')
        p.simulation.positions_pad = [0,0]; 
    end
    % Compute object sizes
    if p.share_object
        p.object_size = ceil(p.asize + max(p.simulation.positions) + p.positions_pad+ p.simulation.positions_pad(1,:));
    else
        for ii = 1:p.numscans
            p.object_size(ii,:) = ceil(p.asize + max(p.simulation.positions(p.scanidxs{ii},:)) + p.positions_pad + p.simulation.positions_pad(min(end,ii),:));
        end
    end
    

    % Generate object using parameters from artificial data template 
    [p.simulation.obj, p.simulation.ref_index] = detector.virtual.create_object(p);

    

    if get_option(p, 'fourier_ptycho')
        warning('FIXME')
        keyboard
        p = fourier_ptycho_data(p);
    end

    
    p.simulation.probe = p.probes; 

    if ~p.simulation.apply_sub_px_shifts
        p.positions = round(p.positions); 
        p.simulation.positions = round(p.simulation.positions);
    end
    
    p.positions = p.simulation.positions;

    if p.simulation.sample_rotation_angles
        % modify the scanning positions to account for the tilted sample geometry 
        p.positions = p.positions .* cosd(p.simulation.sample_rotation_angles([1,2]));
    end
    
    sub_px_shift = p.positions-round(p.positions); 

    if any(p.simulation.sample_rotation_angles)
        % get propagators to the tilted plane or plane rotated around beam axis 
        [fwd_propag_fun, back_propag_fun] = get_tilted_plane_propagators(p.probes, p.simulation.sample_rotation_angles,...
                                        p.lambda, p.dx_spec); 
        
    end

    %% calculate views
    Nlayers = size(p.simulation.obj{p.scanID},4);
    if p.share_object
        obnum = 1;
    else
        obnum = p.scanID;
    end

    verbose(0, 'Creating artificial dataset')
    % auxiliar windows for subpixel shifting 
    win = 0.1+0.9*tukeywin(p.asize(1),0.05) .* tukeywin(p.asize(2), 0.05)';  

    iter = zeros([p.asize, length(p.scanidxs{obnum}), p.probe_modes*p.object_modes], 'like', p.simulation.obj{p.scanID});
    for prmode = 1:p.probe_modes
        for obmode = 1:p.object_modes
            p.simulation.obj{obnum} = single(p.simulation.obj{obnum});
            iter_mode_ind = prmode+(obmode-1)*p.probe_modes;

            probe = p.probes(:,:,1,prmode);
            if check_option(p, 'use_gpu')
                probe = utils.Garray(probe); 
            end
            
            if strcmpi(p.prop_regime, 'farfield')
                % in farfield is the probe and object shift equivalent in nearfield not anymore !!
                probe = imshift_fft(probe,sub_px_shift(p.scanidxs{p.scanID},[2,1]));                
            end
            if Nlayers > 1 && p.simulation.thickness > 0
                % assume that the provided probe is in center plane of the sample
                probe = prop_free_nf(probe, p.lambda, -p.simulation.thickness/2, p.dx_spec(1)); 
            end
            if any(p.simulation.sample_rotation_angles(1:2))
                % propagate the probe to the tilted plane 
                probe = fwd_propag_fun(probe);
            end

            if  p.simulation.thickness == 0
                % thin object 
                obj = prod(p.simulation.obj{obnum}(:,:,obmode,:),4); 
            else
                obj = p.simulation.obj{obnum}(:,:,obmode,:); 
            end

    
            obj_proj = core.get_projections(p, obj(:,:,obmode,1) , p.scanID);
            proj = probe .* obj_proj; 
           

            if  p.simulation.thickness > 0 
                % thick object 
                [~,H] = prop_free_nf(probe, p.lambda, p.simulation.thickness / (Nlayers-1), p.dx_spec); 

                for layer = 2:Nlayers
                    if Nlayers > 2 && utils.verbose >= 0; utils.progressbar(layer-1, Nlayers-1); end
                    proj = ifft2(H.*fft2(proj));
                    obj_proj = core.get_projections(p, p.simulation.obj{obnum}(:,:,obmode,layer), p.scanID, obj_proj);
                    if strcmpi(p.prop_regime, 'nearfield')
                        % in farfield is the probe and object shift equivalent in nearfield not anymore !!
                        obj_proj = imshift_fft(obj_proj .* win,-sub_px_shift(p.scanidxs{p.scanID},[2,1])) ./ win;                
                    end
                    proj = proj .* obj_proj; 
                end
            else
                proj = bsxfun(@times, probe, obj_proj); 
            end
            iter(:,:,:,iter_mode_ind) = proj;
        end
    end
    
    if any(p.simulation.sample_rotation_angles)
        % perform propagation back to the plane parallel with detector 
        iter = back_propag_fun(iter);
    end

    %% create data 
    if check_option(p, 'prop_regime', 'nearfield')
        diffraction = abs(prop_free_nf(iter, p.lambda, p.z, p.dx_spec)).^2;
    elseif p.simulation.prop_from_focus 
        % calculate intensities
        diffraction = abs(prop_free_nf(ifftshift_2D(fft2(fftshift_2D(iter))), p.lambda, -p.simulation.prop_from_focus, p.ds)).^2;     
    else% standard farfield propagation
        diffraction = fftshift_2D(abs(fft2(iter)).^2);
    end
    diffraction = sum(diffraction,4); 
        
    
    %% add incoherence-like blur 
    if p.simulation.incoherence_blur
        diffraction = utils.imgaussfilt3_conv(diffraction, [p.simulation.incoherence_blur,p.simulation.incoherence_blur,0]);
        verbose(3,'- Adding incoherence blur %3.2gpx', p.simulation.incoherence_blur);
    end

    %% add noise 
    if ~isinf(p.simulation.photons_per_pixel)
        illum_sum = sum(diffraction(:)); 
        
        % calculate the total number of photons per scan 
        total_dose = p.simulation.photons_per_pixel *  prod(p.object_size-p.asize); 
        % calculate correction of the intensity 
        corr_ratio = sum(total_dose) / illum_sum; 
        % set identical number of photons to diffr. pattern
        diffraction = diffraction * corr_ratio; 
        p.simulation.probe = p.simulation.probe * sqrt(mean(corr_ratio)); 

        verbose(3,'- Adding Poisson noise');
        % add noise, it is some faster approximation 
        diffraction = randpoisson(diffraction);
    else
        % keep the values close to real X-ray data
        max_value = 1e3; 
        corr_ratio = max_value / max(diffraction(:));
        diffraction = diffraction * corr_ratio; 
        p.simulation.probe = p.simulation.probe * sqrt(corr_ratio); 
        verbose(3,'- No noise added');
    end

    if check_option(p, 'use_gpu')
        diffraction = gather(diffraction); 
        p.simulation.probe = gather(p.simulation.probe); 
        for ii = 1:p.numscans
            p.simulation.obj{ii} = gather(p.simulation.obj{ii}); 
        end
    end
    
    
        
    detStorage.data = double(diffraction);
    detStorage.mask = true(size(diffraction));
    p.positions = positions_0; 
        

end
