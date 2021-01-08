% ASTRA_INITIALIZE Generate inputs needed for astra MEX wrapper
%
% [cfg, vectors] = ASTRA_initialize(Npix, size_projection,angles,lamino_angle, tilt_angle, pixel_scale, rotation_center)
%
% Inputs for angular geometry (!! all angles are expected in degress !!):
%   **Npix - size of tomogram   
%   **size_projection - size of sinogram (Nlayers, width, Nangles)
%   **angles - rotation angles of projections in degrees
% *optional*:
%   **lamino_angle - laminography angle / angles in degrees. lamino_angle ==
%                   90 is standard tomography , default = 90
%   **tilt_angle - tilt of camera with respect to the rotation axis coordinates, in degrees, default = 0
%   **pixel_scale - scale of pixels in tomogram compares to the
%               projection pixel size, default = 1
%   **rotation_center - center of rotation coordinates, default = size_projection/2
%   **skewness_angle - distorsion of parallel axis by [1, sind(alpha); 0, 1]
%
% Inputs for rotation matrix geometry: 
%   **Npix - size of tomogram   
%   **size_projection - size of sinogram (Nlayers, width, Nangles)
%   **rotation_matrix - R is a 3x3xn matrix, for n projections
%
% Outputs:
%   ++cfg  - config structure for ASTRA mex wrapper 
%   ++vectors - parameter vector for each angle 

   
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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





function [cfg, vectors] = ASTRA_initialize(Npix, size_projection, varargin) 
    use_rotmat = length(varargin) == 1 && size(varargin{1},1) == 3 && size(varargin{1},2) == 3 ; 
    if use_rotmat
        % rotation matrices were provided
        rot_mat = varargin{1}; 
        Nangles = size(rot_mat,3);
        r.rotation_center = size_projection/2;  % only centered geometry is supported when rotation matrix is provided
    else
        % angles and other parameters were provided 
        par = inputParser;
        par.KeepUnmatched = true; 
        %% ALL angles are assumed in degrees
        par.addRequired('angles')  % rotation angles of projections in degrees
        par.addOptional('lamino_angle',     90, @isnumeric)   % laminography angle / angles in degrees. lamino_angle == 90 is standard tomography , default = 90
        par.addOptional('tilt_angle',       0, @isnumeric)   % tilt of camera with respect to the rotation axis coordinates, in
        par.addOptional('pixel_scale',      [1,1], @isnumeric)   % scale of pixels in tomogram compares to the projection pixel size, default = 1
        par.addOptional('rotation_center', size_projection/2, @isnumeric)   % center of rotation cooridinates, default = size_projection/2        
        par.addOptional('skewness_angle', 0)   % distorsion of parallel axis by [1, sind(alpha); 0, 1]
        par.addParameter('show_geometry', false)   % plot also a geometry for each projection 
        par.parse(varargin{:})
        r = par.Results;
        Nangles = length(r.angles);
    end
    
    
    % angles should be sorted in order to maximize performance of the astra
    % toolbox -> better use of texture memory
   
    assert(math.isint(Npix), 'Npix is not integer');
    assert(math.isint(size_projection), 'size_projection is not integer');
       
    if isscalar(Npix)
        Npix(2) = Npix;
    end
    if length(Npix) == 2 && all(r.lamino_angle == 90)
        Npix(3) = size_projection(1);  % default behaviour is to have same number of layers in reconstruction and in laminography
    elseif length(Npix) == 2 && any(r.lamino_angle ~= 90)
        error('All three dimensions of the volume size has to be specified for the laminograhy geometry')
    end
   
    cfg.iVolX = Npix(1);
    cfg.iVolY = Npix(2);
    cfg.iVolZ = Npix(3);
    cfg.iProjAngles = Nangles;
    cfg.iProjU = size_projection(2);
    cfg.iProjV = size_projection(1);
    cfg.iRaysPerDet = 1;
    cfg.iRaysPerDetDim = 1;
    cfg.iRaysPerVoxelDim = 1;
    source_distance = 1; % currenlty not implemented in the ASTRA wrapper 


    if use_rotmat
        [vectors] = astra_convert_R_vectors(rot_mat,Nangles);
    else
        % compatibility with iradonfast
        r.angles = r.angles + 90;
        cfg.lamino_angle = r.lamino_angle;
        cfg.pixel_scale = r.pixel_scale;
        cfg.tilt_angle = r.tilt_angle;
        cfg.skewness_angle = r.skewness_angle;
        vectors = ASTRA_get_geometry(r.angles, r.lamino_angle, r.tilt_angle, source_distance,r.pixel_scale,r.skewness_angle,r.show_geometry);
    end

    %%%% apply geometry correction to shift reconstruction into center %%%%%%%%%%%%%%%%%%%%%%%%%%
    vectors(:,4:6) = vectors(:,4:6) -(vectors(:,10:12).*(r.rotation_center(:,1) )+vectors(:,7:9).*(r.rotation_center(:,2) ));

end

function vectors = ASTRA_get_geometry(angles, lamino_angle, tilt_angle, source_distance,pixel_scale,skewness_angle,show)

    Nangles = numel(angles);
    % angles should be sorted in order to maximize performance of the astra
    % toolbox -> better use of texture memory
   
    angles = deg2rad(angles(:));
    lamino_angle = pi/2 - deg2rad(lamino_angle);
    tilt_angle  = deg2rad(tilt_angle);
    skewness_angle  = deg2rad(skewness_angle);

    if isscalar(lamino_angle)
        lamino_angle = lamino_angle .* ones(Nangles,1);
    end
    if isscalar(tilt_angle)
        tilt_angle = tilt_angle .* ones(Nangles,1);
    end
    if isscalar(skewness_angle)
        skewness_angle = skewness_angle .* ones(Nangles,1);
    end  
    if isscalar(pixel_scale) || numel(pixel_scale) == 2
        pixel_scale = bsxfun(@times, pixel_scale , ones(Nangles,2));
    end
       
      % We generate the same geometry as the circular one above.
      vectors = zeros(Nangles, 12);
      % ray direction
      vectors(:,1) = sin(angles).*cos(lamino_angle);
      vectors(:,2) = -cos(angles).*cos(lamino_angle);
      vectors(:,3) = sin(lamino_angle);

      vectors(:,1:3) =  vectors(:,1:3) .*source_distance;
      % center of detector
      vectors(:,4:6) = 0;
      % vector from detector pixel (0,0) to (0,1) 
      vectors(:,7) = cos(angles)./pixel_scale(:,1);
      vectors(:,8) = sin(angles)./pixel_scale(:,1);
      vectors(:,9) = 0/pixel_scale(:,1);

      % vector from detector pixel (0,0) to (1,0) 

      % cross(vectors(i,1:3), vectors(i,7:9))
      % dot(vectors(i,1:3), vectors(i,7:9))

      vectors(:,10) = - sin(lamino_angle).*sin(angles)./pixel_scale(:,2);
      vectors(:,11) =   sin(lamino_angle).*cos(angles)./pixel_scale(:,2);
      vectors(:,12) =   cos(lamino_angle)./pixel_scale(:,2);    

    %  Rodrigues' rotation formula - rotate detector in plane
    %  perpendicular to the beam axis
      if any(tilt_angle ~= 0)
          for i = 1:Nangles
              vectors(i,7:9)=vectors(i,7:9).*cos(tilt_angle(i)) + ...
                             cross(vectors(i,1:3), vectors(i,7:9)).*sin(tilt_angle(i)) + ...
                            (vectors(i,1:3)*dot(vectors(i,1:3),vectors(i,7:9))).*(1-cos(tilt_angle(i)));
              vectors(i,10:12)=vectors(i,10:12).*cos(tilt_angle(i)) + ...
                             cross(vectors(i,1:3), vectors(i,10:12)).*sin(tilt_angle(i)) + ...
                            (vectors(i,1:3).*dot(vectors(i,1:3),vectors(i,10:12))).*(1-cos(tilt_angle(i)));
          end
      end
      
      % search also for skewness => the same as rotation, but rotate
      % only one axis of the detector !!
      if any(skewness_angle ~= 0)
          for i = 1:Nangles
                vectors(i,10:12)=vectors(i,10:12).*cos(skewness_angle(i)/2) + ...
                         cross(vectors(i,1:3), vectors(i,10:12)).*sin(skewness_angle(i)/2) + ...
                        (vectors(i,1:3).*dot(vectors(i,1:3),vectors(i,10:12))).*(1-cos(skewness_angle(i)/2));
          end
      end
         
      %% PLOT THE CURRENT SETUP

      if show
        for i = 1:Nangles
            draw_projection_geometry(vectors(i,:))
        end
      end


   

end

function [vectors] = astra_convert_R_vectors(rot_mat,Nangles)

% R is defined as in the arbitrary projection code by Manuel
% This integrates along z (2nd index) and so this code follows this
% convention.
pixel_scale(1:Nangles,1:2) = [1];


convert_matrix =   [0 0 1
                    0 1 0
                    1 0 0];
    
for ii=1:size(rot_mat,3)
    rot_mat(:,:,ii)=rot_mat(:,:,ii)*convert_matrix;
end

vectors = zeros(Nangles, 12);
    for i = 1:Nangles
        % before starting to mess around: works with a correction matrix
        % with the magnetic contrast, but not for the laminography
          vectors(i,1) = -rot_mat(3,1,i);
          vectors(i,2) = -rot_mat(3,3,i);
          vectors(i,3) = -rot_mat(3,2,i);
                   
          vectors(i,4:6) = 0;
         
          vectors(i,7) = rot_mat(1,1,i)/pixel_scale(i,1);
          vectors(i,8) = rot_mat(1,3,i)/pixel_scale(i,1);
          vectors(i,9) = rot_mat(1,2,i)/pixel_scale(i,1);

          vectors(i,10) = rot_mat(2,1,i)/pixel_scale(i,2);
          vectors(i,11) = rot_mat(2,3,i)/pixel_scale(i,2);
          vectors(i,12) = rot_mat(2,2,i)/pixel_scale(i,2);
    end
end


function draw_projection_geometry(vectors)
    % show geometry saved in the "vectors" matrix
    ray = vectors(1:3);
    c_center = vectors(4:6)-ray;
    c_origin = c_center - vectors( 7:9)/2-vectors( 10:12)/2;

    k = 6;
    n = 2^k-1;
    [x,y,z] = sphere(n);
    c = hadamard(2^k);
    s = 0.5;
    figure(15)
    surf(vectors(4)+s*x,vectors(5)+s*y,vectors(6)+s*z,c);
    shading flat
    colormap([1  1  0; 0  1  1])

    hold all
    plot3d_vec(c_origin, vectors( 7:9), 'r');
    plot3d_vec(c_origin, vectors( 10:12), 'r');
    plot3d_vec(c_origin+vectors( 7:9), vectors( 10:12), 'r');
    plot3d_vec(c_origin+vectors( 10:12), vectors( 7:9), 'r');
    % draw "pixels" on a 10x10 grid
    for x = linspace(0,1,10)
        plot3d_vec(c_origin+x*vectors( 10:12), vectors( 7:9), 'r:');
        plot3d_vec(c_origin+x*vectors( 7:9),vectors( 10:12), 'r:');
    end
    plotting.mArrow3(c_center+ray*2,c_center,  'color', 'blue', 'stemWidth',0.02,'facealpha',0.5);

    hold off
    axis([-1,1,-1,1,-1,1])
    drawnow 
   
end

function h_out = plot3d_vec(x0, vec, varargin)
      h = plot3(x0(1)+[0,vec(1)], ...
          x0(2)+[0,vec(2)], ...
          x0(3)+[0,vec(3)], varargin{:});
      if nargout > 1
          h_out = h;
      end
      
end

