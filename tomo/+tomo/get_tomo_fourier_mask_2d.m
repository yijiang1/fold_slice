function [mask] = get_tomo_fourier_mask_2d( Nproj , N_recon, angles )
%Get fourier mask of the MISSING WEDGE

%Output Image Size
Ny = N_recon(1);
Nx = N_recon(2);

cen_y = floor(Ny/2)+1;
cen_x = floor(Nx/2)+1;

dk_y = 1/Ny;
dk_x = 1/Nx;
Nang = length(angles);
x = (linspace(1,Nproj,Nproj)-cen_x)*dk_x;
x = repmat(x,[Nang,1]);
y = 0;
angles_temp = reshape(angles,[Nang,1]);
angles_temp = repmat(angles_temp,[1,Nproj]);
angles_temp = -angles_temp*pi/180;
y_new = cos(angles_temp)*y + sin(angles_temp).*x;
x_new = -sin(angles_temp)*y + cos(angles_temp).*x;         

mask = zeros(Ny, Nx, 4);

%1
p_y = floor(y_new/dk_y)+cen_y;
p_x = floor(x_new/dk_x)+cen_x;
mask_temp = zeros([Ny*Nx,1]);
mask_temp(p_y+(p_x-1)*Ny)=1;
mask(:,:,1) = reshape(mask_temp,[Ny,Nx]);

%2
p_y = ceil(y_new/dk_y)+cen_y;
p_x = floor(x_new/dk_x)+cen_x;
mask_temp = zeros([Ny*Nx,1]);
mask_temp(p_y+(p_x-1)*Ny)=1;
mask(:,:,2) = reshape(mask_temp,[Ny,Nx]);

%3
p_y = floor(y_new/dk_y)+cen_y;
p_x = ceil(x_new/dk_x)+cen_x;
mask_temp = zeros([Ny*Nx,1]);
mask_temp(p_y+(p_x-1)*Ny)=1;
mask(:,:,3) = reshape(mask_temp,[Ny,Nx]);

%4
p_y = ceil(y_new/dk_y)+cen_y;
p_x = ceil(x_new/dk_x)+cen_x;
mask_temp = zeros([Ny*Nx,1]);
mask_temp(p_y+(p_x-1)*Ny)=1;
mask(:,:,4) = reshape(mask_temp,[Ny,Nx]);

mask = sum(mask,3);
mask = mask==0;
% always set center to 1
%mask(cen_y,cen_x) = 1;
        
%{
%Number of projections
N_ang = max(size(angles));

%Output Image Size
Ny = N_recon(1);
Nx = N_recon(2);

cen_y = floor(Ny/2)+1;
cen_x = floor(Nx/2)+1;
    
mask = zeros(Ny, Nx);


dk_y = 1/Ny;
dk_x = 1/Nx;

for a = 1:N_ang
    %if angles(a)~=0 || angles(a)~=90
    ang = -angles(a)*pi/180;

    %CurrentcN_recon projection
    %p = input(:,a);
    %P = fftshift(fft(ifftshift(p)));
    %{
    if -angles(a) == 90        
        v(:,cen) = v(:,cen) + P(end:-1:1);

        w(:,cen) = w(:,cen) + 1;

        output_F(w~=0) = v(w~=0)./w(w~=0);
    elseif -angles(a) == -90     
        v(:,cen) = v(:,cen) + P;
        w(:,cen) = w(:,cen) + 1;

        output_F(w~=0) = v(w~=0)./w(w~=0);

    else
    %}
   
    for i=1:1:N_proj
        x = (i-cen_x)*dk_x;
        y = 0;

        %Calculate new coordinate 
        y_new = cos(ang)*y + sin(ang)*x;
        x_new = -sin(ang)*y + cos(ang)*x; 

        %Calculate weights
        sy = abs(floor(y_new) - y_new); if sy<1e-6; sy=0; end
        sx = abs(floor(x_new) - x_new); if sx<1e-6; sx=0; end

        %Bilinear Extrapolation
        %P1  
        p_y = floor(y_new/dk_y)+cen_y;
        p_x = floor(x_new/dk_x)+cen_x;
        if p_x >0 && p_x<=Nx && p_y >0 && p_y<=Ny
            mask(p_y,p_x) = 1; 
        end

        %P2
        p_y = ceil(y_new/dk_y)+cen_y;
        p_x = floor(x_new/dk_x)+cen_x;
        if p_x>0 && p_x<=Nx && p_y >0 && p_y<=Ny
            mask(p_y,p_x) = 1; 
        end

        %P3
        p_y = floor(y_new/dk_y)+cen_y;
        p_x = ceil(x_new/dk_x)+cen_x;        
        if p_x>0 && p_x<=Nx && p_y >0 && p_y<=Ny
            mask(p_y,p_x) = 1; 
        end

        %P4
        p_y = ceil(y_new/dk_y)+cen_y;
        p_x = ceil(x_new/dk_x)+cen_x;
        if p_x>0 && p_x<=Nx && p_y >0 && p_y<=Ny
            mask(p_y,p_x) = 1; 
        end
    end

   
end
mask = mask==1;
%}

end

