function [] = adjust_cbed_center(cbed,delta_r,delta_k)
%Declare Global Variables
global cbed_sum cbed_sum_shifted S sx sy ny nx delta_s radius delta_radius;

cbed_sum = sum(sum(cbed,3),4);
cbed_sum_shifted = cbed_sum;
delta_s = delta_k;
delta_radius = delta_r;

Ny = size(cbed_sum,1);
Nx = size(cbed_sum,2);
x = linspace(-floor(Nx/2),ceil(Nx/2)-1,Nx);
y = linspace(-floor(Ny/2),ceil(Ny/2)-1,Ny);
[Y, X] = meshgrid(x,y);
S = sqrt(X.^2+Y.^2);
radius = delta_r;
% Create a figure and specify a callback
%figure('Name', 'Shift Reconstructions Fourier transform','KeyPressFcn',@keypressCallback, 'CloseRequestFcn',@closeCallback)
figure('Name', 'Shift Reconstructions Fourier transform','KeyPressFcn',@keypressCallback)
%Ensure pixels are proper ratio
[nx,ny, ~] = size(cbed_sum);
   
sx = 0;
sy = 0;
%O_recon_f = fftshift(fft2(ifftshift(recon_shifted)));

colormap('jet');
imagesc( abs(cbed_sum.*(S<=radius)));
set(gca,'plotboxaspectratio',[ny nx 1]);
%title( ['sx= ' num2str(sx) ', sy= ' num2str(sy)] );  
colorbar;

%{
IsComplete = 0;

while( ~IsComplete )
    pause();
end
%}
%output = recon_shifted;
end

% Callback subfunction for when a key is pressed
function keypressCallback(src,eventdata)

   %Define global variables
   global cbed_sum cbed_sum_shifted S sx sy ny nx delta_s radius delta_radius;
   if strcmp(eventdata.Key,'w')
       radius = radius + delta_radius;
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar;
       %output = recon_shifted;
   elseif strcmp(eventdata.Key,'s')
       radius = radius - delta_radius;
       if radius<0
            radius=0;
       end
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar;
       %output = recon_shifted;    
   elseif strcmp(eventdata.Key,'rightarrow')
       sx = sx + delta_s;
       cbed_sum_shifted = shift( cbed_sum,1,1,sx,sy ); 
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar;
       %output = recon_shifted;

   elseif strcmp(eventdata.Key,'leftarrow')
       sx = sx - delta_s;
       cbed_sum_shifted = shift( cbed_sum,1,1,sx,sy ); 
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar;
   elseif strcmp(eventdata.Key,'uparrow')
       sy = sy + delta_s;
       cbed_sum_shifted = shift( cbed_sum,1,1,sx,sy ); 
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar; 
   elseif strcmp(eventdata.Key,'downarrow')
       sy = sy - delta_s;
       cbed_sum_shifted = shift( cbed_sum,1,1,sx,sy ); 
       imagesc( abs(cbed_sum_shifted.*(S<=radius)));
       set(gca,'plotboxaspectratio',[ny nx 1]);
       title( ['sx= ' num2str(sx) ', sy= ' num2str(sy) ', radius= ' num2str(radius)] );
       colorbar; 
   %elseif eventdata.Character == 'q'
   %     IsComplete = 1;
   else
       
      fprintf('Null\n');
   end
    
end

%{
function closeCallback(src,eventdata)
    %global IsComplete
    %IsComplete = 1;

    delete(src);

end
%}