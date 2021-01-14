%%% Following a feature
for ii = 1:numel(object)
    obj{ii} = angle(object{ii});
end
pixsize = p.dx_spec(1)*1e6; % Microns

%Inputs
% obj       Is a cell with different objects
% pixsize   Pixel size

% Create useful arrays
for ii = 1:numel(object)
    axisx{ii} = ([1:size(obj{ii},2)]-floor(object_size(2)/2)+1)*pixsize;
    axisy{ii} = ([1:size(obj{ii},1)]-floor(object_size(1)/2)+1)*pixsize;
    xind {ii} = [1:size(obj{ii},2)];
    yind{ii}  = [1:size(obj{ii},1)];  
end


%% Feature characteristics
f.sigma = 1.5; % Feature width in real units
f.contrast = -1;

sigma_ind = f.sigma/pixsize;  % Feature width in pixels

figure(1)
% imagesc(axisx{1},axisy{1},obj{1});
imagesc(obj{1});
axis xy equal tight
colormap bone
xlabel('\mum')
[xinp,yinp] = ginput(1);
xinp = round(xinp);
yinp = round(yinp);

x1 = xind{1}(abs(xind{1}-xinp)<2*sigma_ind);
y1 = yind{1}(abs(yind{1}-yinp)<2*sigma_ind);
x2 = x1;
y2 = y1;

for ii = 1:3
    [X Y] = meshgrid(xind{ii},yind{ii});
    ref = f.contrast*exp(-((X-xinp).^2+(Y-yinp).^2)/(2*sigma_ind^2));
    x1 = xind{1}(abs(xind{1}-xinp)<2*sigma_ind);
    y1 = yind{1}(abs(yind{1}-yinp)<2*sigma_ind);
    x2 = x1;
    y2 = y1;
    [subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(obj{1}, ref, x1, y1, x2, y2, 10, 1, 1);
    % delta is (y,x) correction
    xinp = xinp - delta(2);
    yinp = yinp - delta(1);
    xposobjind(ii) = xinp;
    yposobjind(ii) = yinp;
end

%%
ii = 3;
figure(2)
% imagesc(axisx{1},axisy{1},ref);
% imagesc(ref);
imagesc(obj{ii});
axis xy equal tight
colormap bone
hold on;
plot(xposobjind(ii),yposobjind(ii),'ow')
hold off;
xlabel('pixels')