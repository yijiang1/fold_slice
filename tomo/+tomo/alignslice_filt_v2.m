% [deltastack outputsinogram] =  alignslice_filt(sinogramder,theta,deltastack,pixtol,disp,paramsalign)
% Function to align projections. It relies on having already aligned the
% vertical direction. The code aligns using the consistency before and
% after tomographic combination of projections. Currently it can only deal
% with one slice. It is recommended to use small values of tomographic
% filter for best results.
% WARNING. The code aligns the center of rotation in the ROI defined by
% limsx so that it is consistent (within the window) with the center
% expected by iradon: ceil(size(R,1)/2). Also note that it expects the
% derivative of the sinogram as input
%
%
% Inputs %
% sinogramder   Sinogram derivative, the second index should be the angle
% deltastack    Row array with initial estimates of positions (1,n)
% disp          Display = 0 no images                
%               = 1 Final diagnostic images
%               = 2 Diagnostic images per iteration  (significant overhead)
% pixtol        Tolerance for alignment, it is also used as a search step
%
% Extra parameters (in structure paramsalign)
% paramsalign.interpmeth    = 'sinc' - Shift images with sinc interpolation (default)
%                           = 'linear' - Shift images with linear interpolation (default)
% paramsalign.usecircle     Use a circular mask to eliminate corners of
%                           tomogram
% paramsalign.filtertomo    Frequency cutoff for tomography filter
% paramsalign.cliplow       Minimum value in tomogram
% paramsalign.cliphigh      Maximum value in tomogram
% paramsalign.masklims      Mask in sinograms to evaluate error metric
% paramsalign.binning       Binning of sinograms to increase speed. I am
%                            not fully aware
%
% Outputs %
% deltastack        Object positions
% outputsinogram    Aligned sinogram derivatives (optional)
%
% Manuel Guizar 13 July 2011
% This code is provided as is and without guarantees on its performance
% The algorithm is not yet published. Do not distribute and contact 
% (mguizar@gmail.com) prior to publication of results where this code plays 
% a significant role, or to report a bug.
% Please acknowledge if used.
% v2 Manuel Guizar 2011 09 29
% Added clipping values to tomogram
% Manuel Guizar 2016 03 01 - Adding binning possiblity. Currenlty I suspect
% there may be an offset introduced between different binnings. Beware of
% this while using it and report any behavior that points in this direction.

function [deltastack outputsinogram] =  alignslice_filt_v2(sinogramder,theta,deltastack,pixtol,disp,paramsalign)

import utils.*
import math.projectleg1D_2

%deltastack = round(deltastack);
maxit = 40;

display(['   Registration of sinogram - Single pixel'])
if nargout>1,
    outputsinogram = sinogramder;
end

if exist('paramsalign') == 0,
    paramsalign = 0;
end
%%%%%%%
if isfield(paramsalign,'binning') == 0,
    binning = 0;
else
    binning = paramsalign.binning;
    display(['Using binning = ' num2str(binning)])
end

if (binning ~=0)&&(binning ~=1)
    display(['Using binning on sinogram = ' num2str(binning)]);
    sinogramder_orig = sinogramder;
    sinoaux = sinogramder(1:binning:end-binning+1,:);
    for ii = 2:binning
        sinoaux = sinoaux + sinogramder(ii:binning:end-binning+ii,:);
    end
    sinogramder = sinoaux/binning;
    deltastack = deltastack/binning;
end
%%%%%%

if isfield(paramsalign,'masklims') == 0,
    masklims = [1:size(sinogramder,1)];
else
    display('Error computed using masked values of sinogram')
    masklims = paramsalign.masklims;
end

if isfield(paramsalign,'interpmeth') == 0,
    paramsalign.interpmeth = 'sinc';
    interpmeth = paramsalign.interpmeth;
else 
    interpmeth = paramsalign.interpmeth;
    if (strcmp(interpmeth,'sinc'))||(strcmp(interpmeth,'linear'))
    else
        error('Undefined interpolation method')
    end
end

if isfield(paramsalign,'usecircle') == 0,
    usecircle = true;
else
    usecircle = paramsalign.usecircle;
end

% if isfield(paramsalign,'expshift') == 0,
%     expshift = false;
% else
%     expshift = paramsalign.expshift;
% end

if isfield(paramsalign,'filtertomo') == 0,
    filtertomo = 1;
    display('Using default filter cutoff = 1')
else
    filtertomo = paramsalign.filtertomo;
    display(['Using filtertomo = ' num2str(filtertomo)])
end

if isfield(paramsalign,'cliplow') == 0,
    cliplow = [];
else
    cliplow = paramsalign.cliplow;
    display(['Low limit for tomo values = ' num2str(cliplow)])
    
end

if isfield(paramsalign,'cliphigh') == 0,
    cliphigh = [];
else
    cliphigh = paramsalign.cliphigh;
    display(['High limit for tomo values = ' num2str(cliphigh)])
    warning('Code will use a background "cliphigh", assumes tomogram will go negative from there')
end

if exist('pixtol') == 0,
    pixtol = 1;
end

if exist('deltastack') == 0,
    deltastack = zeros(1,size(sinogramder,2));
end

% Pad sinogram derivatives
padval = 2*round(1/filtertomo);
sinogramder = padarray(sinogramder,[padval 0]);

%%%%%%%%%%%%% Single pixel precision %%%%%%%%%%%%%%
% loop here
clear auxinit;  % Before main loop
domain = 1;
count = 0;
N = size(sinogramder,1);
center = floor((N+1)/2);
xt = [-N/2:N/2-1];
[Xt Yt] = meshgrid(xt,xt);
circulo = 1-radtap(Xt,Yt,10,N/2-10);

filteraux = 1-fract_hanning_pad(N,N,round(N*(1-filtertomo)));
filteraux = repmat(fftshift(filteraux(:,1)),[1 length(theta)]);

sinogramderorig_nofilt = sinogramder;
sinogramderorig = real(ifft(fft(sinogramder).*filteraux));
erroryreg = inf;
while domain == 1,
    count = count+1;
    deltaprev = deltastack;
    erroryregprev = erroryreg;
    for ii = 1:size(sinogramder,2)
        if strcmp(interpmeth,'sinc')
            sinogramder(:,ii) = real(shiftpp2(sinogramderorig(:,ii),deltastack(1,ii),0));
        elseif strcmp(interpmeth,'linear')
            sinogramder(:,ii) = shiftwrapbilinear(sinogramderorig(:,ii),deltastack(1,ii),0);
        end
        if disp>1,
            if strcmp(interpmeth,'sinc')
                sinogramder_nofilt(:,ii) = real(shiftpp2(sinogramderorig_nofilt(:,ii),deltastack(1,ii),0));
            elseif strcmp(interpmeth,'linear')
                sinogramder_nofilt(:,ii) = shiftwrapbilinear(sinogramderorig_nofilt(:,ii),deltastack(1,ii),0);
            end
        end
    end

    %%%%%% Sinogram alignment single pixel precision
    % Compute tomogram with current sinogramder
%     recons = iradonfast_v3(double(sinogramder),...
%         theta,'linear','derivative','Han',size(sinogramder,1),filtertomo);
    recons = tomo.iradonfast_v3(double(sinogramder),...
        theta,'linear','derivative','Ram-Lak',size(sinogramder,1),1);
    if ~isempty(cliplow)
        recons = recons.*(recons>=cliplow) + cliplow*(recons<cliplow);
    end
    if ~isempty(cliphigh)
        recons = recons.*(recons<=cliphigh) + cliphigh*(recons>cliphigh);
        recons = recons-cliphigh;
    end
    if usecircle
        recons = recons.*circulo;
    end
    
    figure(203)
    imagesc(recons);
    axis xy equal tight
    colorbar
    
    % Computed sinogram
    sinogramcomp  = radon(recons,theta);
    Nbig = size(sinogramcomp,1);
    centerbig = floor((Nbig+1)/2);
    
    if strcmp(interpmeth,'sinc')
        sinogramcompder = real(shiftpp2(sinogramcomp,0.5,0)-shiftpp2(sinogramcomp,-0.5,0));
    elseif strcmp(interpmeth,'linear')
        sinogramcompder = real(shiftwrapbilinear(sinogramcomp,0.5,0)-shiftwrapbilinear(sinogramcomp,-0.5,0));
    end
    sinogramcompder = sinogramcompder([1:N]+centerbig-center,:);

    % Compare sinogramder with sinogramcompder
    for ii = 1:size(sinogramder,2)
        %ii,
        errorinit(ii) = sum(abs(sinogramder(padval+masklims,ii)-sinogramcompder(padval+masklims,ii)).^2);
    end
    display(['Initial error metric, E = ' num2str(sum(errorinit))])
    
    
%%% Search for shifts with respect to sinthetic
    for ii = 1:size(sinogramder,2);

        shift = 0;
        %%% Looking bothways
        % compute current shift error
        if strcmp(interpmeth,'sinc')
            auxshift = real(shiftpp2(sinogramder(:,ii),-shift,0));
        elseif strcmp(interpmeth,'linear')
            auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-shift,0));
        end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
        currenterror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
    
        % compute shift forward error
        if strcmp(interpmeth,'sinc')
            auxshift = real(shiftpp2(sinogramder(:,ii),-(shift+1),0));
        elseif strcmp(interpmeth,'linear')
            auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift+1),0));
        end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
        forwarderror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
        
        % compute shift backward error
        if strcmp(interpmeth,'sinc')
            auxshift = real(shiftpp2(sinogramder(:,ii),-(shift-1),0));
        elseif strcmp(interpmeth,'linear')
            auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift-1),0));
        end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
        backwarderror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
    
    mini = min([currenterror backwarderror forwarderror]);

    switch mini
        case currenterror
            deltastack(1,ii) = deltastack(1,ii) - shift;
            auxtempreg(:,ii) = auxshift; %temporary for debugging
            erroryreg(ii) = currenterror;
            dir = 0;
            continue;
        case backwarderror % Looking backward
            dir = -1;
        case forwarderror % Looking forward
            dir = 1;
    end
        
    if dir~=0,
        shift = shift+dir;
        currenterror = mini;
        do = 1;
    end
        while do==1,
            if strcmp(interpmeth,'sinc')
                auxshift = real(shiftpp2(sinogramder(:,ii),-(shift+dir),0));
            elseif strcmp(interpmeth,'linear')
                auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift+dir),0));
            end
%             if rembias,
%                 [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%             end
            nexterror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
        
        %shift,
            if nexterror<=currenterror
                shift = shift+dir;
                currenterror = nexterror;
            else
                do = 0;
                deltastack(1,ii) = deltastack(1,ii) - shift;
                auxtempreg(:,ii) = auxshift;
                erroryreg(ii) = currenterror;
            end
        end
   
        
    end

    %end
    
    
    %deltastack(1,:) = deltastack(1,:) - round(mean(deltastack(1,:)));

    display(['Final error mectric,  E = ' num2str(sum(erroryreg))])
    
    changey = abs(deltaprev(1,:) - deltastack(1,:));
    display(['Max correction = ' num2str(max(abs(changey)))])
    
    if (max(abs(changey))<max(pixtol,1)),
        domain = 0;
    end
    if count >= maxit,
        domain = 0;
        warning('Maximum number of iterations exceeded, increase maxit')
    end
    if sum(erroryregprev)<sum(erroryreg)
        warning('Last iteration made error worse, keeping previous to last positions')
        deltastack = deltaprev;
        domain = 0;
    end
    if disp>1,
    
        figure(200);
        subplot(2,1,1)
        imagesc(sinogramderorig_nofilt);
        axis xy tight
        %colorbar
        title('Initial Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')

        subplot(2,1,2),
        imagesc(sinogramder_nofilt);
        axis xy tight
        %colorbar
        title('Current Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')
        
        figure(201);
        imagesc(sinogramcompder);
        axis xy tight
        %colorbar
        title('Synthetic Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')
       
        figure(2),
        plot(deltastack')
        drawnow,
        title('Object position')

        % figure(205);
        % plot(erroryinit)
        % hold on,
        % plot(erroryreg,'r')
        % hold off
        % legend('Initial error per projection','Current error per projection')

    end
end
% if pixtol >= 1,
%     % Compute the shifted images
%     if nargout == 2,
%         display('Computing aligned images')
%         for ii = 1:size(regstack,3),
%             switch interpmeth
%                 case 'sinc'
%                     regstack(:,:,ii) = real(shiftpp2(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii)));
%                 case 'linear'
%                     regstack(:,:,ii) = shiftwrapbilinear(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii));
%             end
%             if mod(ii,20) == 0,
%                 display(['Image ' num2str(ii) ' of ' num2str(size(regstack,3))]),
%             end
%         end
%     end
% end
%%
if pixtol<1,
    display('   Subpixel refinement')
    clear auxinit;  % Before main loop
    domain = 1;
    count = 0;
    while domain == 1,
        count = count+1;
        deltaprev = deltastack;
        erroryregprev = erroryreg;
        for ii = 1:size(sinogramder,2)
            if strcmp(interpmeth,'sinc')
                sinogramder(:,ii) = real(shiftpp2(sinogramderorig(:,ii),deltastack(1,ii),0));
            elseif strcmp(interpmeth,'linear')
                sinogramder(:,ii) = shiftwrapbilinear(sinogramderorig(:,ii),deltastack(1,ii),0);
            end
            if disp>1,
                if strcmp(interpmeth,'sinc')
                    sinogramder_nofilt(:,ii) = real(shiftpp2(sinogramderorig_nofilt(:,ii),deltastack(1,ii),0));
                elseif strcmp(interpmeth,'linear')
                    sinogramder_nofilt(:,ii) = shiftwrapbilinear(sinogramderorig_nofilt(:,ii),deltastack(1,ii),0);
                end
            end
        end

        %%%%%% Sinogram alignment subpixel precision
        % Compute tomogram with current sinogramder
%         recons = iradonfast_v3(double(sinogramder),...
%             theta,'linear','derivative','Han',size(sinogramder,1),filtertomo);
        recons = tomo.iradonfast_v3(double(sinogramder),...
            theta,'linear','derivative','Ram-Lak',size(sinogramder,1),1);

        if ~isempty(cliplow)
            recons = recons.*(recons>=cliplow) + cliplow*(recons<cliplow);
        end
        if ~isempty(cliphigh)
            recons = recons.*(recons<=cliphigh) + cliphigh*(recons>cliphigh);
            recons = recons-cliphigh;
        end
        if usecircle
            recons = recons.*circulo;
        end
        figure(203)
        imagesc(recons);
        axis xy equal tight
        colorbar
    
        % Computed sinogram
        sinogramcomp  = radon(recons,theta);
        Nbig = size(sinogramcomp,1);
        centerbig = floor((Nbig+1)/2);
    
        if strcmp(interpmeth,'sinc')
            sinogramcompder = real(shiftpp2(sinogramcomp,0.5,0)-shiftpp2(sinogramcomp,-0.5,0));
        elseif strcmp(interpmeth,'linear')
            sinogramcompder = real(shiftwrapbilinear(sinogramcomp,0.5,0)-shiftwrapbilinear(sinogramcomp,-0.5,0));
        end
        sinogramcompder = sinogramcompder([1:N]+centerbig-center,:);

        % Compare sinogramder with sinogramcompder
        for ii = 1:size(sinogramder,2)
            %ii,
            errorinit(ii) = sum(abs(sinogramder(padval+masklims,ii)-sinogramcompder(padval+masklims,ii)).^2);
        end
        display(['Initial error metric, E = ' num2str(sum(errorinit))])
    
    
        %%% Search for shifts with respect to sinthetic
        for ii = 1:size(sinogramder,2);

            shift = 0;
            %%% Looking bothways
            % compute current shift error
            if strcmp(interpmeth,'sinc')
                auxshift = real(shiftpp2(sinogramder(:,ii),-shift,0));
            elseif strcmp(interpmeth,'linear')
                auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-shift,0));
            end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
            currenterror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
    
            % compute shift forward error
            if strcmp(interpmeth,'sinc')
                auxshift = real(shiftpp2(sinogramder(:,ii),-(shift+pixtol),0));
            elseif strcmp(interpmeth,'linear')
                auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift+pixtol),0));
            end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
            forwarderror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
        
            % compute shift backward error
            if strcmp(interpmeth,'sinc')
                auxshift = real(shiftpp2(sinogramder(:,ii),-(shift-pixtol),0));
            elseif strcmp(interpmeth,'linear')
                auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift-pixtol),0));
            end
%         if rembias,
%             [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%         end
            backwarderror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
    
            mini = min([currenterror backwarderror forwarderror]);

            switch mini
                case currenterror
                    deltastack(1,ii) = deltastack(1,ii) - shift;
                    auxtempreg(:,ii) = auxshift; %temporary for debugging
                    erroryreg(ii) = currenterror;
                    dir = 0;
                    continue;
                case backwarderror % Looking backward
                    dir = -pixtol;
                case forwarderror % Looking forward
                    dir = pixtol;
            end
        
            if dir~=0,
                shift = shift+dir;
                currenterror = mini;
                do = 1;
            end
    
            while do==1,
                if strcmp(interpmeth,'sinc')
                    auxshift = real(shiftpp2(sinogramder(:,ii),-(shift+dir),0));
                elseif strcmp(interpmeth,'linear')
                    auxshift = real(shiftwrapbilinear(sinogramder(:,ii),-(shift+dir),0));
                end
%             if rembias,
%                 [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
%             end
                nexterror = sum(abs(auxshift(padval+masklims)-sinogramcompder(padval+masklims,ii)).^2);
        
        %shift,
                if nexterror<=currenterror
                    shift = shift+dir;
                    currenterror = nexterror;
                else
                    do = 0;
                    deltastack(1,ii) = deltastack(1,ii) - shift;
                    auxtempreg(:,ii) = auxshift;
                    erroryreg(ii) = currenterror;
                end
            end
   
        
        end

    
    
    
    %deltastack(1,:) = deltastack(1,:) - round(mean(deltastack(1,:)));

    display(['Final error mectric,  E = ' num2str(sum(erroryreg))])
    
    changey = abs(deltaprev(1,:) - deltastack(1,:));
    display(['Max correction = ' num2str(max(abs(changey)))])
    
    if (max(abs(changey))<pixtol),
        domain = 0;
    end
    if count >= maxit,
        domain = 0;
        warning('Maximum number of iterations exceeded, increase maxit')
    end
    if sum(erroryregprev)<sum(erroryreg)
        warning('Last iteration made error worse, keeping previous to last positions')
        deltastack = deltaprev;
        domain = 0;
    end
        if disp>1,
    
        figure(200);
        subplot(2,1,1)
        imagesc(sinogramderorig_nofilt);
        axis xy tight
        %colorbar
        title('Initial Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')

        subplot(2,1,2),
        imagesc(sinogramder_nofilt);
        axis xy tight
        %colorbar
        title('Current Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')
        
        figure(201);
        imagesc(sinogramcompder);
        axis xy tight
        %colorbar
        title('Synthetic Sinogram')
        ylabel('x [pixels]')
        xlabel('Projection')
       
        figure(2),
        plot(deltastack')
        drawnow,
        title('Object position')
        end
    end
end

%%%%%
if (binning ~=0)&&(binning ~=1)
    deltastack = deltastack*binning;
    sinogramder = sinogramder_orig;
end

%%% Create ouput aligned sinogram sinogramdernopad outputsinogram
if nargout>1,
    for ii = 1:size(sinogramder,2),
        if strcmp(interpmeth,'sinc')
            outputsinogram(:,ii) = real(shiftpp2(outputsinogram(:,ii),deltastack(1,ii),0));
        elseif strcmp(interpmeth,'linear')
            outputsinogram(:,ii) = shiftwrapbilinear(outputsinogram(:,ii),deltastack(1,ii),0);
        end
    end
end



