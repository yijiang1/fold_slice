% [deltastack regstack] =  alignprojections_v4(regstack,limsy,limsx,deltastack,pixtol,rembias,maxorder,disp,paramsalign)
% Function to align projections. It relies on having air on both sides of
% the sample. The code aligns in x using a center of mass and in y by
% aligning the projections (in x) of the sample. It performs a local search
% in y, so convergence issues can be addressed by giving an approximate initial
% guess for a possible drift.
% WARNING. The code aligns the center of rotation in the ROI defined by
% limsx so that it is consistent (within the window) with the center
% expected by iradon: ceil(size(R,1)/2).
%
%
% Inputs %
% restack       Stack of projections
% limsy         Limits of window of interest in y
% limsx         Limits of window of interest in x
% deltastack    Vectors [y;x] of initial estimates for object motion (2,n)
% pixtol        Tolerance for change in registration
% rembias       true -> removal of bias and lower order terms (for y
%               registration). Default false
% maxorder      if rembias is true specify the order of bias removal (e.g.
%               = 1 mean, = 2 linear). Default = 0.
% disp          Display = 0 no images
%               = 1 Final diagnostic images
%               = 2 Diagnostic images per iteration
%
% Extra optional parameters (in structure paramsalign)
% paramsalign.alignx        = true - align x using center of mass (default),
%                           = false - align y only
% paramsalign.expshift      = false - Shift images normally (default)
%                           = true - shift images in phasor space
% paramsalign.interpmeth    = 'sinc' - Shift images with sinc interpolation (default)
%                           = 'linear' - Shift images with linear interpolation (default)
%
% Outputs %
% deltastack    Object positions
% regstack2     Aligned stack of images (discrete sinc interpolation)
%
% Manuel Guizar 23 Sept 2010
% This code is provided as is and without guarantees on its performance
% The algorithm is not yet published. Do not distribute and contact
% (mguizar@gmail.com) prior to publication of results where this code plays
% a significant role, or to report a bug.
% Please acknowledge if used.
% v3 - Option to enable align y only - Guizar Feb 8, 2011
% v4 - Option to change from sinc to bilinear interpolation
%    - Option to shift final images in phasor space
%    - Fixed bug on computation of parabola for empty frame
%      Guizar June 15, 2011

function [deltastack regstack] =  alignprojections_v4(regstack,limsy,limsx,deltastack,pixtol,rembias,maxorder,disp,paramsalign)

import math.*
import utils.*


deltastack = round(deltastack);
maxit = 15;

if exist('paramsalign') == 0,
    paramsalign = 0;
end

if isfield(paramsalign,'interpmeth') == 0,
    paramsalign.interpmeth = 'sinc';
else
    interpmeth = paramsalign.interpmeth;
    if (strcmp(interpmeth,'sinc'))||(strcmp(interpmeth,'linear'))
    else
        error('Undefined interpolation method')
    end
end

if isfield(paramsalign,'expshift') == 0,
    expshift = false;
else
    expshift = paramsalign.expshift;
end

if isfield(paramsalign,'alignx') == 0,
    alignx = true;
else
    alignx = paramsalign.alignx;
end

if exist('pixtol') == 0,
    pixtol = 1;
end

if exist('deltastack') == 0,
    deltastack = zeros(2,size(regstack,3));
end

if exist('rembias') == 0,
    rembias = false;
end

if exist('maxorder') == 0,
    maxorder = 0;
end

display(['   Registration of projections - Single pixel'])
[Ny,Nx] = size(regstack(limsy(1):limsy(2),limsx(1):limsx(2),1));
Xp = [-fix(Nx/2 -0.5):ceil(Nx/2 + 0.5)-1]; % not right for an fft but right for
% iradonfast
Yp = ([-fix(Ny/2):ceil(Ny/2-1)]);
[Xp,Yp] = meshgrid(Xp,Yp);
figure(1);
imagesc(regstack(:,:,1));
axis xy equal tight
colormap bone
hold on
plot([limsx(1) limsx(1)],[limsy(1) limsy(2)],'r')
plot([limsx(2) limsx(2)],[limsy(1) limsy(2)],'r')
plot([limsx(1) limsx(2)],[limsy(1) limsy(1)],'r')
plot([limsx(1) limsx(2)],[limsy(2) limsy(2)],'r')
hold off

%%%%%%%%%%%%% Single pixel precision %%%%%%%%%%%%%%
%Center of mass x
% loop here
clear auxinit;  % Before main loop
domain = 1;
count = 0;
while domain == 1,
    count = count+1;
    deltaprev = deltastack;
    if alignx
        for ii = 1:size(regstack,3),
            mass(ii) = sum(sum(regstack([limsy(1):limsy(2)]+deltastack(1,ii),...
                [limsx(1):limsx(2)]+deltastack(2,ii),ii),1),2);
            center(ii) = squeeze(sum(sum(Xp.*regstack([limsy(1):limsy(2)]+deltastack(1,ii),...
                [limsx(1):limsx(2)]+deltastack(2,ii),ii),1),2));
        end
        center = center./mass;
    else
        center = 0;
    end
    display(['Max correction, center of mass in x =' num2str(max(abs(center)))])
    % Center for iradonfast =  ceil(size(R,1)/2)
    % Correction with mass center
    deltastack(2,:) = deltastack(2,:) + round(center);
    
    %%%%%% Mass distribution registration in y
    clear aux,
    for ii = 1:size(regstack,3)
        %ii
        aux(:,ii) = squeeze(sum(regstack([limsy(1):limsy(2)]+deltastack(1,ii),...
            [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
        if rembias,
            [coeffs auxi] = projectleg1D_2(aux(:,ii),maxorder,Yp(:,1),1);
            aux(:,ii) = auxi.';
        end
    end
    meanaux = mean(aux,2);
    if exist('auxinit')==0,
        auxinit = aux;
        meaninit = meanaux;
        for ii = 1:size(aux,2)
            erroryinit(ii) = sum(abs(aux(:,ii)-meanaux).^2);
        end
        display(['Initial error metric for y, E = ' num2str(sum(erroryinit))])
    end
    
    
    
    %%% Search for shifts with respect to mean
    for ii = 1:size(regstack,3);
        
        shift = 0;
        %%% Looking bothways
        % compute current shift error
        auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-shift+deltastack(1,ii),...
            [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
        if rembias,
            [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
        end
        currenterror = sum(abs(auxshift-meanaux).^2);
        
        % compute shift forward error
        auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift+1)+deltastack(1,ii),...
            [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
        if rembias,
            [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
        end
        forwarderror = sum(abs(auxshift-meanaux).^2);
        
        % compute shift backward error
        auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift-1)+deltastack(1,ii),...
            [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
        if rembias,
            [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
        end
        backwarderror = sum(abs(auxshift-meanaux).^2);
        
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
            
            while do==1,
                auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift+dir)+deltastack(1,ii),...
                    [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
                if rembias,
                    [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
                end
                nexterror = sum(abs(auxshift-meanaux).^2);
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
        
    end
    
    
    deltastack(1,:) = deltastack(1,:) - round(mean(deltastack(1,:)));
    
    display(['Final error mectric for y,  E = ' num2str(sum(erroryreg))])
    
    changey = abs(deltaprev(1,:) - deltastack(1,:));
    display(['Max correction in y = ' num2str(max(abs(changey)))])
    
    if (max(abs(changey))<max(pixtol,1))&&(max(abs(center))<max(pixtol,1)),
        domain = 0;
    end
    if count >= maxit,
        domain = 0;
        warning('Maximum number of iterations exceeded, increase maxit')
    end
    
    if disp>1,
        
        figure(200);
        subplot(2,1,1)
        imagesc(auxinit);
        axis xy tight
        %colorbar
        title(['Initial Integral in x, maxorder = ' num2str(maxorder)])
        ylabel('y [pixels]')
        xlabel('Projection')
        
        subplot(2,1,2),
        imagesc(auxtempreg);
        axis xy tight
        %colorbar
        title('Current Integral in x')
        ylabel('y [pixels]')
        xlabel('Projection')
        
        figure(201)
        subplot(2,1,1),
        plot(auxinit)
        hold on,
        plot(meaninit,'r','Linewidth',2)
        plot(meaninit,'--w')
        hold off,
        title(['Initial Integral in x, maxorder = ' num2str(maxorder)])
        
        subplot(2,1,2),
        %plot(auxtempreg(:,[1 360])) %28
        plot(auxtempreg)
        hold on,
        plot(meanaux,'r','Linewidth',2)
        plot(meanaux,'w')
        mean2 = mean(auxtempreg,2);
        % plot(mean2,'r','Linewidth',2)
        % plot(mean2,'--w')
        hold off,
        title('Current Integral in x')
        
        
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
if pixtol >= 1,
    % Compute the shifted images
    if nargout == 2,
        display('Computing aligned images')
        for ii = 1:size(regstack,3),
            switch interpmeth
                case 'sinc'
                    regstack(:,:,ii) = real(shiftpp2(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii)));
                case 'linear'
                    regstack(:,:,ii) = shiftwrapbilinear(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii));
            end
            if mod(ii,20) == 0,
                display(['Image ' num2str(ii) ' of ' num2str(size(regstack,3))]),
            end
        end
    end
end

if pixtol<1,
    %warning(['Subpixel alignment is to be implemented'])
    display(['    Subpixel alignment'])
    
    % Compute full massx (integral in y) using rounded deltastack for y window
    % Compute full massy (integral in x) using rounded deltastack for x window
    for ii = 1:size(regstack,3),
        massxorig(:,ii) = squeeze(sum(regstack([limsy(1):limsy(2)]+round(deltastack(1,ii)),:,ii),1));
        massyorig(:,ii) = squeeze(sum(regstack(:,[limsx(1):limsx(2)]+round(deltastack(2,ii)),ii),2));
    end
    
    dosubpix = 1;
    count = 0;
    deltay = 1;
    while dosubpix == 1,
        count = count+1;
        deltaprev = deltastack;
        
        % Shift massx subpixel for current deltastack
        % Recompute center of mass
        if alignx
            for ii = 1:size(regstack,3),
                switch interpmeth
                    case 'sinc'
                        massx(:,ii) = real(shiftpp2(massxorig(:,ii),deltastack(2,ii),0));
                    case 'linear'
                        massx(:,ii) = shiftwrapbilinear(massxorig(:,ii),deltastack(2,ii),0);
                end
                mass(ii) = sum(massx([limsx(1):limsx(2)],ii),1);
                center(ii) = sum(Xp(1,:).'.*massx([limsx(1):limsx(2)],ii),1);
            end
            center = center./mass;
        else
            center = 0;
        end
        display(['Max correction, center of mass in x = ' num2str(max(abs(center)))])
        % Center for iradonfast =  ceil(size(R,1)/2)
        % Correction with mass center
        deltastack(2,:) = deltastack(2,:) + center;
        
        % Compute current shift error
        % Shift massy subpixel for current deltastack
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Mass distribution registration in y %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear aux
        %%% Search for shifts with respect to mean
        for ii = 1:size(regstack,3);
            shift = 0;
            %%% Looking bothways
            % compute current shift error
            %             auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-shift+deltastack(1,ii),...
            %                 [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
            switch interpmeth
                case 'sinc'
                    auxshift = real(shiftpp2(massyorig(:,ii),deltastack(1,ii),0));
                case 'linear'
                    auxshift = shiftwrapbilinear(massyorig(:,ii),deltastack(1,ii),0);
            end
            auxshift = auxshift([limsy(1):limsy(2)],:); % Clip to ROI
            if rembias,
                [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
                massy(:,ii) = auxshift;
            end
            currenterror = sum(abs(auxshift-meanaux).^2);
            
            % compute shift forward error
            %             auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift+1)+deltastack(1,ii),...
            %                 [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
            switch interpmeth
                case 'sinc'
                    auxshift = real(shiftpp2(massyorig(:,ii),deltastack(1,ii)-pixtol,0));
                case 'linear'
                    auxshift = shiftwrapbilinear(massyorig(:,ii),deltastack(1,ii)-pixtol,0);
            end
            auxshift = auxshift([limsy(1):limsy(2)],:); % Clip to ROI
            if rembias,
                [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
            end
            forwarderror = sum(abs(auxshift-meanaux).^2);
            
            % compute shift backward error
            %         auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift-1)+deltastack(1,ii),...
            %             [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
            switch interpmeth
                case 'sinc'
                    auxshift = real(shiftpp2(massyorig(:,ii),deltastack(1,ii)+pixtol,0));
                case 'linear'
                    auxshift = shiftwrapbilinear(massyorig(:,ii),deltastack(1,ii)+pixtol,0);
            end
            auxshift = auxshift([limsy(1):limsy(2)],:); % Clip to ROI
            if rembias,
                [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
            end
            backwarderror = sum(abs(auxshift-meanaux).^2);
            
            mini = min([currenterror backwarderror forwarderror]);
            
            switch mini
                case currenterror
                    deltastack(1,ii) = deltastack(1,ii);
                    %                 auxtempreg(:,ii) = auxshift; %temporary for debugging
                    erroryreg(ii) = currenterror;
                    dir = 0;
                    continue;
                case backwarderror % Looking backward
                    dir = -1;
                case forwarderror % Looking forward
                    dir = 1;
            end
            
            if dir~=0,
                shift = shift+dir*pixtol;
                currenterror = mini;
                do = 1;
                
                while do==1,
                    %                 auxshift = squeeze(sum(regstack([limsy(1):limsy(2)]-(shift+dir)+deltastack(1,ii),...
                    %                     [limsx(1):limsx(2)]+deltastack(2,ii),ii),2));
                    switch interpmeth
                        case 'sinc'
                            auxshift = real(shiftpp2(massyorig(:,ii),deltastack(1,ii)-(shift+dir*pixtol),0));
                        case 'linear'
                            auxshift = shiftwrapbilinear(massyorig(:,ii),deltastack(1,ii)-(shift+dir*pixtol),0);
                    end
                    auxshift = auxshift([limsy(1):limsy(2)],:); % Clip to ROI
                    if rembias,
                        [coeffs auxshift] = projectleg1D_2(auxshift,maxorder,Yp(:,1),1);
                        massy(:,ii) = auxshift;
                    end
                    nexterror = sum(abs(auxshift-meanaux).^2);
                    %shift,
                    if nexterror<=currenterror
                        shift = shift+dir*pixtol;
                        currenterror = nexterror;
                    else
                        do = 0;
                        deltastack(1,ii) = deltastack(1,ii) - shift;
                        auxtempreg(:,ii) = auxshift;
                        erroryreg(ii) = currenterror;
                    end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Up to here it obtained the next estimate %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Evaluate if pixel tolerance is already met
        changey = abs(deltaprev(1,:) - deltastack(1,:));
        display(['Max correction in y = ' num2str(max(abs(changey)))])
        
        if (max(abs(center))<pixtol)&&(max(abs(changey))<pixtol)
            %if (max(abs(changey))<1)&&(max(abs(center))<1),
            dosubpix = 0;
        end
        if count >= maxit,
            dosubpix = 0;
            warning('Maximum number of iterations exceeded, increase maxit')
        end
        
        if disp>1,
            
            figure(200);
            
            subplot(2,1,2),
            imagesc(massy);
            axis xy
            %colorbar
            title('Current Integral in x')
            ylabel('y [pixels]')
            xlabel('Projection')
            
            
            figure(201)
            subplot(2,1,2),
            %plot(auxtempreg(:,[1 360])) %28
            plot(massy)
            hold on,
            plot(meanaux,'r','Linewidth',2)
            plot(meanaux,'w')
            hold off,
            title('Current Integral in x')
            
            
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
    
    
    
    % Compute the shifted images
    if nargout == 2,
        if expshift == 0,
            display('Computing aligned images')
            for ii = 1:size(regstack,3),
                switch interpmeth
                    case 'sinc'
                        regstack(:,:,ii) = real(shiftpp2(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii)));
                    case 'linear'
                        regstack(:,:,ii) = shiftwrapbilinear(regstack(:,:,ii),deltastack(1,ii),deltastack(2,ii));
                end
                if mod(ii,20) == 0,
                    display(['Image ' num2str(ii) ' of ' num2str(size(regstack,3))]),
                end
            end
        elseif expshift == 1,
            display('Computing aligned images in phasor space')
            for ii = 1:size(regstack,3),
                switch interpmeth
                    case 'sinc'
                        regstack(:,:,ii) = angle(shiftpp2(exp(1i*regstack(:,:,ii)),deltastack(1,ii),deltastack(2,ii)));
                    case 'linear'
                        regstack(:,:,ii) = angle(shiftwrapbilinear(exp(1i*regstack(:,:,ii)),deltastack(1,ii),deltastack(2,ii)));
                end
                if mod(ii,20) == 0,
                    display(['Image ' num2str(ii) ' of ' num2str(size(regstack,3))]),
                end
            end
        end
    end
    
    if disp>=1,
        figure(200)
        subplot(2,1,2),
        imagesc(massy);
        axis xy fill
        %colorbar
        title('Current Integral in x')
        ylabel('y [pixels]')
        xlabel('Projection')
        
        figure(201)
        subplot(2,1,2),
        %plot(auxtempreg(:,[1 360])) %28
        plot(massy)
        hold on,
        plot(meanaux,'r','Linewidth',2)
        plot(meanaux,'w')
        
        hold off,
        title('Current Integral in x')
        
        
        figure(2),
        plot(deltastack')
        drawnow,
        title('Object position')
        legend('dy','dx')
    end
end


