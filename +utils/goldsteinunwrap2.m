% Implementation of Goldstein unwrap algorithm based on location of
% residues and introduction of branchcuts.
% R. M. Goldstein, H. A. Zebker and C. L. Werner, Radio Science 23, 713-720
% (1988).
% Inputs
% fase      Phase in radians, wrapped between (-pi,pi)
% disp      (optional) = 1 to show progress (will slow down code)
%           will also display the branch cuts
% start     (optional) [y,x] position to start unwrapping. Typically faster 
%           at the center of the array
% Outputs
% faserecon Unwrapped phase ( = fase where phase could not be unwrapped)
% shadow    = 1 where phase could not be unwrapped 
% 31 August, 2010 - Acknowledge if used

% Modified 20 Sept 2010 - Find a safe area to unwrap around the first point

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function [faserecon shadow] = goldsteinunwrap2(fase,disp,start)


display('Unwrapping with Goldstein algorithm')
[nr nc] = size(fase);
if nargin < 2
    disp = 0;
end
if nargin <3
    nrstart = round(nr/2);
    ncstart = round(nc/2);
else
    nrstart = start(1);
    ncstart = start(2);
end

residues =            wrapToPi(fase(2:end,1:end-1)   - fase(1:end-1,1:end-1));
residues = residues + wrapToPi(fase(2:end,2:end)     - fase(2:end,1:end-1));
residues = residues + wrapToPi(fase(1:end-1,2:end)   - fase(2:end,2:end));
residues = residues + wrapToPi(fase(1:end-1,1:end-1) - fase(1:end-1,2:end));
residues = residues/(2*pi);
%%% Find residues
[posr,posc] = find(round(residues)==1);
respos = [posr posc ones(length(posr),1)];
[posr,posc] = find(round(residues)==-1);
resneg = [posr posc -ones(length(posr),1)];
%[posr,posc] = find(round(residues)~=0);
%res = [posr posc];
%res = [respos;resneg];
nres = length(respos(:,1))+length(resneg(:,1));
display(['Found ' num2str(nres) ' residues'])

if nres == 0,
    faserecon = unwrap(unwrap(fase')');
    shadow = faserecon*0;
    return;
end

%%% Find minimum length walls
%currentwall = residues*0;
currentwall = zeros(nr+2,nc+2);
%currentwallcharge = 0;
%wallsegdone = 0;

%currentwall(res(1,1)+1,res(1,2)+1) = 1;


for ii = 1:min(length(respos(:,1)),length(resneg(:,1))),
    dist = (respos(1,1) - resneg(:,1)).^2 + (respos(1,2)-resneg(:,2)).^2;
    ind = find(dist == min(dist),1,'first');
    if sqrt(dist(ind)) < min(nc,nr)/4,%/4
        
        currentwall( respos(1,1)+1,min(respos(1,2),resneg(ind,2))+1 : max(respos(1,2),resneg(ind,2))+1 ) = 1;
        currentwall(min(resneg(ind,1),respos(1,1))+1:max(resneg(ind,1),respos(1,1))+1,resneg(ind,2)+1) = 1;
    
        respos = respos(2:end,:); % Remove from respos
        resaux = resneg(1:ind-1,:);
        resaux = [resaux;resneg(ind+1:end,:)];
        resneg = resaux; 
    else  % Wall too long between them, send to window edge
        % for respos
        distedges = [nr-respos(1,1) respos(1,1) nc-respos(1,2) respos(1,2)]; %upper, lower, right, left
        switch min(distedges)
            case distedges(1) %upper
                currentwall(respos(1,1)+1:nr+2, respos(1,2)+1 ) = 1;
            case distedges(2) %lower
                currentwall(1:respos(1,1)+1, respos(1,2)+1) = 1;  
            case distedges(3); %right
                currentwall(respos(1,1)+1, respos(1,2)+1:nc+2) = 1;
            case distedges(4); %left
                currentwall(respos(1,1)+1, 1:respos(1,2)+1) = 1;
        end
        
        % for resneg
        distedges = [nr-resneg(ind,1) resneg(ind,1) nc-resneg(ind,2) resneg(ind,2)]; %upper, lower, right, left
        switch min(distedges)
            case distedges(1) %upper
                currentwall(resneg(ind,1)+1:nr+2, resneg(ind,2)+1 ) = 1;
            case distedges(2) %lower
                currentwall(1:resneg(ind,1)+1, resneg(ind,2)+1) = 1;  
            case distedges(3); %right
                currentwall(resneg(ind,1)+1, resneg(ind,2)+1:nc+2) = 1;
            case distedges(4); %left
                currentwall(resneg(ind,1)+1, 1:resneg(ind,2)+1) = 1;
        end   
            
    end
end 
% else
%     error('Need to implement for unbalanced charge residues')
% end


% Branch cuts for unpaired residues
res = [respos; resneg];
display([num2str(length(res(:,1))) ' unpaired residues'])
for ii = 1:length(res(:,1)),
    distedges = [nr-res(1,1) res(1,1) nc-res(1,2) res(1,2)]; %upper, lower, right, left
    switch min(distedges)
        case distedges(1) %upper
            currentwall(res(1,1)+1:nr+2, res(1,2)+1 ) = 1;
        case distedges(2) %lower
            currentwall(1:res(1,1)+1, res(1,2)+1) = 1;  
        case distedges(3); %right
            currentwall(res(1,1)+1, res(1,2)+1:nc+2) = 1;
        case distedges(4); %left
            currentwall(res(1,1)+1, 1:res(1,2)+1) = 1;
    end
    res = res(2:end,:);
end
if disp == 1,
figure(4);
imagesc(currentwall);
colorbar
axis xy
title('Branch cuts')
colormap gray
drawnow
end

%% Safe unwrap from start position (this could be made faster)
% Only defined the maximum square, could be made faster by defining a
% rectangle for example
[wallposy wallposx] = find(currentwall == 1); % finds wall positions
%distnearest = (wallposy-nrstart).^2+(wallposx-ncstart).^2;
%distnearest = abs(wallposy-nrstart)+abs(wallposx-ncstart);
distnearest = max(abs(wallposy-nrstart),abs(wallposx-ncstart));
indi = find(distnearest == min(distnearest),1);
longi = min(distnearest)-2;
%longi = min([longi nrstart-1 ncstart-1 nc-ncstart-1 nr-nrstart-1]);

% figure(100); 
% plot(wallposx,wallposy,'o'); 
% hold on, 
% plot(ncstart,nrstart,'or'), 
% plot([-longi longi]+ncstart,[-longi -longi]+nrstart,'-r'); 
% plot([-longi longi]+ncstart,[longi longi]+nrstart,'-r'); 
% plot([-longi -longi]+ncstart,[longi -longi]+nrstart,'-r'); 
% plot([longi longi]+ncstart,[longi -longi]+nrstart,'-r'); 
% hold off,
%%


faserecon = fase*0;
shadow = faserecon+1; % not unwrapped yet

% faserecon(nrstart,ncstart) = fase(nrstart,ncstart);
% shadow(nrstart,ncstart) = 0;
% counter = 0;

% faserecon(nrstart+[-longi:longi],ncstart+[-longi:longi]) ...
%     = unwrap(unwrap(  fase(nrstart+[-longi:longi],ncstart+[-longi:longi])')');
% shadow(nrstart+[-longi:longi],ncstart+[-longi:longi]) = 0;

xmask = [max(1,ncstart-longi):min(nc,ncstart+longi)];
ymask = [max(1,nrstart-longi):min(nr,nrstart+longi)];   
faserecon(ymask,xmask) = unwrap(unwrap(  fase(ymask,xmask)')');
shadow(ymask,xmask) = 0;



counter = 0;


% Start unwrapping
maxiter = 2*max(nr,nc);
wallvert = currentwall(1:end-1,:)&currentwall(2:end,:); % prevents horizontal integration
%wallvert = [zeros(1,nc-1);wallvert;zeros(1,nc-1)];
wallhor = currentwall(:,1:end-1)&currentwall(:,2:end); % prevents horizontal integration
%wallhor = [zeros(nr-1,1) wallhor zeros(nr-1,1)];


while (counter <maxiter)&&(max(shadow(:))==1);
    shadowprev = shadow;
    %%%%% Step right
    %newrec =  [zeros(nr,1) shadow(:,2:end)-shadow(:,1:end-1)] == 1;
    newrec =  [false(nr,1) shadow(:,2:end)&not(shadow(:,1:end-1))];
    %prev = [zeros(nr,1) shadow(:,2:end)-shadow(:,1:end-1)]   == -1;
    % Block forbiden paths here
    newrec(:,2:end) = newrec(:,2:end)&(1-wallvert(1:end-1,2:end-2));
    deltafase = [zeros(nr,1) fase(:,2:end)-faserecon(:,1:end-1)].*newrec;
    %faserecon = faserecon + (fase - round(deltafase/(2*pi))*2*pi).*newrec;
    faserecon(newrec) = fase(newrec) - round(deltafase(newrec)/(2*pi))*2*pi;
    shadow(newrec) = 0;
    
    %%%%% Step left
    %newrec =  [shadow(:,1:end-1)-shadow(:,2:end) zeros(nr,1)] == 1;
    newrec =  [shadow(:,1:end-1)&not(shadow(:,2:end)) false(nr,1)];
    %prev = [zeros(nr,1) shadow(:,2:end)-shadow(:,1:end-1)]   == -1;
    % Block forbiden paths here
    newrec(:,1:end-1) = newrec(:,1:end-1)&(1-wallvert(1:end-1,2:end-2));
    deltafase = [fase(:,1:end-1)-faserecon(:,2:end) zeros(nr,1)].*newrec;
    %faserecon = faserecon + (fase - round(deltafase/(2*pi))*2*pi).*newrec;
    faserecon(newrec) = fase(newrec) - round(deltafase(newrec)/(2*pi))*2*pi;
    shadow(newrec) = 0;
        
    %%%%% Step up (positive y)
    %newrec =  [zeros(1,nc) ; shadow(2:end,:)-shadow(1:end-1,:)] == 1;
    newrec =  [false(1,nc) ; shadow(2:end,:)&not(shadow(1:end-1,:))];
    %prev = [zeros(nr,1) shadow(:,2:end)-shadow(:,1:end-1)]   == -1;
    % Block forbiden paths here
    newrec(2:end,:) = newrec(2:end,:)&(1-wallhor(2:end-2,1:end-1));
    deltafase = [zeros(1,nc) ; fase(2:end,:)-faserecon(1:end-1,:)].*newrec;
    %faserecon = faserecon + (fase - round(deltafase/(2*pi))*2*pi).*newrec;
    faserecon(newrec) = fase(newrec) - round(deltafase(newrec)/(2*pi))*2*pi;
    shadow(newrec) = 0;
    
    %%%%% Step down (negative y)
    %newrec =  [shadow(1:end-1,:)-shadow(2:end,:) ; zeros(1,nc)] == 1;
    newrec =  [shadow(1:end-1,:)&not(shadow(2:end,:)) ; false(1,nc)];% Logical input does not seeem to help with computing time
    %prev = [zeros(nr,1) shadow(:,2:end)-shadow(:,1:end-1)]   == -1;
    % Block forbiden paths here
    newrec(1:end-1,:) = newrec(1:end-1,:)&(1-wallhor(2:end-2,1:end-1));
    deltafase = [fase(1:end-1,:)-faserecon(2:end,:) ; zeros(1,nc)].*newrec;
    %faserecon = faserecon + (fase - round(deltafase/(2*pi))*2*pi).*newrec;
    faserecon(newrec) = fase(newrec) - round(deltafase(newrec)/(2*pi))*2*pi;
    shadow(newrec) = 0;
    
    counter = counter+1;
    
    if any(not(shadow(:)==shadowprev(:))) == 0,
        warning('Not all points are accessible for integration')
        faserecon(shadow==1) = fase(shadow==1);
        break;
    end
    
    
if disp == 1,
    figure(5);
    imagesc(faserecon);
    colorbar
    axis xy
    title('Reconstructed phase')
    colormap jet
    drawnow;
end
    
end

if counter == maxiter,
    warning('Maximum number of iterations exceeded for unwrapping. Increase maxiter.'),
end