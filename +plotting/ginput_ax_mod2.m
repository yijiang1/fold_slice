function varargout = ginput_ax_mod2(ha,n)
if nargin<2
    n=1;
end
k = 0;
button = 0;

%%// Tolerance so that in the inifnity case, this could act as
%%// the thresholding distance below which the
%%// input extracting operation must be terminated
TOL = 0.01;

%%// Placeholders for X-Y and button type could be stored 
button1 = [];
xy = [];

hf = get(ha,'parent');
figure(hf);
set(hf,'WindowButtonMotionFcn',@changepointer)
set(ha,'ButtonDownFcn',@getpoints)
hp = get(ha,'children');
ht = get(hp,'hittest');
set(hp,'hittest','off')
axlim = get(ha,'Position');
fglim = get(hf,'Position');
x1 = axlim(1)*fglim(3) + fglim(1);
x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
y1 = axlim(2)*fglim(4) + fglim(2);
y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
waitfor(hf,'WindowButtonMotionFcn',[])
if iscell(ht)
    for jj=1:length(ht)
        set(hp(jj),'hittest',ht{jj})
    end
else
    set(hp,'hittest',ht)
end
selType = get(hf,'SelectionType');

% Mouse-Button recognition...
if(strcmp(button, 'normal'))
    button = 1; % left
elseif(strcmp(button, 'extend'))
    button = 2; % right
elseif(strcmp(button, 'alt'))
    button = 3; % middle
else
    button = 4; % double click any mousebutton
end

if nargout==3
    varargout{1} = xy(:,1);
    varargout{2} = xy(:,2);
    varargout{3} = button1(:,1);
elseif nargout==2
    varargout{1} = xy(:,1);
    varargout{2} = xy(:,2);
else
    varargout{1} = xy;
end
    function changepointer(~,~)
        pntr = get(0,'PointerLocation');
        if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
            set(hf,'Pointer','crosshair')
        else
            set(hf,'Pointer','arrow')
        end
    end
    function getpoints(src,evnt)
        cp = get(src,'CurrentPoint');
        button = get(hf, 'SelectionType');
        k = k+1;

        if k==1
            xy = [xy ;cp(1,1:2)];
            button1 = [button1; {button}];
        end

        if k>=2
            if pdist2(cp(1,1:2),xy(k-1,:))<TOL && isinf(n)
                k = n;
            else
                xy = [xy ;cp(1,1:2)];
                button1 = [button1; {button}];
            end
        end
        if k==n
            set(hf,'Pointer','arrow')
            set(hf,'WindowButtonMotionFcn',[])
            set(ha,'ButtonDownFcn',[])
            return;
        end
    end
end