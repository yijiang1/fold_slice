function [a,f] = brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin)
% one-dimensional minimization by parabolic interpolation & golden
% section (does not use the gradient)
%
% [a,f] = brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin)
%
% func = string name of objective function
% x0 = starting point of linesearch
% dx = direction of linesearch
% a1,a,a2 = bracketing triplet of steplengths (a1<a<a2)
% f1,f,f2 = objective function at steplengths a1, a, & a2
% a = output as final steplength
% f = output as objective function at final steplength

gold = 0.3819660; % golden ratio
itmax = 5;
tol = 0.5;
% check order of bracket
if (a1>a)|(a2<a), error('brent called with bracket in wrong order'), end
% initialize
v = a;fv = f; % middle point on step before last
w = a;fw = f; % middle point on last step
e = 0; % distance moved on step before last
% iterations
for it = 1:itmax
   am = 0.5*(a1+a2);
   tol1 = tol*abs(a)+eps;
   tol2 = 2*tol1;
   % test for convergence
   if abs(a-am)<=(tol2-0.5*(a2-a1)), return, end
   % choose next point
   if abs(e)>tol1 % construct a trial parabolic fit
       r = (a-w)*(f-fv);
       q = (a-v)*(f-fw);
       p = (a-v)*q-(a-w)*r;
       q = 2*(q-r);
       if q>0, p = -p; end
       q = abs(q);
       etemp = e;
       e = d;
       % check acceptability of parabolic fit
       ok = ~(abs(p)>=abs(0.5*q*etemp) | p<=q*(a1-a) | p>=q*(a2-a));
       if ok % take parabolic step
           d = p/q;
           u = a+d;
           if (u-a1)<tol2 | (a2-u)<tol2, d = sign(am-a)*tol1; end
       else % take golden section step
           if a>=am
               e = a1-a;
           else
               e = a2-a;
           end
           d = gold*e;
       end
   else % take golden section step
       if a>=am
           e = a1-a;
       else
           e = a2-a;
       end
       d = gold*e;
   end
   % arrive here with d computed either from
   % parabolic fit or else from golden section
   if abs(d)>=tol1
       u = a+d;
   else
       u = a+sign(d)*tol1;
   end
   fu = feval(func,x0+u*dx,varargin{:}); % one function evaluation per iteration
   if fu<=f
       if u>=a
           a1 = a;
       else
           a2 = a;
       end
       v = w; fv = fw;
       w = a; fw = f;
       a = u; f = fu;
   else
       if u<a
           a1 = u;
       else
           a2 = u;
       end
       if fu<=fw | w==a
           v = w; fv = fw;
           w = u; fw = fu;
       elseif fu<=fv | v==a | v==w
           v = u; fv = fu;
       end
   end
end
disp('exceeded maximum number of iterations')
return
end

