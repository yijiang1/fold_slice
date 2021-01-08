function [x, p] = cgmin1(func,x,itmax,ftol,xtol,varargin)
% conjugate-gradient optimization routine
% NOTE: linesearch subroutines do not use the gradient
%
% [x] = cgmin1(func,x,itmax,ftol,xtol,varargin)
%
%     func = string name of objective function which returns both the
%            objective function value and the gradient
%        x = input as initial starting point and output as final point
%    itmax = maximum number of iterations (empty for default = 50)
%     ftol = relative function tolerance (empty for default = 1e-3)
%     xtol = absolute solution tolerance (empty for default = 1e-3)
% varargin = extra variables required by objective function
%
% DISCLAIMER: This code is not intended for distribution. I have many
% versions of this code and am constantly revising it. I believe this
% version is working properly. However, I will not vouch for the code.
% Anyone using the code for thesis research has a responsibility to go
% through the code line-by-line and read relevant references to understand
% the code completely. In my opinion, you have two options if you want to
% publish results obtained with the code: (i) go through the code line-by-
% line and read relevent references to understand how the code works and make
% sure it is working properly for your application, or (ii) I can sit down
% with you an go through this code and the additional code that you have
% written to go along with it and make sure it is working properly. Option
% (i) is preferred, and I ask that you do NOT acknowledge me in print (first,
% it would be more appropriate for you to reference "Numerical Recipes",
% and second, I prefer not to be named in a paper with which I do not have
% detailed knowledge). If you decide to go with option (ii), I would expect
% to learn the details of your research and be included in the author list.
%
% Sam Thurman, May 9, 2005

import utils.*

if isempty(itmax), itmax = 50; end
if isempty(ftol), ftol = 1e-3; end
if isempty(xtol), xtol = 1e-3; end
% loop
flg = 0; % use steepest descent for first iteration
step = 0; % to guess at initial steplength
for it = 1:itmax
   % function evaluation
   [f,grad,p] = feval(func,x,varargin{:});
%    disp(f)
   % check for feasibility
   if isinf(f), error('encountered an infeasible solution'), end
   if norm(grad(:))==0, return, end % done if gradient is zero (unlikely)
   % pick search direction
   if (flg==1) & (rem(it,25)~=0) % linesearch found a minimum -> use cg equations
       gg = g(:)'*g(:);
%        dgg = grad(:)'*grad(:); % this statement for Fletcher-Reeves
       dgg = (grad(:)+g(:))'*grad(:);  % this statement for Polak-Ribiere
       ga = dgg/gg;
       g = -grad;
       h = g+ga*h;
       dx = h/norm(h(:));
       df = grad(:)'*dx(:);
   end
   if (flg==0) | (rem(it,25)==0) | (df>0) % revert to steepest decent
       g = -grad;
       h = g;
       dx = h/norm(h(:));
       df = grad(:)'*dx(:);
   end
   % initial steplength guess
   if step == 0
       step = max(0.001,min([1,2*abs(f/(grad(:)'*dx(:)))])); % same as fminusub.m (line 124) in optim toolbox
   else % oterwise use previous steplength
        step = step/10;
   end
   % linesearch
   [x,fvalue,step,flg] = engines.ML.linesearch(func,x,f,df,dx,step,varargin{:});
   % test for convergence
   if (2*abs(f-fvalue)<=ftol*(abs(f)+abs(fvalue)+ftol)) & (step*norm(dx(:))<=xtol) & (it~=1) % normal return
       return
   end
end
verbose(3, 'Maximum number of iterations exceeded.')
return
end

