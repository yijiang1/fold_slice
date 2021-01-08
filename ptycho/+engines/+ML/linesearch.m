% linesearch routine (does not use the gradient)
%
% [x,f,a,flg] = linesearch(func,x0,f0,df0,dx,a,varargin)
%
% func = string name of objective function
% x0 = starting point of search
% f0 = objective function at x0
% df0 = derivative of objective function along dx at x0
% dx = direction of linesearch
% a = steplength input as guess output as taken
% x = final point of search
% f = objective function at final point
% flg = indicates how step was determined (0 for Armijo step, 1 for
%        bracketing and refining a minimum)

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [x,f,a,flg] = linesearch(func,x0,f0,df0,dx,a,varargin)

% check if descent direction
if df0>=0
   warning('linesearch called w/o descent direction')
   x = x0; f = f0; fcount = [0,0]; flg = 0;
   return
end
% first point
a1 = 0; f1 = f0;
% try initial steplength
f = feval(func,x0+a*dx,varargin{:}); % no gradient returned
% keep track of old values & hopefully bracket a minimum
a2 = a; f2 = f;
% make sure initial step is feasible
while isinf(f)
   a2 = a;
   a = 0.25*a;
   f = feval(func,x0+a*dx,varargin{:});
end
% decide what to do next based on Armijo condition
b = 0; % parameter in Armijo condition (>=0)
% if f does not satisfy Armijo condition (steplength may be too
% large) -> decrease steplength until Armijo is satisfied
if f>f0+a*b*df0
   while f>f0+a*b*df0
       if isinf(f2)||(f<=f2)
           a2 = a; f2 = f;
       end
       a = 0.25*a; % decrease steplength
       f = feval(func,x0+a*dx,varargin{:});
   end
end
% arrive here with a1=a0=0 and f<f1 (at least)
tmp = 1;
while f2<=f % try doubling a2
   a2 = 2*a2;
   f2 = feval(func,x0+a2*dx,varargin{:});
   if f2<f
       a1 = a; f1 = f;
       a = a2; f = f2;
   end
   if tmp==5; break; else tmp = tmp+1; end
end
% case where we should have a bracket, but f2 is infinite
while isinf(f2)
%    disp('should have a bracket but f2 is infinite')
   u = a+0.25*(a2-a); % point between a and a2
   fu = feval(func,x0+u*dx,varargin{:});
   if fu<f
       a1 = a; f1 = f;
       a = u;  f = fu;
   else
       a2 = u; f2 = fu;
   end
end
% last steps
if (f<f1)&&(f<f2)&&isfinite(f2) % bracketing successful -> refine minimum
   [a,f] = engines.ML.brent(func,x0,dx,a1,a,a2,f1,f,f2,varargin{:});
   flg = 1; % use conjugate gradient next loop
else % bracketing unsuccessful -> stop
   flg = 0; % use steepest descent next loop
end
x = x0+a*dx;
return
end
