% halfprecision converts IEEE 754 floating point to half precision IEEE 754r
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    halfprecision
%  Filename:    halfprecision.c
%  Programmer:  James Tursa
%  Version:     1.0
%  Date:        March 3, 2009
%  Copyright:   (c) 2009 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
% 
%  halfprecision converts the input argument to/from a half precision floating
%  point bit pattern corresponding to IEEE 754r. The bit pattern is stored in a
%  uint16 class variable. Please note that halfprecision is *not* a class. That
%  is, you cannot do any arithmetic with the half precision bit patterns.
%  halfprecision is simply a function that converts the IEEE 754r half precision
%  bit pattern to/from other numeric MATLAB variables. You can, however, take
%  the half precision bit patterns, convert them to single or double, do the
%  operation, and then convert the result back manually.
% 
%  1 bit sign bit
%  5 bits exponent, biased by 15
%  10 bits mantissa, hidden leading bit, normalized to 1.0
% 
%  Special floating point bit patterns recognized and supported:
% 
%  All exponent bits zero:
%  - If all mantissa bits are zero, then number is zero (possibly signed)
%  - Otherwise, number is a denormalized bit pattern
% 
%  All exponent bits set to 1:
%  - If all mantissa bits are zero, then number is +Infinity or -Infinity
%  - Otherwise, number is NaN (Not a Number)
% 
%  Building:
% 
%  halfprecision requires that a mex routine be built (one time only). This
%  process is typically self-building the first time you call the function
%  as long as you have the files halfprecision.m and halfprecision.c in the
%  same directory somewhere on the MATLAB path. If you need to manually build
%  the mex function, here are the commands:
%
%  >> mex -setup
%    (then follow instructions to select a C / C++ compiler of your choice)
%  >> mex halfprecision.c
%
%  If you have an older version of MATLAB, you may need to use this command:
%
%  >> mex -DDEFINEMWSIZE halfprecision.c
% 
%  Syntax
% 
%  B = halfprecision(A)
%  C = halfprecision(B,S)
%      halfprecision(B,'disp')
% 
%  Description
% 
%  A = a MATLAB numeric array, char array, or logical array.
%
%  B = the variable A converted into half precision floating point bit pattern.
%      The bit pattern will be returned as a uint16 class variable. The values
%      displayed are simply the bit pattern interpreted as if it were an unsigned
%      16-bit integer. To see the halfprecision values, use the 'disp' option, which
%      simply converts the bit patterns into a single class and then displays them.
%
%  C = the half precision floating point bit pattern in B converted into class S.
%      B must be a uint16 or int16 class variable.
%
%  S = char string naming the desired class (e.g., 'single', 'int32', etc.)
%      If S = 'disp', then the floating point bit values are simply displayed.
% 
%  Examples
%  
%  >> a = [-inf -1e30 -1.2 NaN 1.2 1e30 inf]
%  a =
%  1.0e+030 *
%      -Inf   -1.0000   -0.0000       NaN    0.0000    1.0000       Inf
% 
%  >> b = halfprecision(a)
%  b =
%  64512  64512  48333  65024  15565  31744  31744
% 
%  >> halfprecision(b,'disp')
%      -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
% 
%  >> halfprecision(b,'double')
%  ans =
%      -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
% 
%  >> 2^(-24)
%  ans =
%  5.9605e-008
% 
%  >> halfprecision(ans)
%  ans =
%      1
% 
%  >> halfprecision(ans,'disp')
%  5.9605e-008
% 
%  >> 2^(-25)
%  ans =
%  2.9802e-008
% 
%  >> halfprecision(ans)
%  ans =
%      1
% 
%  >> halfprecision(ans,'disp')
%  5.9605e-008
% 
%  >> 2^(-26)
%  ans =
%   1.4901e-008
% 
%  >> halfprecision(ans)
%  ans =
%      0
% 
%  >> halfprecision(ans,'disp')
%     0
% 
%  Note that the special cases of -Inf, +Inf, and NaN are handled correctly.
%  Also, note that the -1e30 and 1e30 values overflow the half precision format
%  and are converted into half precision -Inf and +Inf values, and stay that
%  way when they are converted back into doubles.
% 
%  For the denormalized cases, note that 2^(-24) is the smallest number that can
%  be represented in half precision exactly. 2^(-25) will convert to 2^(-24)
%  because of the rounding algorithm used, and 2^(-26) is too small and underflows
%  to zero.
% 
%**************************************************************************

function varargout = halfprecision(varargin)
disp(' ');
disp('You must build the mex routine before you can use halfprecision.');
disp('Attempting to do so now ...');
disp(' ');
mname = mfilename('fullpath');
cname = [mname '.c'];
if( isempty(dir(cname)) )
    disp('Cannot find the file halfprecision.c in the same directory as the');
    disp('file halfprecision.m. Please ensure that they are in the same');
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error('Unable to compile halprecision.c');
else
    disp(['Found file halfprecision.c in ' cname]);
    disp(' ');
    disp('Now attempting to compile ...');
    disp('(If prompted, please press the Enter key and then select any C/C++');
    disp('compiler that is available, such as lcc.)');
    disp(' ');
    disp(['mex(''' cname ''')']);
    disp(' ');
    try
        mex(cname);
        disp('mex halfprecision.c build completed ... you may now use halfprecision.');
        disp(' ');
    catch
        disp(' ');
        disp('Well, *that* didn''t work ... now trying it with mwSize defined ...');
        disp(' ');
        try
            disp(' ');
            disp(['mex(''-DDEFINEMWSIZE'',''' cname ''')']);
            disp(' ');
            mex('-DDEFINEMWSIZE',cname);
            disp('mex halfprecision.c build completed ... you may now use halfprecision.');
            disp(' ');
        catch
            disp('Hmmm ... That didn''t work either.');
            disp(' ');
            disp('The mex command failed. This may be because you have already run');
            disp('mex -setup and selected a non-C compiler, such as Fortran. If this');
            disp('is the case, then rerun mex -setup and select a C/C++ compiler.');
            disp(' ');
            error('Unable to compile halprecision.c');
        end
    end
end
if false
    varargout = varargin; % Get rid of the lint message
end
end
