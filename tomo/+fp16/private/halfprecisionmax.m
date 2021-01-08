% halfprecisionmax returns IEEE 754r bit pattern of max half precision value
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    halfprecisionmax
%  Filename:    halfprecisionmax.m
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
% Type 'help halfprecision' to get details of this bit pattern
%
%******************************************************************************

function h = halfprecisionmax
if nargin ~= 0
    error(nargchk(0, 0, nargin));
end
if nargout > 1
    error(nargoutchk(0, 1, nargout));
end
h = uint16(hex2dec('7BFF'));
end
