%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GoldsteinUnwrap2D implements 2D Goldstein branch cut phase unwrapping algorithm.
%
% References::
% 1. R. M. Goldstein, H. A. Zebken, and C. L. Werner, �Satellite radar interferometry:
%    Two-dimensional phase unwrapping,� Radio Sci., vol. 23, no. 4, pp. 713�720, 1988.
% 2. D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
%    Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
% Posted by Bruce Spottiswoode on 22 December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2008, Bruce Spottiswoode
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% * Neither the name of the University of Cape Town nor the names
% of its contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [unph] = goldstein_unwrap2D(a, max_box_radius)



IM=a;
IM_mask=ones(size(IM));                     %Mask (if applicable)

IM_mag=abs(IM);                             %Magnitude image
IM_phase=angle(IM);                         %Phase image

% Unwrap
residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

unph=IM_unwrapped;

end



