% DISTMAT Compute a Distance Matrix for One or Two Sets of Points
%
% 
% 
% Copyright (c) 2015, Joseph Kirk
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
% 

% Filename: distmat.m
%
% Description: Computes a matrix of pair-wise distances between points in
%              A and B, using one of {euclidean,cityblock,chessboard} methods
%
% Author:
%       Joseph Kirk
%       jdkirk630@gmail.com
%
% Date: 02/27/15
%
% Release: 2.0
%
% Inputs:
%     A      - (required) MxD matrix where M is the number of points in D dimensions
%     B      - (optional) NxD matrix where N is the number of points in D dimensions
%                 if not provided, B is set to A by default
%     METHOD - (optional) string specifying one of the following distance methods:
%                 'euclidean'                       Euclidean distance (default)
%                 'taxicab','manhattan','cityblock' Manhattan distance
%                 'chebyshev','chessboard','chess'  Chebyshev distance
%                 'grid','diag'                     Diagonal grid distance
%
% Outputs:
%     DMAT   - MxN matrix of pair-wise distances between points in A and B
%
% Usage:
%     dmat = distmat(a)
%       -or-
%     dmat = distmat(a,b)
%       -or-
%     dmat = distmat(a,method)
%       -or-
%     dmat = distmat(a,b,method)
%
% Example:
%     % Pairwise Euclidean distances within a single set of 2D points
%     xy = 10*rand(25,2);  % 25 points in 2D
%     dmat = distmat(xy);
%     figure; plot(xy(:,1),xy(:,2),'.');
%     for i=1:25, text(xy(i,1),xy(i,2),[' ' num2str(i)]); end
%     figure; imagesc(dmat); colorbar
%
% Example:
%     % Pairwise Manhattan distances within a single set of 2D points
%     xy = 10*rand(25,2);  % 25 points in 2D
%     dmat = distmat(xy,'cityblock');
%     figure; plot(xy(:,1),xy(:,2),'.');
%     for i=1:25, text(xy(i,1),xy(i,2),[' ' num2str(i)]); end
%     figure; imagesc(dmat); colorbar
%
% Example:
%     % Pairwise Chebyshev distances within a single set of 2D points
%     xy = 10*rand(25,2);  % 25 points in 2D
%     dmat = distmat(xy,'chebyshev');
%     figure; plot(xy(:,1),xy(:,2),'.');
%     for i=1:25, text(xy(i,1),xy(i,2),[' ' num2str(i)]); end
%     figure; imagesc(dmat); colorbar
%
% Example:
%     % Inter-point Euclidean distances for 2D points
%     xy = 10*rand(15,2);  % 15 points in 2D
%     uv = 10*rand(25,2);  % 25 points in 2D
%     dmat = distmat(xy,uv);
%     figure; plot(xy(:,1),xy(:,2),'.');
%     for i=1:15, text(xy(i,1),xy(i,2),[' ' num2str(i)]); end
%     figure; plot(uv(:,1),uv(:,2),'.');
%     for i=1:25, text(uv(i,1),uv(i,2),[' ' num2str(i)]); end
%     figure; imagesc(dmat); colorbar
%
% See also:
%
function dmat = distmat(a,varargin)
    
    
    % Set defaults
    method = 'euclidean';
    b = a;
    
    % Error check primary input
    if ~isnumeric(a)
        error('Expecting a matrix of floating point values for A input.');
    end
    
    % Process optional inputs
    for var = varargin
        arg = var{1};
        if ischar(arg)
            method = arg;
        elseif ~isempty(arg)
            b = arg;
        end
    end
    
    % Check input dimensionality
    [na,aDims] = size(a);
    [nb,bDims] = size(b);
    if (aDims ~= bDims)
        error('Input matrices must have the same dimensionality.');
    end
    
    % Create index matrices
    [j,i] = meshgrid(1:nb,1:na);
    
    % Compute array of inter-point differences
    delta = a(i,:) - b(j,:);
    
    % Compute distance by specified method
    dmat = zeros(na,nb);
    switch lower(method)
        case {'euclidean','euclid'}
            % Euclidean distance
            dmat(:) = sqrt(sum(delta.^2,2));
        case {'cityblock','city','block','manhattan','taxicab','taxi'}
            % Cityblock distance
            dmat(:) = sum(abs(delta),2);
        case {'chebyshev','cheby','chessboard','chess'}
            % Chebyshev distance
            dmat(:) = max(abs(delta),[],2);
        case {'grid','diag'}
            dmat(:) = max(abs(delta),[],2) + (sqrt(2) - 1)*min(abs(delta),[],2);
        otherwise
            error('Unrecognized distance method %s',method);
    end
    
end

