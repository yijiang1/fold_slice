% wrapper class around sharedmatrix function , it tries to go around so bugs and
% crashes of the code and safely use it for fast interprocess sharing of
% data 
% Use: 
%   self = shm(protected = false, shm_key = randi(1e9))
%
% COMPILE: 	mex -R2017b @shm/private/sharedmatrix.c -output @shm/private/sharedmatrix
% CLEAN all allocated memory :
%  !ipcs -m | cut -d' ' -f2 | grep '^[0-9]' | while read x; do ipcrm -m $x; done
%
%     % EXAMPLE:
%     x = randn(10,'single')+1i;
%     s = shm();
%     s.upload(x)
%     s.detach;   % locally forget the data 
%     [s, data] = s.attach % load them back from shared memory 
%     disp(s)     % show the SHM class 
%     disp(data)  
%     clear s    % or s.detach
%     disp(data-x)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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



classdef shm  < handle

    properties (GetAccess = 'public', SetAccess = 'private')
        shm_key;
    end
    properties (GetAccess = 'public', SetAccess = 'public')
        protected; % dont delete shared variable when cleared 
        data_handle;
    end
    properties (GetAccess = 'private', SetAccess = 'private')
        isComplex; 
    end
   
    %% METHODS:

    methods ( Access = 'public' )
        function self = shm(varargin)
            % generate a key for the shared memory, use random number to
            % avoid collisions between several processes 
            if nargin < 2 
                self.shm_key = randi(1e9);
            else
                self.shm_key = varargin{2};
            end
            self.isComplex = false;
            
            % make sure that @shm is under all curcumstancess seen by
            % matlab, otherwise the share memory may result in matlab crash
            fpath =  mfilename('fullpath'); 
            fpath = fpath(1:end-8); 
            if ~contains(path, fpath)
                addpath(fpath); 
            end
            
            if nargin < 1
                self.protected = false; 
            else
                self.protected = varargin{1}; 
            end
        end
        function disp(self)
            if self.isattached()
                fprintf(['Size of stored data : %i', repmat('x%i',1,ndims(self.data_handle{1})-1),'\n'], size(self.data_handle{1}))
                fprintf('Class of stored data: %s \n', class(self.data_handle{1}))
            else
               disp('No attached data') 
            end
            fprintf('Protected: %i \n', self.protected)
            fprintf('Complex: %i \n', self.isComplex)

            !ipcs -m
            % clear all allocated memory 
            % !ipcs -m | cut -d' ' -f2 | grep '^[0-9]' | while read x; do ipcrm -m $x; done
        end

        function upload(self, data)
            % only create the memory space and upload there the data 
%             if isnumeric(data)
                assert(isnumeric(data) || islogical(data), 'Only numeric arrays are tested to work safely')
                assert(~isa(data, 'double') || isscalar(data), 'Use single precision to make it fast/memory effecient')
                assert(numel(data)*4 < 100e9, 'Datasets than 100GB may result in failure, TESTME')
                % gather the data from GPU 
                if isa(data, 'gpuArray')
                    data = gather(data);
                end
                if  isscalar(data) && isa(data, 'double')
                    data = single(data);
                end
%             end
            self.free_safe(); % clear the memory if it was used 
            self.create_shm(data, true);  % upload to storage 
        end
        function allocate(self, data)
            % only create the memory space, do not write any data 
            assert(isnumeric(data), 'Only numeric arrays are tested to work safely')
            assert(~isa(data, 'double') , 'Use single/integer precision to make it fast/memory effecient')
            assert(numel(data)*4 < 100e9, 'Datasets than 100GB may result in failure, TESTME')

            self.free_safe(); % clear the memory if it was used 
            self.create_shm(data, false); % use only the array to allocate the storage 
        end
        
        function [self, data] = attach(self)
            assert(nargout < 3 && nargout > 0, 'Number of output argument has to be 1 or 2')
            if ~self.isattached()
                try
                    self.data_handle=sharedmatrix('attach',self.shm_key);
                    
                catch err 
                    if strcmpi(err.identifier,'MATLAB:sharedmatrix:attach')
                        error('No shared data available at address: %i', self.shm_key)
                    end
                    rethrow(err)
                end
            end
            if ~isempty(self.data_handle{2})
                %% very bad option because memcpy is enforced, but the only solution till sharedmatrix.c is properly fixed 
                data = complex(self.data_handle{:});
            else
                data = self.data_handle{1};         
            end
        end
        function detach(self)
            if self.isattached()
                try
                    sharedmatrix('detach',self.shm_key,self.data_handle);
                catch err
                    warning(err.message)
                end
                self.data_handle = []; 
            end
        end
        function free(self)
            % detach the shared memory
            self.detach();
            % and immediatelly delete the memory to avoid pilling it up 
            try sharedmatrix('free',self.shm_key);  end
        end
    end
    
    methods ( Access = 'private' )

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% use MEX file  %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        varargout = sharedmatrix(varargin)
        
        %%%%%%%% other auxiliary functions %%%%%%

        function shmsiz = create_shm(self, data, clone)
            if nargin < 3
                clone = false; % clone == false -> only create the memory space, do not write any data 
            end
            if isnumeric(data) && ~isreal(data)
                %warning('SHM:complex_numbers', 'the sharematrix.c code was modified for matlab newer than R2018b, resulting in very poor performance in complex-valued operations')
                %warning('off', 'SHM:complex_numbers');
                data_handle = {real(data), imag(data)}; 
            else
                data_handle = {data, []}; 
            end
            try
                shmsiz=sharedmatrix('clone',self.shm_key,data_handle, ~clone); 
            catch err
               if strcmpi(err.identifier, 'MATLAB:mex:ErrInvalidMEXFile')
                   % recompile the MEX file 
                   mex @shm/private/sharedmatrix.c -output @shm/private/sharedmatrix
                   shmsiz=sharedmatrix('clone',self.shm_key,data_handle, ~clone); 
               else
                   disp(err)
                   shmsiz = 0; 
               end
            end
                        
            if shmsiz < numel(data)*class2byte(data)
                % it was not possible to allocate enough memory for data 
                self.free_safe();
                disp('Availiable memory:')
                ! free -h
                fprintf('Required memory: %3.2fGB \n', numel(data)*class2byte(data)/1e9)
                ! ipcs -m
                disp(' Use following command to clean all allocated memory ')
                disp('! ipcs -m | cut -d'' '' -f2 | grep ''^[0-9]'' | while read x; do ipcrm -m $x; done')
                error('Shared memory allocation failed, maybe free memory was too low')
            end
            self.isComplex = ~isreal(data);
        end 
        function status = isattached(self)
           % return true is the data are attached 
           status = ~(isempty(self.data_handle) || (iscell(self.data_handle) && isempty(self.data_handle{1}))); 
        end
        function delete(self)
            self.free_safe(); 
        end
        function free_safe(self)
            % detach the shared memory
            self.detach();
            if ~self.protected
                % and immediatelly delete the memory to avoid pilling it up 
                try sharedmatrix('free',self.shm_key);  end
                %try sharedmatrix('free',self.shm_key+1);  end
            end
        end
    end
end

function out = class2byte(in)
    numclass = {'double'; 'single'; 'int8'; 'int16'; 'int32'; 'int64'; 'uint8'; 'uint16'; 'uint32'; 'uint64'};
    numbytes = [NaN;8;4;1;2;4;8;1;2;4;8];

    [~,loc]  = ismember(class(in),numclass);
    out   = numbytes(loc+1);
    if ~isreal(in)
       out = out * 2;  
    end
end
