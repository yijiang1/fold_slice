%ASTRA_GPU_WRAPPER wrapper around the astra toolkit, it automatically
%recompiles the wrapper if some problems with the MEX file are detected 
%
% varargout = ASTRA_GPU_wrapper(direction, input_array, cfg, vectors,output_array,varargin)
% 
%
% ** direction     'fp' or 'bp' - forward or backward projection operator 
% ** input_array   either projected volume or backprojected projection array 
% ** cfg           cfg structed created by ASTRA_initialize     
% ** vectors       projection geometry created by ASTRA_initialize 
%
% optional:
% ** output_array         either reconstructed volume or projected array
% ** deformation_fields   3x2 or 3x1 cell array if deformation vector fields for nonrigid deformation tomography 
% 
%
% returns:
% ++ output    resulting reconstruction, if output_array ~= [], result will
%               be written to output_array directly to avoid memory allocation 



function varargout = ASTRA_GPU_wrapper(direction, input_array, cfg, vectors,varargin)
    varargout = cell(nargout,1); 
    assert(ismember(direction, {'fp', 'bp'}), 'Wrong option')
    
    try
        % call mex function 
        [varargout{:}] = ASTRA_GPU_wrapper(direction, input_array, cfg, vectors,varargin{:});
    catch err 
        warning(err.identifier,  'ASTRA wrapper returned the following error: %s', err.message)
        if any(strcmp(err.identifier, { 'MATLAB:UndefinedFunction','MATLAB:mex:ErrInvalidMEXFile'}))
            path = replace(mfilename('fullpath'), mfilename, ''); 
            utils.verbose(0, 'Trying to recompile the MEX function ... ')
            
            mexcuda('-outdir',fullfile(path, 'private'), ...
                                        fullfile(path, 'ASTRA_GPU_wrapper/ASTRA_GPU_wrapper.cu'), ...
                                        fullfile(path, 'ASTRA_GPU_wrapper/util3d.cu'),  ... 
                                        fullfile(path, 'ASTRA_GPU_wrapper/par3d_fp.cu'), ...
                                        fullfile(path, 'ASTRA_GPU_wrapper/par3d_bp.cu'));

            [varargout{:}] = ASTRA_GPU_wrapper(direction, input_array, cfg, vectors,varargin{:});
        else
           utils.report_GPU_usage
           rethrow(err) 
        end
    end
end
