% REGULATION_MULTILAYERS try to avoid ambiguity in the multilayer reconstruction by weakly forcing missing cone
% values towards zero 
%
% self = regulation_multilayers(self, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ self        self-like structure with final reconstruction
%

function self = regulation_multilayers(self, par, cache)
    import engines.GPU_MS.GPU_wrapper.*
    %% Added by ZC. Use CPU is object size is too big
    % par.obj_size_limit_on_gpu is maximum object size (in MB) allowed on gpu.
    % Automatically use cpu if exceed the limit.
    Obj_size_limit = par.obj_size_limit_on_gpu / 8 * 2^20 ;
    %%
    Npix = [self.Np_o, par.Nlayers]; % -1, Not last inf layer, by Zhen Chen
    for i = 1:3
        grid{i} = ifftshift((-fix(Npix(i)/2):ceil(Npix(i)/2)-1))'/Npix(i);
        grid{i} = shiftdim(grid{i},1-i);
    end
    % calculate force of regularization based on the idea that DoF = resolution^2/lambda
    W = 1-atan(( par.regularize_layers * abs(grid{3}) ./ sqrt(grid{1}.^2+grid{2}.^2+1e-3)).^2) / (pi/2); 
    relax = 1; 
    alpha = 1; 
    Wa = W.*exp(-alpha*(grid{1}.^2 + grid{2}.^2));
    
    for kk = 1:size(self.object,1)
        %% Added by ZC. use CPU to save GPU memory if object is too big
        if  numel(self.object{kk,1}) * size(self.object,2) > Obj_size_limit && par.use_gpu
            object=cellfun(@gather, self.object(kk,:),'UniformOutput',false);
            obj = cat(3, object{:});
        else
            obj = cat(3, self.object{kk,:});
        end
        %%
        % find correction for amplitude 
        aobj = abs(obj); 
        fobj = fftn(aobj); 
        fobj = fobj .* Wa; 
        aobj_upd = ifftn(fobj); 
        % push towards zero 
        aobj_upd = 1+0.9*(aobj_upd-1); 
        % find correction for phase 
        Wphase = min(1, 10*(cache.illum_sum_0{kk}/cache.MAX_ILLUM(kk))); 
        if numel(obj) > Obj_size_limit && par.use_gpu % Added by ZC to save memory
        	Wphase=gather(Wphase);          
        end
        pobj = math.unwrap2D_fft2(obj,[],0,Wphase,-1);
        fobj = (fftn((pobj))); 
        fobj = fobj .* Wa; 
        pobj_upd = ifftn(fobj); 
        if numel(obj) > Obj_size_limit && par.use_gpu % Added by ZC to save memory
          obj_upd =  regulation_multilayers_kernel(obj, aobj,aobj_upd, pobj, pobj_upd, Wphase, relax); 
        else
          obj_upd =  Gfun(@regulation_multilayers_kernel,obj, aobj,aobj_upd, pobj, pobj_upd, Wphase, relax); 
        end
        for ii = 1:par.Nlayers  % -1, Not last inf layer, by Zhen Chen
          self.object{kk,ii} = Garray(obj_upd(:,:,ii)); 
        end
    end
end
function [obj,corr] =  regulation_multilayers_kernel(obj, aobj,aobj_upd, pobj, pobj_upd, weights, relax)
  aobj_upd = (real(aobj_upd) - aobj); 
  pobj_upd = weights.*(real(pobj_upd) - pobj); 
  corr = (1+relax*aobj_upd) .* exp(1i*relax*pobj_upd); 
  obj = obj .* corr; 
end        
  