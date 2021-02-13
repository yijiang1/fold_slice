%MATLAB_POS calculate scan parameters based on the values set in the
%template

function [ p ] = matlab_pos( p )

for ii = 1:p.numscans
    positions_real = zeros(0,2); 
    switch p.scan.type
        case 'raster'
            scan_order_x = 1:p.scan.nx;
            scan_order_y = 1:p.scan.ny;
            % Added by ZC: flip positions similar to eng.custom_data_flip in GPU engines
            if isfield(p.scan, 'custom_flip') && any(p.scan.custom_flip) 
                warning('Applying custom scan flip: %i %i %i ', p.scan.custom_flip(1), p.scan.custom_flip(2), p.scan.custom_flip(3))
                if p.scan.custom_flip(1)
                    scan_order_x = fliplr(scan_order_x); 
                end
                if p.scan.custom_flip(2)
                    scan_order_y = fliplr(scan_order_y);
                end
            end

            for iy=1:length(scan_order_y) %modified by YJ. seems odd to begin with 0...
                for ix=1:length(scan_order_x)
                    xy = [scan_order_y(iy) * p.scan.step_size_y, scan_order_x(ix) * p.scan.step_size_x] + ...
                            randn(1,2).*p.scan.step_randn_offset.*[ p.scan.step_size_y,  p.scan.step_size_x];
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
            
            if isfield(p.scan, 'custom_flip') && p.scan.custom_flip(3) % switch x/y by ZC
            	positions_real=fliplr(positions_real);
            end
            
        case 'round'
            dr = (p.scan.radius_out - p.scan.radius_in)/ p.scan.nr;
            for ir=1:p.scan.nr+1
                rr = p.scan.radius_in + ir*dr;
                dth = 2*pi / (p.scan.nth*ir);
                for ith=0:p.scan.nth*ir-1
                    th = ith*dth;
                    xy = rr * [sin(th), cos(th)];
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
            
        case 'round_roi'
            rmax = sqrt((p.scan.lx/2)^2 + (p.scan.ly/2)^2);
            nr = 1 + floor(rmax/p.scan.dr);
            for ir=1:nr+1
                rr = ir*p.scan.dr;
                dth = 2*pi / (p.scan.nth*ir);
                for ith=0:p.scan.nth*ir-1
                    th = ith*dth;
                    xy = rr * [sin(th), cos(th)];
                    if( abs(xy(1)) >= p.scan.ly/2 || (abs(xy(2)) > p.scan.lx/2) )
                        continue
                    end
                    positions_real(end+1,:) = xy; %#ok<AGROW>
                end
            end
            
        case 'fermat'
            % this should be changed to have the same variable
            % conventions as in its spec implementation
            phi=2*pi*((1+sqrt(5))/2.) + p.scan.b*pi;
            start = 1;
            if ~isempty(p.scan.lx)
                for ir=start:p.scan.n_max
                    r=p.scan.step*0.57*sqrt(ir);
                    if abs(r*sin(ir*phi))> p.scan.ly/2
                        continue
                    end
                    if abs(r*cos(ir*phi))> p.scan.lx/2
                        continue
                    end
                    xy  = [r*sin(ir*phi)+p.scan.cenxy(1) r*cos(ir*phi)+p.scan.cenxy(2)];
                    positions_real(end+1,:) = xy;
                end
            else
                for ir=start:p.scan.n_max
                    r=p.scan.step*0.57*sqrt(ir);
                    xy  = [r*sin(ir*phi)+p.scan.cenxy(1) r*cos(ir*phi)+p.scan.cenxy(2)];
                    positions_real(end+1,:) = xy;
                end
            end
            
        case 'custom' %for PSI's data
            fn_splt = strsplit(p.scan.custom_positions_source,'.');
            if length(fn_splt)>1
                % file already has an extension
                ext = fn_splt(end);
                if strcmp(ext, 'm')
                    [~, positions_real, ~] = p.scan.custom_positions_source(p);
                elseif strcmp(ext, 'mat')
                    posi = load(p.scan.custom_positions_source, 'pos');
                    positions_real = posi.pos;
                    clear posi;
                else
                    error('File extenstion %s is not supported.', ext)
                end
                
            else
                % file does not have an extension
                if exist([p.scan.custom_positions_source '.m'], 'file')
                    [~, positions_real, ~] = p.scan.custom_positions_source(p);
                elseif exist([p.scan.custom_positions_source '.mat'], 'file')
                    posi = load(p.scan.custom_positions_source, 'pos');
                    positions_real = posi.pos;
                    clear posi;
                else
                    error('Could not find function or data file %s', p.scan.custom_positions_source);
                end
            end
            
        case 'custom_GPU' %added by YJ for customized GPU engines' output
            if ~isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = p.scan.custom_positions_source;
            else
            	error('Position file is not given');
            end
            
            try
                r_output = load(pos_file,'outputs');
                r_p = load(pos_file,'p');
            	ppX = r_output.outputs.probe_positions(:,1)*r_p.p.dx_spec(1);
                ppY = r_output.outputs.probe_positions(:,2)*r_p.p.dx_spec(2);
                ppX = ppX(:);
                ppY = ppY(:);
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY;
                positions_real(:,2) = -ppX;                
            catch
                error('Failed to load positions from %s', pos_file);
            end
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real];
end
    
end

