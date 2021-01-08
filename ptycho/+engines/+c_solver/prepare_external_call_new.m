% Call the external reconstruction program

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
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
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
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
    
function [p, fdb] = prepare_external_call(p, fdb)
import utils.check_cpu_load
import utils.verbose

    if isempty(p.reconstruction_program)
        
        if ispc
            error('C_solver engine is not supported for Windows systems')
        end
        if p.single_prec
            prec = 'singlePrec';
        else
            prec = 'doublePrec';
        end

        if isempty(p.caller_suffix)
            suffix = '';
        else
            suffix = ['_' p.caller_suffix];
        end
        
        if isfield(p, 'ra_reservation') && ~isempty(p.ra_reservation)
            reservation_flag = [' --reservation=' p.ra_reservation];
        else
            reservation_flag = ' ';
        end        
        
        solver_path = fullfile(p.ptycho_package_path,'/+engines/+c_solver/');
        
        if ~p.use_gpu
            num_threads = sprintf('OMP_NUM_THREADS=%d', p.threads);
            
            [status, hostname] = system('hostname');
            fdb.status = core.engine_status(status);
            host_pre = strsplit(hostname, '-');
            
            if ~isfield(p, 'hybrid')
                p.hybrid = 1;
            end
            if p.number_iterations+p.opt_iter > 5000
                slurm_partition='week';
            else
                slurm_partition='day';
            end
            switch host_pre{1}
                case 'ra'
                    if p.ra_nodes == 0
                        multiproc = 'OMP';
                        p.reconstruction_program = [num_threads ' ', solver_path 'ptycho_' prec '_' multiproc suffix];
                    elseif p.ra_nodes==1
                        multiproc = 'OMP';
                        p.reconstruction_program = ['OMP_PROC_BIND=spread ' num_threads ' salloc -p ' slurm_partition reservation_flag ' --ntasks-per-node=1 --job-name=ptycho_recons -N ' num2str(p.ra_nodes) ' $(which mpirun) -x OMP_NUM_THREADS -x OMP_PROC_BIND --bind-to none --map-by ppr:1:node  ', solver_path 'ptycho_' prec '_' multiproc suffix];
                    else
                        if p.hybrid
                            multiproc = 'hybrid';
                            p.reconstruction_program = ['OMP_PROC_BIND=spread ' num_threads ' salloc -p ' slurm_partition reservation_flag ' --ntasks-per-node=1 --job-name=ptycho_recons -N ' num2str(p.ra_nodes) ' $(which mpirun) -x OMP_NUM_THREADS -x OMP_PROC_BIND --bind-to none --map-by core=' num2str(p.threads) ' --map-by ppr:1:node  ', solver_path 'ptycho_' prec '_' multiproc suffix];
                        else
                            multiproc = 'MPI';
                            p.reconstruction_program = ['salloc -p ' slurm_partition reservation_flag ' --ntasks-per-node=1 -N ' num2str(p.ra_nodes) ' $(which mpirun) -x --map-by core=' num2str(p.threads) ' --map-by ppr:1:node  ', solver_path 'ptycho_' prec '_' multiproc suffix];
                        end
                    end
                case 'x12sa'
                    if isempty(p.beamline_nodes)
                        multiproc = 'OMP';
                        p.reconstruction_program = [num_threads ' ', solver_path 'ptycho_' prec '_' multiproc suffix];
                        
                    else
                        %                     if p.check_cpu_load
                        %                         fprintf('Checking CPU load...\n');
                        %                         cn_usage = check_cpu_load(p.beamline_nodes);
                        %                         if any(cn_usage>=15)
                        %                             fprintf('Usage: \n');
                        %                             check_cpu_load;
                        %                             pause(10);
                        %                             fprintf('Okay, I will try to squeeze in...\n');
                        %                         end
                        %                     end
                        
                        hosts = join(cellstr(p.beamline_nodes), ',');
                        Nhosts = length(p.beamline_nodes);
                        
                        if p.hybrid
                            multiproc = 'hybrid';
                            p.reconstruction_program = [num_threads ' $(which mpirun) --bind-to none -x OMP_NUM_THREADS -x LD_LIBRARY_PATH -H ' hosts{1} ' -np ' num2str(Nhosts) ' ', solver_path 'ptycho_' prec '_' multiproc suffix];
                        else
                            multiproc = 'MPI';
                            p.reconstruction_program = ['$(which mpirun) -x LD_LIBRARY_PATH -H ' hosts{1} ' -np ' num2str(Nhosts) ' ', solver_path '/ptycho_' prec '_' multiproc suffix];
                        end
                        
                    end
                otherwise
                    multiproc = 'OMP';
                    p.reconstruction_program = [num_threads ' ', solver_path 'ptycho_' prec '_' multiproc suffix];
                    verbose(2, 'Unknown host. Please make sure that the correct libraries are loaded or run it on the DaaS / X12SA beamline nodes.');
            end
        else
            if p.num_gpus==1
                
                p.reconstruction_program = [solver_path '/ptycho_' prec '_singleGPU' suffix ' --cuda_device ' num2str(p.gpu_id-1)];
                disp(p.reconstruction_program )
            else
                p.reconstruction_program = ['$(which mpirun) -np ' num2str(p.num_gpus) ' ' solver_path '/ptycho_' prec '_multiGPU' suffix ' --cuda_device ' num2str(p.gpu_id-1)];
            end
        end
        
    end
end

