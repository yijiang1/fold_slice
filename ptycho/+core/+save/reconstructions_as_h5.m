%RECONSTRUCTIONS_AS_MAT
% Save reconstruction into a h5 file 
% 
% ** p              p structure
% ** final          boolean; false for intermediate saving routines
%
% returns:
% ++ p          p structure
%
% see also: core.save.save_results

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

function p = reconstructions_as_h5(p, final)
    import utils.verbose
    import utils.relative_path
    import io.HDF.*


    if final
        p.plot.extratitlestring = sprintf(' (%dx%d) - Final', p.asize(2), p.asize(1));
    end

    % check if last engine was c_solver
    if strcmpi(p.engines{p.current_engine_id}.name, 'c_solver')
        append2file = true;
    else
        append2file = false;
    end

    s = core.save.extract4saving(p, append2file);


    for ii=1:p.numscans
        if p.share_object
            obnum = 1;
        else
            obnum = ii;
        end
        if p.share_probe
            prnum = 1;
        else
            prnum = ii;
        end
        s.reconstruction.object = ['int_soft:/reconstruction/p/objects/object_' num2str(obnum-1)];
        s.reconstruction.probes = ['int_soft:/reconstruction/p/probes/probe_' num2str(prnum-1)];
        s.reconstruction.Attributes.obnum = obnum;
        s.reconstruction.Attributes.prnum = prnum;
        if isfield(p, 'recon_filename')
            filename_with_path = p.recon_filename{ii};
            if p.queue.isreplica
                [~, fname, ext] = fileparts(p.recon_filename{ii});
                filename_with_path = fullfile(p.save_path{ii}, [fname ext]);
            end
        else
            recons_filename = sprintf('%s_recons.%s',p.run_name, p.save.output_file);
            filename_with_path = fullfile(p.save_path{ii}, recons_filename);

            if exist(filename_with_path, 'file')
                verbose(3,'File %s exists!', filename_with_path);
                alt_filename = filename_with_path;
                [~, fbase,f2] = fileparts(filename_with_path);
                append_number = 0;
                while exist(alt_filename, 'file')
                    f1 = sprintf('%s_%02d', fbase, append_number);
                    alt_filename = fullfile(p.save_path{ii}, [f1 f2]);
                    append_number = append_number + 1;
                end
                filename_with_path = alt_filename;
            end
        end

        if ii==1
            % If the last engine was c_solver, we can use the already existing
            % h5 file.
            if append2file
                movefile(p.recon_filename_c, filename_with_path);
            end
            root_file = filename_with_path;
            s.measurement.data = ['ext:' relative_path(filename_with_path, [p.prepare_data_path p.prepare_data_filename]) ':/'];

        else
            hdf5_cp_file(relative_path(filename_with_path, root_file), filename_with_path, 'groups', {'/measurement/data'; '/measurement/meta_all'; '/reconstruction/p'});
        end

        s.measurement.meta = ['int_soft:/measurement/meta_all/meta_all_' num2str(ii-1)];
        save2hdf5(filename_with_path, s, 'comp', p.io.file_compression);

        try
            fsz = dir(filename_with_path);
            verbose(0, 'Reconstructed scan S%05d: %s', p.scan_number(ii), filename_with_path)
            verbose(2, 'File size: %0.4f MB', fsz.bytes/1e6)
        catch
            verbose(0, 'Saved reconstruction to file %s.', filename_with_path);
        end

        s = [];

    end
        
end