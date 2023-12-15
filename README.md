# fold_slice

This is Yi Jiang's customized code for X-ray/electron ptychography and tomography/laminography.

The package is built upon the Matlab code developed by the Science IT and the coherent X-ray scattering (CXS) groups at Paul Scherrer Institut, Switzerland:
https://www.psi.ch/en/sls/csaxs/software. Copyright and license issues should follow the agreements (see below) and/or refer to their website.

# Get started
1. Check the official [documentation](https://doi.org/10.1107/S1600576720001776) of the PtychoShelves package to see its requirements. You need the following Matlab toolbox to use all the features: Parallel computing, Curve Fitting, Image processing, Optimization, and Signal processing.

2. For ptychography, try the data preparation and reconstruction scripts in /fold_slice/ptycho/examples to get familiar with the data format and reconstruction parameters.

# Resources
1. The files in /fold_slice/ptycho/notes/ could help you understand the overall code structure. Warning: some notes might be outdated.
2. [Dr. Chia-Hao Lee](https://sites.google.com/view/chiahao-lee/home?pli=1) wrote a great [blog](https://chiahao-blog.super.site/posts/theory-algorithm-and-code-structure-of-ptychoshelves) that details the algorithms and code structure of PtychoShelves.
3. We have a weekly study group to discuss novel computational imaging techniques in electron microscopy. Some tutorial lectures can be found [here](https://anl.box.com/s/f7lk410lf62rnia70fztd5l7n567btyv).

# Major differences from the PtychoShelves package
1. Some data and reconstruction I/O conventions have been changed to accommodate for electron ptychography. See the example scripts for more details.

2. A modified least-squares maximum likelihood multi-slice ptychography algorithm is added as a new engine: GPU_MS

It's based on the GPU engine (written by Michal Odstrcil) with improvements such as multiple probe modes and bug fixes. Usage of the code should include additional citations:

Z. Chen, Y. Jiang, Y. Shao, M. E. Holtz, M. Odstrčil, M. Guizar-Sicairos, I. Hanke, S. Ganschow, D. G. Schlom, D. A. Muller, Electron ptychography achieves atomic-resolution limits set by lattice vibrations. Science 372 (6544), 826-831.

3. We developed an automatic parameter tuning workflow for ptychography using Bayesian optimization with Gaussian processes: https://doi.org/10.1038/s41598-022-16041-5. See the example scripts for more details.

4. A non-exhaustive list of new features in the GPU and GPU_MS engines: 

| Features  | GPU         |  GPU_MS | 
| :---         |     :---:      |  :---: |
| Mixed-states + multislice ptychography  | :heavy_multiplication_x:  | :heavy_check_mark:  |
| Dynamic multislice reconstruction | :heavy_multiplication_x:  | :heavy_check_mark:  |
| Advanced arbitrary-path fly-scan ptychography| :heavy_check_mark: | :heavy_multiplication_x:|
| Multi-scan reconstruction | :heavy_check_mark: | :heavy_check_mark: |
| TV regularization on object phase| :heavy_check_mark: | :heavy_check_mark: |
| Grid artifact removal| :heavy_check_mark: | :heavy_multiplication_x: |
| Automatic parameter selection| :heavy_check_mark: | :heavy_check_mark: |
| Account for detector blur with a Gaussian kernel| :heavy_check_mark: | :heavy_check_mark: |

# Other ptychography software
If you don't own Matlab or want to explore other ptychography software. Here are some public repositories:

Adorym: https://github.com/mdw771/adorym

Ptycho_gui: https://github.com/NSLS-II/ptycho_gui

Ptychodus: https://github.com/AdvancedPhotonSource/ptychodus

PtychoNN: https://github.com/mcherukara/PtychoNN

Ptychopy: https://github.com/kyuepublic/ptychopy

Py4DSTEM: https://github.com/py4dstem/py4DSTEM

PyNX: http://ftp.esrf.fr/pub/scisoft/PyNX/doc/

Tike: https://github.com/tomography/tike


# Academic License Agreement

Source Code

Introduction 

This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").

Terms and Conditions of the LICENSE
1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions hereinafter set out and until termination of this license as set forth below.
2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the LICENSEE's responsibility to ensure its proper use and the correctness of the results.
3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged in the commercial use, application or exploitation of works similar to the PROGRAM.
5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into another computing language:
"Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul Scherrer Institut, Switzerland."

Additionally, any publication using the package, or any translation of the code into another computing language should cite

(for PtychoShelves) K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)


(for difference map) P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379-382 (2008). (doi: 10.1126/science.1158573).

(for maximum likelihood) P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). (doi: 10.1088/1367-2630/14/6/063004).

(for mixed coherent modes) P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68-71 (2013). (doi: 10.1038/nature11806).

(and/or for multislice) E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089-29108 (2016). (doi: 10.1364/OE.24.029089).

6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
@ All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before the courts of Zurich, Switzerland. 
