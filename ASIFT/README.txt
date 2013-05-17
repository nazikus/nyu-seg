Demo software for Affine-SIFT (ASIFT) image matching
-------------------------------------------------------------------------
-------------------------------------------------------------------------
Jean-Michel Morel (Jean-Michel Morel <morel@cmla.ens-cachan.fr>)
Guoshen Yu (yu@cmap.polytechnique.fr)

Version 2.2, April. 10, 2010

This directory contains the Windows software for ASIFT, a fully affine 
invariant image matching algorithm. 

Software installation and usage are detailed in this manual. If you have 
any problem using the this program, please contact Guoshen Yu 
yu@cmap.polytechnique.fr
				   
For more information about ASIFT, please see the web page at 
http://www.ipol.im/pub/algo/my_affine_sift/.
The source code is available on the web page. You can also try ASIFT using 
the online demo. The online demo allows testing ASIFT with your own images 
without installing the program. 

If you use the ASIFT software, please cite the following paper: 
J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 

-------------------------------------------------------------------------
- The provided Windows executable demo_ASIFT.exe has been compiled by the Intel C++ 
compiler on 32-bit Windows. It is executable on both 32-bit and 64-bit Windows, 
although it is not optimized for the latter. 

- The executable has not been extensively tested. If you have any problem using it,
please contact Guoshen Yu yu@cmap.polytechnique.fr.


-------------------------------------------------------------------------

Usage:
1. Install the program.
Double click the file demo_ASIFTsetp.exe. A small library distributed by Microsoft 
(Microsoft Visual C++ 2010 Redistributable Package) will be installed to your PC. 
The ASIFT software will be installed to C:\Program Files\demo_ASIFT

2. Run ASIFT. 
Run a Dos command prompt (you find it in Start > All Programs > Accessories 

> Command Prompt)
- Go to the ASIFT directory by typing
cd C:\Program Files\demo_ASIFT.
- Type demo_ASIFT for syntax.
- Example: 
demo_ASIFT adam1.png adam2.png imgOutVert.png imgOutHori.png matchings.txt 
keys1.txt keys2.txt

You can of course move the ASIFT directory C:\Program Files\demo_ASIFT to 
wherever that is more convenient. 

-------------------------------------------------------------------------
Troubleshooting 
1. If you are able to run the program but the results are not written to 
the output files, check and make sure that you have the write file permission 
in the demo_ASIFT directory.

2. Microsoft Visual Studio is NOT required to run the program. However, 
in case you cannot run the program (for example some library or dll is missing), 
you may try installing a Microsoft Visual Studio and then running again 
the program. 

-------------------------------------------------------------------------
MATLAB INTERFACE (OPTIONAL)
Run ASIFT via Matlab: Open test_demo_ASIFT.m in Matlab and execute the script. 
The Matlab interface reads most standard image formats.
-------------------------------------------------------------------------
Licensing conditions

This software is being made available for research purposes only. See the
file LICENSE in this directory for conditions of use.
