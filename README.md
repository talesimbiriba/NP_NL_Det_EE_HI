************************************
************************************

This Package has all programs and third part software to run the simulations presented in [1]. 
All third part software used here has free distribution and noncomercial use licenses. 
The routine setup.m sets up the matlab path for our and third party software and other 
configuration needed. 
All third party software is listed below. 

The use and redistribution of this code is allowed for noncomercial purposes as long as the 
copyrights presented here are replicated. 

Author: Tales Imbiriba.

Date: December 2015. 

Ref.:

[1]  T. Imbiriba, J.C.M. Bermudez, C. Richard, J.-Y. Tourneret, "Nonparametric Detection 
of Nonlinearly Mixed Pixels and Endmember Estimation in Hyperspectral Images", to appear 
in the IEEE Transaction on Image Processing Magazine, 2015.


************************************
************************************

Main Files:

setup.m 
	This program adds all the needed code (ours and third part) to matlabs path. Furthermore
	it also creates runOnAllJavaMonitorCP.m function needed by the java pregression Monitor
	used with 'parfor' to speed up simulations.

allSamplesGPLSRatioNonlinearMixtureDetection.m
	This program call the detectors and generates the ROCs and histograms of the statistics. 

detectAndUnmix_IteratEMEst_SynthData_RecError.m
	This script detect nonlinearly mixed pixels, using different detectors, and unmix the 
	data using FCLS or SK-Hype. It also has the option to estimate the endmembers using the 
	proposed method calling the function iterativeEndmemberEstAndNonlinDetect.m

IndianPines_IterativeMVSEndmemberEstimationAndDetection.m
	This program Unmix the Indian Pines image using the proposed endmember estimation and 
	the detect-then-unmix strategy as described in [1]. 

Cuprite_Simulation.m
	This program executes the simulation with the Cuprite Mining field. 



************************************
************************************


Third part software:


	SK-Hype (Matlab code avaiable at http://cedric-richard.fr/pub.html)

	MVES - Minimum volume enclosing simplex (http://mx.nthu.edu.tw/~tsunghan/Source%20codes.html)

	ParforProgMon version 2  - progress monitor in parfor loops (copyright below)

	CCRLB - Nonlinearity detector considering the PPNMM. Provided by Yoann Altmann. 

	SEDUMI - Matlab Toolbox for Linear, Second order, Semidefinite or mixed	problems. 
		It can handle free variables and complex data as well. (http://sedumi.ie.lehigh.edu)
	
	GPML Version 3.3 - Gaussian Process package. (Copyright (c) by Carl Edward Rasmussen and Hannes
      			   Nickisch, 2013-10-19.) avaiable at http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html



************************************
************************************

**** Copyrights ****** 



ParforProgMonv2: 

Copyright (c) 2011, Willem-Jan de Goeij
Copyright (c) 2009, The MathWorks, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
    * Neither the name of the The MathWorks, Inc. nor the names 
      of its contributors may be used to endorse or promote products derived 
      from this software without specific prior written permission.
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.




GPML: 


 GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION Toolbox version 3.6
    for GNU Octave 3.2.x and Matlab 7.x

The code is released under the FreeBSD License.

Copyright (c) 2005-2015 Carl Edward Rasmussen & Hannes Nickisch. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY CARL EDWARD RASMUSSEN & HANNES NICKISCH ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CARL EDWARD RASMUSSEN & HANNES NICKISCH OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of Carl Edward Rasmussen & Hannes Nickisch.

The code and associated documentation is available from
http://gaussianprocess.org/gpml/code.



