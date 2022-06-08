# HologramsGenration

The purpose of this library is to generate optical tweezers hologram for a spatial light modulator (SLM).
The optical tweezers are a set of points on the image plane of a microscope. The SLM is located on the Fourier plane of the microscope.

The fastphase algorithm is a naive Fourier transform of the given points.

The dsphase algorithm takes fastphase as its starting points, and adds optimization according to the Direct Search algorithm. The details of direct search can be found here: 
-  David G. Grier and Yael Roichman, "Holographic optical trapping," Appl. Opt. 45, 880-887 (2006) 
-  Marco Polin, Kosta Ladavac, Sang-Hyuk Lee, Yael Roichman, and David G. Grier, "Optimized holographic optical traps," Opt. Express 13, 5831-5845 (2005) 
