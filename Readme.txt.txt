I.  High Level programs (in the "Main_Calls" folder) call the FDTR and TDTR calculations to accomplish specific tasks.
*THESE ARE THE PROGRAMS USERS NORMALLY EDIT/RUN FOR THEIR SPECIFIC NEEDS.  IF ONE DOESN'T EXIST FOR YOUR TASK, YOU SHOULD WRITE ONE.*
These are the included "Main_Calls" right now:

	A.  Load/Fit standard TDTR data ( ratio vs td)
	B.  Estimate Error bars
	C.  Calculate FWHM of beam-offset TDTR experiments -> ("TDTR_Ani3D_FWHM_Map.m")
	D.  Calculate sensitivity plots for
		i)  standard TDTR, ratio vs. td
		ii)  beam offset TDTR, FWHM vs. freq

%%%%%%%%%%%%%%%%%
How the scripts/functions are structured:

There are several programs that, together, allow you to simulate TDTR or FDTR with arbitrary spot shape and arbitrary thermal conductivity tensor.  

NOTE:  The programs utilize multicore and multiCPU capability when software and hardware permit.  This distribution is particularly enhanced by the parallel computing toolbox, although it will use multicore regardless.
%%%%%%%%%%%%%%%%%

II.  Calculating the Frequency Response (i.e FDTR)

At the core of the program is the ability to calculate the frequency domain response (i.e. the "sensed temperature" of a temporaly contant intensity sensing beam when a sinusoidally oscillating "pump" beam heats the surface) for arbitrary spot shapes and thermal conductivity tensors.  For TDTR, the signal is simulated as a near-infinite sum over all the relavent frequency responses.  

The required M-files are:
"rombint2D.m"
"trap2D.m"
"TDTR_3DAni_getTintegral_wOffset_v2.m"
"TDTR_3DAni_getG_savespace.m"

To calculate the frequency response:  

- A 2D integration is performed using "rombint2D.m" which uses Romberg integration in 2D. The integration is over a temperature "Kernal" which is basically the Fourier transforms of the beam intensities(which we denote S(xi,eta) and P(xi,eta)), and the Fourier domain Green's function for the system thermal response (which we denote G(xi,eta)).

- "trap2D.m" is the 2D implementation of the trapezoidal rule and is called by "rombint2D.m"

- The Kernal is evaluated using the function call:

[Kernal] = TDTR_3DAni_getTintegral_wOffset_v2(xi,eta,f,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset)

- "TDTR_3DAni_getTintegral_wOffset_v2.m" calls the function "TDTR_3DAni_getG_savespace.m"  The "savespace" refers to the fact the the program reuses some of the matrices in order to save memory.  This was a major issue in early version of the program, but may not be necessary now.  As a result of the "savespace" the script is difficult to read, (consider revising in a future version)  

[G] = TDTR_3DAni_getG_savespace(xi,eta,f,k,C,h)

%%%%%%%%%%%%%%%%
III.  Calculating the TDTR Signal

The required M-files are:
(Every file required for FDTR) +
"TDTR_3DAni_getVout_wOffset_v2.m"

To calculate the TDTR signal:

- the call to "TDTR_3DAni_getVout_wOffset_v2.m" computes the TDTR signal by summing over all relavant frequency domain responses.  The call format is:

[Vout] = TDTR_3DAni_getVout_wOffset_v2(tdelay,k,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset)

where "tdelay" is a vector of time delays.  If TCR=1, then "Vout" is the Out-of-phase Temperature, if TCR is the thermoreflectance coefficient, then "Vout" is the out-of-phase reflectivity signal. 
%%%%%%%%%%%%%%%%