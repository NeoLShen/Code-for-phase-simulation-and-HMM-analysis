# Codes for Phase Separation Simulation and HMM Analysis
These are codes for phase separation simulation and HMM analysis that are coded with MATLAB.


The HMM analysis folder contains 7 scripts. The imported data should contain the information of x and y coordinates as "Xm" and "Ym", and information of frame number as "Frame". All three variables should be a vector.

"TwoPhaseTwoStateProc.m" was used for automatically determining the optimized search range in condensed phase and dilute phase, and generating trajectory information that contains [X coordinates, Y coordinates, Frame number, Phase label, Displacement, Mean squared displacement] in variable 'Trajsio'.

The "CreateMask.m" was quoted in "TwoPhaseTwoStateProc.m" first to process the imported data and generate a condensed phase mask based on local density.

"trajasin.m" and "trajasinout.m" were quoted in "TwoPhaseTwoStateProc.m" to process specifically condensed phase trajectories in default settings and all trajectories in both condensed and dilute phase regions simultaneously.

"paraest.m" was quoted in "TwoPhaseTwoStateProc.m", and was used for performing HMM analysis on selected trajectories.


The Phase simulation folder contains 3 scripts.

"phasesimulator1.m" was used for generating simulated molecule motion followed 2-state motion switch model with experimental condensed phase mask in a simulated region with periodic boundaries. The script generates [X,Y,P,M,t] data, contains x/y coordinates, and labels indicate whether the molecule is in condensed phase (P, 1 for in condensed phase and 0 for dilute phase)/motion state (M, 1 for mobile state and 0 for confined state)/time stamp (t).

"phasesimulator2.m" was used for generating simulated molecule motion following a 2-state motion switch model with manually set phase mask in a simulated region with periodic boundaries. The script generates [X,Y,P,M,t] data, contains x/y coordinates, and labels indicate whether the molecule is in condensed phase (P, 1 for in condensed phase and 0 for dilute phase)/motion state (M, 1 for mobile state and 0 for confined state)/time stamp (t). Parallel computation was required, if you don't want to use parallel computation, please replace 'parfor' with 'for'.


If you encounter any technical problems, please contact the author by email: "zshenai@connect.ust.hk" or "shenzeyu@ust.hk"
