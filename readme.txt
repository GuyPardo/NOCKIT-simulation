written by Guy 2020_09_25

there are 3 currently usable files in the linear simulation
1. SolveArrayScriptV3 -  solves the problem for a single frequency, calculates the voltage and current along the lines and plots the power propagation.
2. SolveArrayFunV3 - solves the same kind of problem, without reconstruction of physical quantities and without plotting. since it's a function it is more comfortable to use for parameter scans.
3. freq_scan_V3 - uses SolveArrayFunV3 in a loop on frequencies, and plots transmittance and reflectance as a function of frequency 

I think the code is pretty well documented inside, but please let me know if you have questions (or if you find bugs...)

Guy.
guy.pardo@mail.huji.ac.il