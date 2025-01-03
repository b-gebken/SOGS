% This folder contains all the scripts required for reproducing the
% benchmark in Section 5.1 of our article. The code for the other methods
% is not included and has to be downloaded and placed into respective
% folder manually. (In case of LMBM, it also has to be compiled.) To this
% end, see:
%   Gradsamp (gradsamp): https://cs.nyu.edu/~overton/papers/gradsamp/alg/
%   DGS (DGS-main):      https://github.com/b-gebken/DGS
%   HANSO (hanso3_0):    https://cs.nyu.edu/~overton/software/hanso/
%   SLQPGS (SLQPGS):     https://github.com/frankecurtis/SLQPGS
%   LMBM (LMBM):         https://napsu.karmitsa.fi/lmbm/
% For Gradsamp, HANSO and LMBM, the file 'my_fun' in the respective
% folders have to be kept. For SLQPGS, the file "my_slqpgs_problem" in the
% "src" subfolder has to be kept. (Furthermore, for Gradsamp, the function
% "isinteger" has to be renamed due to a conflict with a MATLAB built-in
% function.)
% 
% "script_benchmark.m" applies all six solvers to the 20 test problems.
% (Note that since this includes all points at which subgradients were
% evaluated, this file will be ~1.2GB in size.)
%
% "script_conv_cutoff_table" produces Table 1 and Table 2. Since our
% artificial termination criterion requires the objective values in all
% points at which subgradients were evaluated, this takes some time.
%
% "script_conv_cutoff_runtime_table.m" produces Table 3. (Here, the cutoff
% threshold is only used to check for convergence.)