# A practical max-min fair resource allocation algorithm for rate-splitting multiple access

This is a code package related to the following [paper](https://ieeexplore.ieee.org/document/10304220):

F. Luo and Y. Mao, "A Practical Max-Min Fair Resource Allocation Algorithm for Rate-Splitting Multiple Access," in _IEEE Communications Letters_, vol. 27, no. 12, pp. 3285-3289, Dec. 2023, doi: 10.1109/LCOMM.2023.3329149.

# Content of Code Package

Here is a detailed description of the package: 
- The code in all packages are implemented in MATLAB environment with CVX toolbox assisted. 
- The procedure for Fig. 1:
    - Run 'channel_generate.m' and you will get a 'channel.mat' file. 
    - Copy the 'channel.mat' file to every subfolders corresponding to each method.
    - Run every '\*main.m' MATLAB script and you will get corresponding '\*.mat' files.
    - Copy the '\*.mat' files to corresponding subfolders that is in 'data_processing' folder.
    - Run 'data_processing.m' MATLAB script.
- Fig. 2 of the above paper will be reproduced by running MATLAB script 'optimalT.m'. By changing the variable 'SNRdB', you can reproduce Fig. 2(a-c). 
- Fig. 3 of the above paper will be reproduced by running MATLAB script 'relativeGain.m'.

# Abstract of the Article

This letter introduces a novel resource allocation algorithm for achieving max-min fairness (MMF) in a rate-splitting multiple access (RSMA) empowered multi-antenna broadcast channel. Specifically, we derive the closed-form solution for the optimal allocation of the common rate among users and the power between the common and private streams for a given practical low-complexity beamforming direction design. Numerical results show that the proposed algorithm achieves 90% of the MMF rate on average obtained by the conventional iterative optimization algorithm while only takes an average of 0.1 millisecond computational time, which is three orders of magnitude lower than the conventional algorithm. It is therefore a practical resource allocation algorithm for RSMA.

# License and Referencing
This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

# Acknowledgements

This work has been supported in part by the National Natural Science Foundation of China under Grant 62201347; and in part by Shanghai Sailing Program under Grant 22YF1428400.