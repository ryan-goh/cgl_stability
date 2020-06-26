README.txt

%%Software:
For continuation runs used AUTO07p

For numerical spectrum computations used MATLAB2019b


Description of directories and code contained:

/Pulled  --> Numerics for the equation A_t = (1+i \alpha)A_{xx} +\chi(x -ct)A - (1+i\gamma)A|A|^2
    /Parameter1 -->Numerics for parameters \alpha = -0.1, \gamma = -0.2,
    /Parameter2 -->Numerics for parameters \alpha = -0.2, \gamma = -0.9,

        Inside each directory:
            /AUTO_Existence --> Script, FORTRAN files, and constant files to numerically continue front solutions. To run, open AUTO07p in a terminal, and copy and run the uncommented script text. Run this before running spectral code

            spec_wgt.m --> Numerical eigenvalue computation of linearized operator about the front solution, draws solution profiles from /AUTO_Existence data,  Set gpuon = 1  (only in Parameter 1) to convert the sparse linearization matrix into a non-sparse/full matrix and use "eig" function to find eigenvalues using a GPU. Note to do this one must have an NVIDIA gpu with the CUDA toolkit installed.

            plot_spec.m --> Plot figure for the manuscript (as well as others)
                        --> Set makespecvid = 1 (only in Parameter 1) to create a video showing how eigenvalues move as c is varied.

/Pushed
    /AUTO_Existence --> Numerics for the equation A_t = (1+i \alpha)A_{xx} +\chi(x -ct)A + (\mu+i\gamma)A|A|^2 - (1+i*\beta)A|A|^4
    To run, open AUTO07p in a terminal, and copy and run the uncommented script text.

    pushed_spec.m -->  Numerical eigenvalue computation of linearized operator about the front solution, draws solution profiles from /AUTO_Existence data,  Set gpuon = true  to convert the sparse linearization matrix into a non-sparse/full matrix and use "eig" function to find eigenvalues using a GPU. Note to do this one must have an NVIDIA gpu with the CUDA toolkit installed.


    plot_spec.m --> Plots figure for the manuscript (as well as others)
