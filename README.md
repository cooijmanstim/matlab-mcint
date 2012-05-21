matlab-mcint
============

MEX interface for calling the GNU Scientific Library's Monte Carlo integration functions from MATLAB.

The interface is:

    [y e] = mcint(algorithm, dim, A, B, f, calls, 'PropertyName', PropertyValue);

where algorithm is 'plain', 'miser' or 'vegas'; dim is the dimensionality of the domain of integration; A and B are vectors specifying, for each dimension, the lower and upper bounds of said domain; f is a handle to a function that takes a column vector and returns a scalar; calls is the number of function evaluations; y is the result of the integral; and e is the estimated error.

For example:

    >> [y e] = mcint('vegas', 3, -[0 1 2], [3 4 5], @(x) x'*x, 1000)
    y = 1.4364e+03
    e = 2.6920

which, you know, could easily be correct.

The extra properties that can be passed are listed in the GSL documentation pages for the three algorithms: http://www.gnu.org/software/gsl/manual/html_node/Monte-Carlo-Integration.html .

To build it, you will need to run "mex -setup" in MATLAB and then figure out how to get it to pass the right flags to your compiler.  For me (GCC on 32-bit linux on x86), I had to modify ~/.matlab/R2010a/mexopts.sh and modify the glnx86 section, adding -std=c99 to CFLAGS and -lgsl -lgslcblas to CLIBS.  Then make sure the C file is in your path and run "mex mcint.c".

This was written by Tim Cooijmans in 2012, and released into the public domain.  No warranty.

