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

If you're lucky, you won't need to build it.  There's a binary for 32-bit x86 linux, mcint.mexglc, and a very hard-earned binary for 64-bit x86 windows, mcint.mexw64.

To build it, you will need to run "mex -setup" in MATLAB and then figure out how to get it to pass the right flags to your compiler.  For me (GCC on 32-bit linux on x86), I had to modify ~/.matlab/R2010a/mexopts.sh and modify the glnx86 section, adding -std=c99 to CFLAGS and -lgsl -lgslcblas to CLIBS.  Then make sure the C file is in your path and run "mex mcint.c".

Building on Windows was hell for me; the successful try roughly went like this: install the Microsoft Windows SDK 7.1 package as recommended by MathWorks at http://www.matchworks.com/support/compilers/R2011b/win64.html , start MATLAB and run "mex -setup" in it, select the SDK compiler.  Modify the mexopts.bat script in the directory returned by the MATLAB prefdir command, specifying the right include and library directories and the libraries to link with (libgsl.a, libgslcblas.a).  Run "mex mcint.c" and pray.

If you don't have suitable GSL binaries, here's how to build them.  Download Cygwin, install the packages starting with mingw64-x86_64- (you don't need all of them, but it won't hurt), as well as the package make.  Download the GSL source, add /usr/x86_64-w64-mingw32/bin to your PATH, cd to the GSL source directory, run CC=x86_64-w64-mingw32-gcc ./configure; make.  Pray.  Now libgsl.a should be in .libs/ and libgslcblas.a should be in cblas/.libs.

This was written by Tim Cooijmans in 2012, and released into the public domain.  No warranty.

