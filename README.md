This repo holds the Julia program used to compute Feigenbaum's alpha
to thousands of digits.  This program runs with no warnings in Julia
0.5.  It also runs under Julia 0.6 with a deprecation warning which
you can ignore (for now).

To run the program, execute the function
iterate_alpha().  The results are dumped to the screen, and are also
appended to the file called alpha.dat.

The algorithm is based on solving the self-similarity equation

          g(z, a) = alpha*g(g(z/alpha, a), a)

to get the function g itself.  The computation works like this:

1.  The function g(z) is expanded in a power series:

         g(z) = 1 + a1*x^2 + a2*x^4 + a3*x6 + ....

2.  The function is sampled at a bunch of grid points zi between zero
and one.  

3.  Newton's method is used to compute the ai coefficients: 
We take the gradient of the function at each point zi.  The
gradient is taken w.r.t. the expansion coefficients ai.
The gradients are concatenated to build the Jacobian.  The
function values are also concatenated to build a vector of function
values.  Newton's method finds an approximation to the coefficents ai
by doing a linear solve using the Jacobian and function values the
usual way.  I don't need to do any backtracking or other
special tricks. 

4.  After Newton's method finds the coefficients a1 to the desired
number of digits, alpha is computed by alpha = g(1).

5.  The coefficients are returned to an outer loop (iterate_alpha())
which asks for increasingly more digits.  The coefficients are from
the last run of Newtwon's method are used to seed the next run.  This
helps keep the algorithm from wandering away from the desired
solution, and also provides a series of increasingly precise alpha
values which can be analyzed for convergence.

I have performed some timing measurements, and it (empircally) looks
as if the algorithm asymptotes to O(n^4).  This makes sense because
the linear solve in Newton's method is O(n^3), where n is the number
of sample points zi.  The underlying BigFloat stuff, MPFR, scales as
O(n) with number of digits n.  Since convergence requires the number
of sample points is proportional to number of digits requested, the
overall complexity is O(n^4).  

In very rough numbers, I find that getting around 1000 digits takes on
the order of a day using a ca. 2000 HP laptop.  Note that I have made
no attempt to time individual parts of the code, nor to optimize the
code. 

The solution method is described in a number of publications,
including:

* "Introduction to universality and renormalization group techniques", https://arxiv.org/abs/1210.2262

* "A simpler derivation of Feigenbaum's renormalization group equation
for the period doubling sequence", http://chaosbook.org/library/Copper98Feig.pdf

---

A subdirectory of results is included here.  The files are:

alpha_me_*.txt -- My computed results from Oct 2017.  The numeric
                  suffixes are the number of digits requested by
		  iterate_alpha(), not necessariy the number of
		  correct digits.  

alpha_oeis.txt -- All 1018 digits of alpha available from 
                  http://www.plouffe.fr/simon/constants/feigenbaum.txt

verify.py -- a python program which compares each digit from two input
             files and reports the number of digits of agreement.  Use
	     this to analyze convergence of the computed alphas.


Last update: Stuart Brorson  10.14.2017.


# Instructions:
```
# Run one time to get the Julia ForwardDiff package:
Pkg.add("ForwardDiff")
```

```
# To run the program:
include("compute_alpha.jl")
iterate_alpha();
```
