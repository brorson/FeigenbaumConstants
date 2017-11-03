# A high precision calculation of Feigenbaum's alpha using Julia

Stuart Brorson, Department of Mathematics, Northeastern University, November 2017.

## Background

The theory of dynamical systems has been an area of active research
amongst physicists and mathematicians for several decades.  Perhaps
one of the most interesting objects of study in dynamical systems is
the phenomenon of chaos.  Chaos
(https://en.wikipedia.org/wiki/Chaos_theory) is a phenomenon observed
in some deterministic systems where the system's time-dependent
trajectory neither diverges, nor converges to a limit point, nor moves
periodically.    Rather, the system moves on an orbit which looks
random.  

A simple model system which evidences chaos is the so-called logistic
map (https://en.wikipedia.org/wiki/Logistic_map),

where $\lambda$ is a parameter allowed to vary between 0 4, and xn is the
dynamical variable which varies between 0, 1.  The idea is to choose a
value of lambda, then start with an arbitrary x0, and insert it into
the equation to get a new value, x1.  Then insert x1 into the equation
to get an update, x2.  Repeat this process for large numbers of
iterations and for many values of lambda  the value xn converges to a
value (or set of values).  This value is called the “fixed point” of
the iteration, typically written x*.  The value of x*  depends upon
lambda. 

A plot of the fixed point vs. Lambda is shown below.  What's
interesting are different types of behavior obtained for different
values of lambda.  For lambda < 3, x converges to a single value.  But
starting at lambda = 3, a new behavior emerges:  instead of converging
to one value, the variable x hops from one value to a second one, then
back.  For example, when lambda = 3.2, x hops between approximately
0.7995 and 0.5130.  This is called a “period two” orbit, since x
visits two values, alternating with each iteration.  Moreover,
increasing lambda leads to a situation where x visits 4 values, then 8
then 16, and eventually at a particular value of lambda, the  period
becomes infinity, meaning that x has no periodic orbit.  In this case,
x wanders around the unit interval quasi-randomly – the behavior we
call chaos.

This transition from the period-one orbit to chaos is called the
“period-doubling” route to chaos.  It attracted the attention of many
physicists in the 1970s and 80s since it offered the promise of
illustrating a mechanism by which we could understand the development
of complicated chaotic phenomena such as turbulence in fluids.
Indeed, several physical systems were identified which evidenced a
period-doubling route to chaos, including Duffing's oscillator
[Holmes, P., Whitley, D. On the attracting set for Duffing’s
equation. Physica D, 111–123 (1983)], a dripping faucet[A. D’Innocenzo
and L. Renna, Modeling leaky faucet dynamics, Phys. Rev. E, 55, 6676
(1997).], and some simple electronic circuits["Period Doubling and
Chaotic Behavior in a Driven Anharmonic Oscillator", Paul
Lindsay,. Physical Review Letters, 47 1349 (1981).].  Unfortunately,
the larger ambition of finally getting a grasp on turbulence by
studying the logistic equation did not pan out.  Nonetheless, some
very interesting mathematics was discovered in the process.

In the late 1970s, Mitchell Feigenbaum, a mathematician at Los Alamos
research laboratory, was playing around with the logistic equation
using a hand calculator.  He noticed some interesting numbers
characterize the period doubling transition to chaos in the logistic
map.  He observed: 