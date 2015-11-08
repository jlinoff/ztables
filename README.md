# ztables
Ever wonder how Standard Normal and Student t-Distribution z-tables are generated?

This program shows how to generate z-tables for Standard Normal
Distributions and Student-t Distributions. It also shows how to
generate the z-values for specific probabilities from those
distributions without using any special libraries.

The goal is to show how to solve the probability density functions
(PDFs) and use those solutions to approximate the area under a curve
in order to obtain the z-value or, in the case of the z-value lookup,
show how to do a binary search the probabilities to find the desired
z-value.

It is meant to be a learning tool for those that are interested in
how the tables are generated.

I explicitly chose to use the trapezoidal rule to approximate the
definite integral representing the area under the curve because it is
simple, reasonably fast and accurate enough. You could, of course,
choose any other approximation method (like Simpson's Rule).

I also chose the Lanczos approximation for solving the gamma function
for the t-distribution PDF. The algorithm for the approximation was
obtained from page 214 of "Numerical Recipes in C, 2nd Edition".

Here is an example that shows you how to generate the z-table that
contains the probabilities for the standard normal distribution using
the default parameters. This is the table that you often see in text
books.

```
$ ztables.py -s

z-Table for Standard Normal Distribution (10,000)

  z     .00     .01     .02     .03     .04     .05     .06     .07     .08     .09
=====  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
 -3.4  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0002
 -3.3  0.0005  0.0005  0.0005  0.0004  0.0004  0.0004  0.0004  0.0004  0.0004  0.0003
 -3.2  0.0007  0.0007  0.0006  0.0006  0.0006  0.0006  0.0006  0.0005  0.0005  0.0005
 -3.1  0.0010  0.0009  0.0009  0.0009  0.0008  0.0008  0.0008  0.0008  0.0007  0.0007
 -3.0  0.0013  0.0013  0.0013  0.0012  0.0012  0.0011  0.0011  0.0011  0.0010  0.0010
  .
  .

  3.0  0.9987  0.9987  0.9987  0.9988  0.9988  0.9989  0.9989  0.9989  0.9990  0.9990
  3.1  0.9990  0.9991  0.9991  0.9991  0.9992  0.9992  0.9992  0.9992  0.9993  0.9993
  3.2  0.9993  0.9993  0.9994  0.9994  0.9994  0.9994  0.9994  0.9995  0.9995  0.9995
  3.3  0.9995  0.9995  0.9995  0.9996  0.9996  0.9996  0.9996  0.9996  0.9996  0.9997
  3.4  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9998
```

Here is an example that shows how to determine the z values for 95%,
98% and 99% probabilities using the Standard Normal Distribution and
various Student-t distributions. It shows how the distributions
converge on the SND values as the DOF increases.

```
$ ./ztables.py -s -t 10 -t 20 -t 30 -t 100 -t 200 -p 0.95 -p 0.98 -p 0.99

Probabilities to z-values Table

                     t-dist  t-dist  t-dist  t-dist  t-dist
  Probability  SND     10      20      30     100     200
  ===========  ====  ======  ======  ======  ======  ======
     95.00%    1.96    2.23    2.09    2.04    1.98    1.97
     98.00%    2.33    2.77    2.53    2.46    2.37    2.34
     99.00%    2.58    3.17    2.84    2.75    2.62    2.60
```

Comments and corrections greatly appreciated.

Enjoy!
