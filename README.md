# ztables
Ever wonder how Standard Normal and Student-t Distribution z-tables are generated?

## Introduction
This program shows how to generate z-tables for Standard Normal
Distributions and Student-t Distributions. It also shows how to
generate the z-values for specific probabilities from those
distributions without using any special libraries.

I have tested the program using Python 2.7.10 and 3.5.

## Discussion

The goal is to show how to solve the probability density functions
(PDFs) and use those solutions to approximate the area under a curve
in order to obtain the z-value or, in the case of the z-value lookup,
show how to do a binary search over the probabilities to find the
desired z-value.

It is meant to be a learning tool for those that are interested in
how the tables are generated. For production code use the scipy
and numpy packages. They are optimized and much better debugged.

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

In the sections below some of the code fragments are provided as well
as the program options and some sample outputs.

## Program Options
These are the program options.

Option        | Long                      | Description
------------- | ------------------------- | -----------
-h            | --help                    | Help.
-i INT        | --intervals INT           | The number intervals used for the area calculation.
-l&nbsp;FLOAT | --lower-bound&nbsp;FLOAT  | The leftmost z-value to report.
-m&nbsp;FLOAT | --minimum FLOAT           | z-value at minus infinity.
-p&nbsp;FLOAT | --probability FLOAT       | Find the z-value for this probability. Must be in the range 0.001 to .9999.
-s            | --snd                     | Generate the SND table unless -p specified. If -p is specified use the SND to find the z-value.
-t DOF        | --tdist DOF               | Generate the t-dist table with DOF degress of freedom unless -p specified. If -p is specified use the t-dist to find the z-value.
-u&nbsp;FLOAT | --upper-bound&nbsp;FLOAT  | The rightmost z-value to report.
-v            | --verbose                 | Increase the level of verbosity. Only used for debugging -p.
-V            | --version                 | Print the program version and exit.

## SND PDF
This is the code for the probability density function (PDF) for the
standard normal distribution (SND).

```python
def pdf_snd(x):
    '''
    Calculate the probability density function (PDF) for a standard
    normal distribution.

    s = standard deviation (1 for a standard normal distribution)
    u = mean (0 for a standard normal distribution)

    This is the height of the curve at x.

    It is exactly the same as pdf_nd(x, 1, 0) but is somewhat more
    efficient.
    '''
    dx2 = float(x) ** 2
    den = math.sqrt(2 * math.pi)
    exp = math.e ** - (dx2 / 2)
    y =  exp / den
    return y
```

## PDF for the normal distribution
This is the code for the probability density function (PDF) for the
normal distribution where the mean and the standard deviation are
explicitly specified.

```python
def pdf_nd(x, s=1.0, u=0.0):
    '''
    Calculate the probability density function (PDF) for a normal
    distribution.

    s = standard deviation (1 for a standard normal distribution)
    u = mean (0 for a standard normal distribution)

    This is the height of the curve at x.
    '''
    dx = float(x) - float(u)
    dx2 = dx * dx
    xden = 2 * (s ** 2)
    den = s * math.sqrt(2 * math.pi)
    exp = math.e ** ( -dx2 / xden )
    y =  exp / den
    return y
```

## T-dist PDF
This is the code for the t-distribution PDF.

```python
def pdf_t(x, dof):
    '''
    Calculate the probability density function (PDF) at x for a
    student-t distribution with dof degrees of freedom.

    This is basically the height of the curve at x.
    '''
    assert dof > 2

    x1 = gamma((float(dof) + 1.0) / 2.0)
    x2 = math.sqrt(dof * math.pi) * gamma((float(dof) / 2.0))
    x3 = 1.0 + (float((x ** 2)) / float(dof))
    x4 = float((dof + 1)) / 2.0
    x5 = x3 ** -x4

    y = (x1 * x5) / x2
    return y
```

## Gamma Function
This is the code that estimates the gamma function for t-distribution calculations.

```python
def gamma(x):
    '''
    Gamma function.

    Uses the Lanczos approximation and natural logarithms.

    For integer values of x we can use the exact value of (x-1)!.

       gamma(1/2) = 1.77245385091
       gamma(3/2) = 0.886226925453
       gamma(5/2) = 1.32934038818
       gamma(7/2) = 3.32335097045
       gamma(4)   = 6.0
    '''
    if (x - int(x)) == 0:
        # Optimization for integer values: (x-1)!.
        return reduce(lambda a, b: a * b, [float(i) for i in range(1, int(x))])

    # Lanczos approximation, page 214 of Numerical Recipes in C.
    c = [76.18009172947146,
         -86.50532032941677,
         24.01409824083091,
         -1.231739572450155,
         0.1208650973866179e-2,
         -0.5395239384953e-5,
    ]
    c0 = 1.000000000190015
    c1 = 2.5066282746310005
    x1 = float(x) + 5.5
    x2 = (float(x) + 0.5) * math.log(x1)
    x3 = x1 - x2
    x4 = c0
    x5 = float(x)
    for i in range(6):
        x5 += 1.0
        x4 += c[i] / x5
    x6 = math.log((c1 * x4) / float(x))
    x7 = -x3 + x6  # ln(gamma(x))
    g = math.exp(x7)
    return g
```

## Estimate the Area Under the Curve using the Trapezoidal Rule
This function shows how to estimate the area under the curve using the
Trapezoidal Rule. There many other methods available calculating area
under the curve for definite integrals but this one is reasonably fast
and sufficiently accurate for this purpose.


```python
def area_under_curve(x1, x2, intervals, fct, *args, **kwargs):
    '''
    Calculate the approximate area under a curve using trapezoidal
    approximation.

    It breaks the interval between x1 and x2 into trapezoids whose
    width is fixed (proportional to how the interval is sliced). The
    height of each rectangle is the pdf function value for x at the
    start of the interval. The accumulation of the areas provides an
    estimate of the area under the curve.

    The greater the number of intervals the better the estimate is at
    the cost of performance.
    '''
    assert x2 > x1  # just a sanity check
    assert intervals > 1  # another sanity check
    
    total_area = 0.0
    width = (float(x2) - float(x1)) / float(intervals)
    
    x = float(x1)
    py = float(fct(x, *args, **kwargs))
    for i in range(intervals):
        y = float(fct(x, *args, **kwargs))
        rectangle_area = width * y  # area of rectangle at x with height y
        triangle_area = ((y - py) * width) / 2.0  # adjustment based on height change
        total_area += rectangle_area + triangle_area  # trapezoid area
        x += width  # advance to the next edge
        py = y  # remember the previous height
        
    return total_area
```

## Here is the SND Table.
Here is the Standard Normal Distribution Table generated by this program using the command: '`./ztables.py -s`'.

```
z-Table for Standard Normal Distribution (10,000)

  z     .00     .01     .02     .03     .04     .05     .06     .07     .08     .09
=====  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
 -3.4  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0003  0.0002
 -3.3  0.0005  0.0005  0.0005  0.0004  0.0004  0.0004  0.0004  0.0004  0.0004  0.0003
 -3.2  0.0007  0.0007  0.0006  0.0006  0.0006  0.0006  0.0006  0.0005  0.0005  0.0005
 -3.1  0.0010  0.0009  0.0009  0.0009  0.0008  0.0008  0.0008  0.0008  0.0007  0.0007
 -3.0  0.0013  0.0013  0.0013  0.0012  0.0012  0.0011  0.0011  0.0011  0.0010  0.0010
 -2.9  0.0019  0.0018  0.0018  0.0017  0.0016  0.0016  0.0015  0.0015  0.0014  0.0014
 -2.8  0.0026  0.0025  0.0024  0.0023  0.0023  0.0022  0.0021  0.0021  0.0020  0.0019
 -2.7  0.0035  0.0034  0.0033  0.0032  0.0031  0.0030  0.0029  0.0028  0.0027  0.0026
 -2.6  0.0047  0.0045  0.0044  0.0043  0.0041  0.0040  0.0039  0.0038  0.0037  0.0036
 -2.5  0.0062  0.0060  0.0059  0.0057  0.0055  0.0054  0.0052  0.0051  0.0049  0.0048
 -2.4  0.0082  0.0080  0.0078  0.0075  0.0073  0.0071  0.0069  0.0068  0.0066  0.0064
 -2.3  0.0107  0.0104  0.0102  0.0099  0.0096  0.0094  0.0091  0.0089  0.0087  0.0084
 -2.2  0.0139  0.0136  0.0132  0.0129  0.0125  0.0122  0.0119  0.0116  0.0113  0.0110
 -2.1  0.0179  0.0174  0.0170  0.0166  0.0162  0.0158  0.0154  0.0150  0.0146  0.0143
 -2.0  0.0228  0.0222  0.0217  0.0212  0.0207  0.0202  0.0197  0.0192  0.0188  0.0183
 -1.9  0.0287  0.0281  0.0274  0.0268  0.0262  0.0256  0.0250  0.0244  0.0239  0.0233
 -1.8  0.0359  0.0351  0.0344  0.0336  0.0329  0.0322  0.0314  0.0307  0.0301  0.0294
 -1.7  0.0446  0.0436  0.0427  0.0418  0.0409  0.0401  0.0392  0.0384  0.0375  0.0367
 -1.6  0.0548  0.0537  0.0526  0.0516  0.0505  0.0495  0.0485  0.0475  0.0465  0.0455
 -1.5  0.0668  0.0655  0.0643  0.0630  0.0618  0.0606  0.0594  0.0582  0.0571  0.0559
 -1.4  0.0808  0.0793  0.0778  0.0764  0.0749  0.0735  0.0721  0.0708  0.0694  0.0681
 -1.3  0.0968  0.0951  0.0934  0.0918  0.0901  0.0885  0.0869  0.0853  0.0838  0.0823
 -1.2  0.1151  0.1131  0.1112  0.1093  0.1075  0.1056  0.1038  0.1020  0.1003  0.0985
 -1.1  0.1357  0.1335  0.1314  0.1292  0.1271  0.1251  0.1230  0.1210  0.1190  0.1170
 -1.0  0.1587  0.1562  0.1539  0.1515  0.1492  0.1469  0.1446  0.1423  0.1401  0.1379
 -0.9  0.1841  0.1814  0.1788  0.1762  0.1736  0.1711  0.1685  0.1660  0.1635  0.1611
 -0.8  0.2119  0.2090  0.2061  0.2033  0.2005  0.1977  0.1949  0.1922  0.1894  0.1867
 -0.7  0.2420  0.2389  0.2358  0.2327  0.2296  0.2266  0.2236  0.2206  0.2177  0.2148
 -0.6  0.2743  0.2709  0.2676  0.2643  0.2611  0.2578  0.2546  0.2514  0.2483  0.2451
 -0.5  0.3085  0.3050  0.3015  0.2981  0.2946  0.2912  0.2877  0.2843  0.2810  0.2776
 -0.4  0.3446  0.3409  0.3372  0.3336  0.3300  0.3264  0.3228  0.3192  0.3156  0.3121
 -0.3  0.3821  0.3783  0.3745  0.3707  0.3669  0.3632  0.3594  0.3557  0.3520  0.3483
 -0.2  0.4207  0.4168  0.4129  0.4090  0.4052  0.4013  0.3974  0.3936  0.3897  0.3859
 -0.1  0.4602  0.4562  0.4522  0.4483  0.4443  0.4404  0.4364  0.4325  0.4286  0.4247
  0.0  0.5000  0.5040  0.5080  0.5120  0.5160  0.5199  0.5239  0.5279  0.5319  0.5359
  0.1  0.5398  0.5438  0.5478  0.5517  0.5557  0.5596  0.5636  0.5675  0.5714  0.5753
  0.2  0.5793  0.5832  0.5871  0.5910  0.5948  0.5987  0.6026  0.6064  0.6103  0.6141
  0.3  0.6179  0.6217  0.6255  0.6293  0.6331  0.6368  0.6406  0.6443  0.6480  0.6517
  0.4  0.6554  0.6591  0.6628  0.6664  0.6700  0.6736  0.6772  0.6808  0.6844  0.6879
  0.5  0.6915  0.6950  0.6985  0.7019  0.7054  0.7088  0.7123  0.7157  0.7190  0.7224
  0.6  0.7257  0.7291  0.7324  0.7357  0.7389  0.7422  0.7454  0.7486  0.7517  0.7549
  0.7  0.7580  0.7611  0.7642  0.7673  0.7704  0.7734  0.7764  0.7794  0.7823  0.7852
  0.8  0.7881  0.7910  0.7939  0.7967  0.7995  0.8023  0.8051  0.8078  0.8106  0.8133
  0.9  0.8159  0.8186  0.8212  0.8238  0.8264  0.8289  0.8315  0.8340  0.8365  0.8389
  1.0  0.8413  0.8438  0.8461  0.8485  0.8508  0.8531  0.8554  0.8577  0.8599  0.8621
  1.1  0.8643  0.8665  0.8686  0.8708  0.8729  0.8749  0.8770  0.8790  0.8810  0.8830
  1.2  0.8849  0.8869  0.8888  0.8907  0.8925  0.8944  0.8962  0.8980  0.8997  0.9015
  1.3  0.9032  0.9049  0.9066  0.9082  0.9099  0.9115  0.9131  0.9147  0.9162  0.9177
  1.4  0.9192  0.9207  0.9222  0.9236  0.9251  0.9265  0.9279  0.9292  0.9306  0.9319
  1.5  0.9332  0.9345  0.9357  0.9370  0.9382  0.9394  0.9406  0.9418  0.9429  0.9441
  1.6  0.9452  0.9463  0.9474  0.9484  0.9495  0.9505  0.9515  0.9525  0.9535  0.9545
  1.7  0.9554  0.9564  0.9573  0.9582  0.9591  0.9599  0.9608  0.9616  0.9625  0.9633
  1.8  0.9641  0.9649  0.9656  0.9664  0.9671  0.9678  0.9686  0.9693  0.9699  0.9706
  1.9  0.9713  0.9719  0.9726  0.9732  0.9738  0.9744  0.9750  0.9756  0.9761  0.9767
  2.0  0.9772  0.9778  0.9783  0.9788  0.9793  0.9798  0.9803  0.9808  0.9812  0.9817
  2.1  0.9821  0.9826  0.9830  0.9834  0.9838  0.9842  0.9846  0.9850  0.9854  0.9857
  2.2  0.9861  0.9864  0.9868  0.9871  0.9875  0.9878  0.9881  0.9884  0.9887  0.9890
  2.3  0.9893  0.9896  0.9898  0.9901  0.9904  0.9906  0.9909  0.9911  0.9913  0.9916
  2.4  0.9918  0.9920  0.9922  0.9925  0.9927  0.9929  0.9931  0.9932  0.9934  0.9936
  2.5  0.9938  0.9940  0.9941  0.9943  0.9945  0.9946  0.9948  0.9949  0.9951  0.9952
  2.6  0.9953  0.9955  0.9956  0.9957  0.9959  0.9960  0.9961  0.9962  0.9963  0.9964
  2.7  0.9965  0.9966  0.9967  0.9968  0.9969  0.9970  0.9971  0.9972  0.9973  0.9974
  2.8  0.9974  0.9975  0.9976  0.9977  0.9977  0.9978  0.9979  0.9979  0.9980  0.9981
  2.9  0.9981  0.9982  0.9982  0.9983  0.9984  0.9984  0.9985  0.9985  0.9986  0.9986
  3.0  0.9987  0.9987  0.9987  0.9988  0.9988  0.9989  0.9989  0.9989  0.9990  0.9990
  3.1  0.9990  0.9991  0.9991  0.9991  0.9992  0.9992  0.9992  0.9992  0.9993  0.9993
  3.2  0.9993  0.9993  0.9994  0.9994  0.9994  0.9994  0.9994  0.9995  0.9995  0.9995
  3.3  0.9995  0.9995  0.9995  0.9996  0.9996  0.9996  0.9996  0.9996  0.9996  0.9997
  3.4  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9997  0.9998
```

## Here is a Student-t Table with 20 DOF.
Here is the Student-t Table with 20 degress of freedom (DOF) generated by this program using the command '`ztables.py -t 20`'.

```
z-Table for Student-t Distribution (10,000, 20 DOF)

  z     .00     .01     .02     .03     .04     .05     .06     .07     .08     .09
=====  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
 -3.4  0.0014  0.0014  0.0014  0.0013  0.0013  0.0013  0.0012  0.0012  0.0012  0.0012
 -3.3  0.0018  0.0017  0.0017  0.0017  0.0016  0.0016  0.0016  0.0015  0.0015  0.0015
 -3.2  0.0022  0.0022  0.0021  0.0021  0.0021  0.0020  0.0020  0.0019  0.0019  0.0018
 -3.1  0.0028  0.0028  0.0027  0.0026  0.0026  0.0025  0.0025  0.0024  0.0024  0.0023
 -3.0  0.0035  0.0035  0.0034  0.0033  0.0032  0.0032  0.0031  0.0030  0.0030  0.0029
 -2.9  0.0044  0.0043  0.0042  0.0041  0.0040  0.0040  0.0039  0.0038  0.0037  0.0036
 -2.8  0.0055  0.0054  0.0053  0.0052  0.0051  0.0049  0.0048  0.0047  0.0046  0.0045
 -2.7  0.0069  0.0067  0.0066  0.0065  0.0063  0.0062  0.0060  0.0059  0.0058  0.0057
 -2.6  0.0086  0.0084  0.0082  0.0080  0.0079  0.0077  0.0075  0.0074  0.0072  0.0070
 -2.5  0.0106  0.0104  0.0102  0.0100  0.0097  0.0095  0.0093  0.0091  0.0089  0.0088
 -2.4  0.0131  0.0129  0.0126  0.0123  0.0121  0.0118  0.0116  0.0113  0.0111  0.0108
 -2.3  0.0162  0.0158  0.0155  0.0152  0.0149  0.0146  0.0143  0.0140  0.0137  0.0134
 -2.2  0.0199  0.0195  0.0191  0.0187  0.0183  0.0179  0.0176  0.0172  0.0169  0.0165
 -2.1  0.0243  0.0238  0.0234  0.0229  0.0224  0.0220  0.0215  0.0211  0.0207  0.0203
 -2.0  0.0296  0.0291  0.0285  0.0279  0.0274  0.0269  0.0263  0.0258  0.0253  0.0248
 -1.9  0.0360  0.0353  0.0346  0.0340  0.0333  0.0327  0.0320  0.0314  0.0308  0.0302
 -1.8  0.0435  0.0427  0.0419  0.0411  0.0403  0.0396  0.0388  0.0381  0.0374  0.0367
 -1.7  0.0523  0.0514  0.0504  0.0495  0.0486  0.0477  0.0468  0.0460  0.0451  0.0443
 -1.6  0.0626  0.0615  0.0604  0.0594  0.0583  0.0573  0.0563  0.0552  0.0543  0.0533
 -1.5  0.0746  0.0733  0.0721  0.0708  0.0696  0.0684  0.0672  0.0661  0.0649  0.0638
 -1.4  0.0884  0.0870  0.0855  0.0841  0.0827  0.0813  0.0799  0.0786  0.0772  0.0759
 -1.3  0.1042  0.1025  0.1009  0.0992  0.0976  0.0960  0.0945  0.0929  0.0914  0.0899
 -1.2  0.1221  0.1202  0.1183  0.1165  0.1147  0.1129  0.1111  0.1093  0.1076  0.1059
 -1.1  0.1422  0.1401  0.1380  0.1359  0.1339  0.1319  0.1299  0.1279  0.1259  0.1240
 -1.0  0.1646  0.1623  0.1600  0.1577  0.1554  0.1531  0.1509  0.1487  0.1465  0.1443
 -0.9  0.1894  0.1868  0.1843  0.1817  0.1792  0.1767  0.1743  0.1718  0.1694  0.1670
 -0.8  0.2166  0.2137  0.2109  0.2082  0.2054  0.2027  0.2000  0.1973  0.1947  0.1920
 -0.7  0.2460  0.2429  0.2399  0.2369  0.2339  0.2310  0.2281  0.2251  0.2223  0.2194
 -0.6  0.2776  0.2744  0.2711  0.2679  0.2647  0.2615  0.2584  0.2553  0.2521  0.2491
 -0.5  0.3113  0.3078  0.3044  0.3010  0.2976  0.2942  0.2908  0.2875  0.2842  0.2809
 -0.4  0.3467  0.3431  0.3395  0.3359  0.3323  0.3288  0.3252  0.3217  0.3182  0.3147
 -0.3  0.3836  0.3799  0.3761  0.3724  0.3687  0.3650  0.3613  0.3576  0.3540  0.3503
 -0.2  0.4217  0.4179  0.4141  0.4102  0.4064  0.4026  0.3988  0.3950  0.3912  0.3874
 -0.1  0.4607  0.4568  0.4528  0.4489  0.4450  0.4411  0.4372  0.4334  0.4295  0.4256
  0.0  0.5000  0.5039  0.5079  0.5118  0.5158  0.5197  0.5236  0.5276  0.5315  0.5354
  0.1  0.5393  0.5432  0.5472  0.5511  0.5550  0.5589  0.5628  0.5666  0.5705  0.5744
  0.2  0.5782  0.5821  0.5859  0.5898  0.5936  0.5974  0.6012  0.6050  0.6088  0.6126
  0.3  0.6164  0.6201  0.6239  0.6276  0.6313  0.6350  0.6387  0.6424  0.6460  0.6497
  0.4  0.6533  0.6569  0.6605  0.6641  0.6677  0.6712  0.6748  0.6783  0.6818  0.6853
  0.5  0.6887  0.6922  0.6956  0.6990  0.7024  0.7058  0.7092  0.7125  0.7158  0.7191
  0.6  0.7224  0.7256  0.7289  0.7321  0.7353  0.7385  0.7416  0.7447  0.7478  0.7509
  0.7  0.7540  0.7570  0.7601  0.7631  0.7661  0.7690  0.7719  0.7748  0.7777  0.7806
  0.8  0.7834  0.7863  0.7891  0.7918  0.7946  0.7973  0.8000  0.8027  0.8053  0.8080
  0.9  0.8106  0.8132  0.8157  0.8183  0.8208  0.8233  0.8257  0.8282  0.8306  0.8330
  1.0  0.8354  0.8377  0.8400  0.8423  0.8446  0.8469  0.8491  0.8513  0.8535  0.8557
  1.1  0.8578  0.8599  0.8620  0.8641  0.8661  0.8681  0.8701  0.8721  0.8741  0.8760
  1.2  0.8779  0.8798  0.8817  0.8835  0.8853  0.8871  0.8889  0.8907  0.8924  0.8941
  1.3  0.8958  0.8975  0.8991  0.9008  0.9024  0.9040  0.9055  0.9071  0.9086  0.9101
  1.4  0.9116  0.9130  0.9145  0.9159  0.9173  0.9187  0.9201  0.9214  0.9228  0.9241
  1.5  0.9254  0.9267  0.9279  0.9292  0.9304  0.9316  0.9328  0.9339  0.9351  0.9362
  1.6  0.9374  0.9385  0.9396  0.9406  0.9417  0.9427  0.9437  0.9448  0.9457  0.9467
  1.7  0.9477  0.9486  0.9496  0.9505  0.9514  0.9523  0.9532  0.9540  0.9549  0.9557
  1.8  0.9565  0.9573  0.9581  0.9589  0.9597  0.9604  0.9612  0.9619  0.9626  0.9633
  1.9  0.9640  0.9647  0.9654  0.9660  0.9667  0.9673  0.9680  0.9686  0.9692  0.9698
  2.0  0.9704  0.9709  0.9715  0.9721  0.9726  0.9731  0.9737  0.9742  0.9747  0.9752
  2.1  0.9757  0.9762  0.9766  0.9771  0.9776  0.9780  0.9785  0.9789  0.9793  0.9797
  2.2  0.9801  0.9805  0.9809  0.9813  0.9817  0.9821  0.9824  0.9828  0.9831  0.9835
  2.3  0.9838  0.9842  0.9845  0.9848  0.9851  0.9854  0.9857  0.9860  0.9863  0.9866
  2.4  0.9869  0.9871  0.9874  0.9877  0.9879  0.9882  0.9884  0.9887  0.9889  0.9892
  2.5  0.9894  0.9896  0.9898  0.9900  0.9903  0.9905  0.9907  0.9909  0.9911  0.9912
  2.6  0.9914  0.9916  0.9918  0.9920  0.9921  0.9923  0.9925  0.9926  0.9928  0.9930
  2.7  0.9931  0.9933  0.9934  0.9935  0.9937  0.9938  0.9940  0.9941  0.9942  0.9943
  2.8  0.9945  0.9946  0.9947  0.9948  0.9949  0.9951  0.9952  0.9953  0.9954  0.9955
  2.9  0.9956  0.9957  0.9958  0.9959  0.9960  0.9960  0.9961  0.9962  0.9963  0.9964
  3.0  0.9965  0.9965  0.9966  0.9967  0.9968  0.9968  0.9969  0.9970  0.9970  0.9971
  3.1  0.9972  0.9972  0.9973  0.9974  0.9974  0.9975  0.9975  0.9976  0.9976  0.9977
  3.2  0.9978  0.9978  0.9979  0.9979  0.9979  0.9980  0.9980  0.9981  0.9981  0.9982
  3.3  0.9982  0.9983  0.9983  0.9983  0.9984  0.9984  0.9984  0.9985  0.9985  0.9985
  3.4  0.9986  0.9986  0.9986  0.9987  0.9987  0.9987  0.9988  0.9988  0.9988  0.9988
```

## Final
Comments and corrections greatly appreciated.

Enjoy!

