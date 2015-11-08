#!/usr/bin/env python
'''
This program generates z-tables for Standard Normal Distributions and
Student-t Distributions. It also generates the z-values for specific
probabilities from those distributions without using any special
libraries.

The goal is to show how to solve the probability density functions
(PDFs) and use those solutions to approximate the area under a curve
in order to obtain the z-value or, in the case of the z-value lookup,
show how to do a binary search the over probabilities to find the
desired z-value.

It is meant to be a learning tool for those that are interested in
how the tables are generated. For production code using the scipy
and numpy packages. They are optimized and much better debugged.

I explicitly chose to use the trapezoidal rule to approximate the
definite integral representing the area under the curve for two points
because it is simple, reasonably fast and accurate enough. You
could, of course, choose any other approximation method (like
Simpson's Rule).

I also chose the Lanczos approximation for solving the gamma function
for the t-distribution PDF. The algorithm for the approximation was
obtained from page 214 of "Numerical Recipes in C, 2nd Edition".

Here is an example that shows you how to generate the z-table that
contains the probabilities for the standard normal distribution using
the default parameters. This is the table that you often see in text
books.

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

Here is an example that shows how to determine the z values for 95%,
98% and 99% probabilities using the Standard Normal Distribution and
various Student-t distributions. It shows how the distributions
converge on the SND values as the DOF increases.

   $ ./ztables.py -s -t 10 -t 20 -t 30 -t 100 -t 200 -p 0.95 -p 0.98 -p 0.99
   
   Probabilities to z-values Table
   
                        t-dist  t-dist  t-dist  t-dist  t-dist
     Probability  SND     10      20      30     100     200  
     ===========  ====  ======  ======  ======  ======  ======
        95.00%    1.96    2.23    2.09    2.04    1.98    1.97
        98.00%    2.33    2.77    2.53    2.46    2.37    2.34
        99.00%    2.58    3.17    2.84    2.75    2.62    2.60
'''
import argparse
import math
import os
import sys


VERSION = '0.9.0'  # pre-release


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

    
def ztable(title, lower_bound, upper_bound, minval, intervals, fct, *args):
    '''
    Create the z values table for the standard normal distribution.
    '''
    sys.stdout.write('''
{0}

  z     .00     .01     .02     .03     .04     .05     .06     .07     .08     .09
=====  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
'''.format(title))

    write = sys.stdout.write
    flush = sys.stdout.flush

    z = lower_bound
    while z <= upper_bound:
        write('{0:>5.1f}'.format(z))

        z1 = z
        for i in range(10):
            probability = area_under_curve(minval, float(z1), intervals, fct, *args)
            write('  {0:.4f}'.format(probability))
            flush()
            z1 += 0.01 if z >= 0 else -0.01

        z += 0.1
        write('\n')

    write('\n')


def binary_search_for_z(probability, tolerance, maxtop, minval, iterations, verbose, fct, *args):
    '''
    Get the z value that matches the specified percentage.
    '''
    # Binary search to find the closest value.
    if verbose:
        print('')
    z = 0.0
    adjustment = float(maxtop) / 2.0
    top = maxtop
    bot = 0.0
    diff = tolerance * 2  # start the loop
    while diff > tolerance:
        mid = bot + ((top - bot) / 2.0)
        z = mid - adjustment
        q = area_under_curve(minval, z, iterations, fct, *args)
        current_probability = 1.0 - (2.0 * (1.0 - q))
        diff = abs(current_probability - probability)
        
        if verbose:
            print('    p={0:>.4f},  z={1:>.2f},  d={2:>.4f}, q={3:>.4f}, m={4}'.format(current_probability,
                                                                                       z,
                                                                                       diff,
                                                                                       q,
                                                                                       maxtop))
            
        if probability < current_probability:
            # It is to the right.
            top = mid
        elif probability > current_probability:
            # It is to the left.
            bot = mid
        else:
            break
        
        assert top <= maxtop
        assert bot >= 0

    return z


def probability_table(opts):
    '''
    List the z-values for specific probabilities. This is the
    inverse of the z-table operations. It uses a binary search
    to find the values.
    '''
    write = sys.stdout.write
    flush = sys.stdout.flush
    maxtop = 2 * round(abs(opts.lower_bound) + opts.upper_bound + 0.5, 0)

    write('''
Probabilities to z-values Table
''')

    # Create the header.
    h0 = '             '
    h1 = '  Probability'
    h2 = '  ==========='

    if opts.snd is True or opts.tdist is None:
        h0 += '      '
        h1 += '  SND '
        h2 += '  ===='
        
    if opts.tdist:
        for tdist in opts.tdist:
            h0 += '  t-dist'
            h1 += '  {:^6}'.format(tdist)
            h2 += '  ======'
                     
    write(h0 + '\n')
    write(h1 + '\n')
    write(h2 + '\n')
            
    for probability in opts.probability:
        write('  {0:^11.2%}'.format(probability))
        flush()

        if opts.snd is True or opts.tdist is None:
            z = binary_search_for_z(probability, 0.0001, maxtop, opts.minimum, opts.intervals, opts.verbose, pdf_snd)
            write('  {0:4.2f}'.format(z))
            flush()

        if opts.tdist:
            for tdist in opts.tdist:
                z = binary_search_for_z(probability, 0.0001, maxtop, opts.minimum, opts.intervals, opts.verbose, pdf_t, tdist)
                write('  {0:6.2f}'.format(z))
                flush()

        write('\n')

    write('\n')


def getopts():
    '''
    Get the command line options.
    '''

    def p_opt(val):
        '''
        The probability.
        Must be between .001 and .9999.
        '''
        try:
            num = float(val)
            if num < 0.001 or num > 0.9999:
                raise argparse.ArgumentTypeError('Probability is out of range: [0.001..0.9999].')
        except ValueError:
            raise argparse.ArgumentTypeError('Probability must be a floating point number.')
        return float(val)
    
    # Trick to capitalize the built-in headers.
    # Unfortunately I can't get rid of the ":" reliably.
    def gettext(s):
        lookup = {
            'usage: ': 'USAGE:',
            'optional arguments': 'OPTIONAL ARGUMENTS',
            'show this help message and exit': 'Show this help message and exit.\n ',
        }
        return lookup.get(s, s)
    
    argparse._ = gettext  # to capitalize help headers
    base = os.path.basename(sys.argv[0])
    name = os.path.splitext(base)[0]
    usage = '\n  {0} [OPTIONS]'.format(base)
    desc = 'DESCRIPTION:{0}'.format('\n  '.join(__doc__.split('\n')))
    epilog = r'''
EXAMPLES:
  $ # ================================
  $ # Example 1: help
  $ {0} --help

  $ # ================================
  $ # Example 2: Generate the z-table for the Standard Normal Distribution (SND)
  $ {0} -s

  $ # ================================
  $ # Example 3: Generate the z-table for the Student-t Distribution
  $ #            with 20 degrees of freedom
  $ {0} -t -d 20

  $ # ================================
  $ # Example 4: Generate the z-values for 90%, 95% and 98% using the SND
  $ {0} -p 0.90 -p 0.95 -p 0.98 -s

  $ # ================================
  $ # Example 5: Generate the z-values for 90%, 95% and 98% using the
  $ #            t-distribution with 20 DOF
  $ {0} -p 0.90 -p 0.95 -p 0.98 -t 20

  $ # ================================
  $ # Example 6: Generate the z-values for 95% and 98% using the SND
  $ #            and the t-distribution with 200 DOF
  $ {0} -s -t 200 -p 0.95 -p 0.98
 '''.format(base)
    afc = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=afc,
                                     description=desc[:-2],
                                     usage=usage,
                                     epilog=epilog)

    parser.add_argument('-i', '--intervals',
                        action='store',
                        type=int,
                        metavar=('INT'),
                        default=10000,
                        help='''The number of intervals in the area under the curve being
analyzed.

Each interval is a trapezoid whose width is fixed and whose
heights (y1, y2) are the PDF values of x1 and x2. The area
of the trapezoids are accumulated to obtain the total area
(Trapezoidal Rule).

The PDF (probability density function) represents the height
of the curve (y) for each x value.

Default=%(default)s. This was determined empirically.
 ''')

    parser.add_argument('-l', '--lower-bound',
                        action='store',
                        type=float,
                        metavar=('FLOAT'),
                        default=-3.4,
                        help='''The left most z-value to report. This is the number of
standard deviations to the left of the center of the curve.
Default=%(default)s.
 ''')
    
    parser.add_argument('-m', '--minimum',
                        action='store',
                        type=float,
                        metavar=('FLOAT'),
                        default=-7.0,
                        help='''The minimum value used as the left most bound for
calculations.  This is the number of standard deviations to
the left of the center of the curve that is used to
represent minus infinity. It needs be large enough that a
reasonable value can be obtained for the lower bound. It
must be smaller than the lower-bound.
Default=%(default)s. This was deemed sufficient empirically.
 ''')

    parser.add_argument('-p', '--probability',
                        action='append',
                        type=p_opt,
                        metavar=('FLOAT'),
                        help='''The probablity to generate a z-value for. It must be number
between 0.001 and 0.999. Generate the z-value for this
probability using the SND.
 ''')

    parser.add_argument('-s', '--snd',
                        action='store_true',
                        help='''Generate the z-table of probabilities for the Standard
Normal Distribution unless -p is specified.  In that case
generate the z-values for the probabilities using the SND.
 ''')

    parser.add_argument('-t', '--tdist',
                        action='append',
                        type=int,
                        metavar=('DOF'),
                        help='''Generate the z-table of probabilities for a Student-t
distribution with DOF degrees of freedom.
 ''')

    parser.add_argument('-u', '--upper-bound',
                        action='store',
                        type=float,
                        metavar=('FLOAT'),
                        default=3.49,
                        help='''The right most z-value to report. This is the number of
standard deviations to the left of the center of the curve.
Default=%(default)s.
 ''')

    parser.add_argument('-v', '--verbose',
                        action='count',
                        help='''Increase the level of verbosity.
 ''')

    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s v{0}'.format(VERSION),
                        help="""Show program's version number and exit.
 """)
    
    opts = parser.parse_args()
    return opts


def main():
    '''
    Main
    '''
    opts = getopts()
    if opts.probability is None:
        # Generate the z-tables.
        if opts.snd is True:
            title = 'z-Table for Standard Normal Distribution ({0:,})'.format(opts.intervals)
            ztable(title,
                   opts.lower_bound,
                   opts.upper_bound,
                   opts.minimum,
                   opts.intervals,
                   pdf_snd)
            
        if opts.tdist is not None:
            for tdist in opts.tdist:
                title = 'z-Table for Student-t Distribution '
                title += '({0:,}, {1} DOF)'.format(opts.intervals, tdist)
                ztable(title,
                       opts.lower_bound,
                       opts.upper_bound,
                       opts.minimum,
                       opts.intervals,
                       pdf_t,
                       tdist)
                
        if opts.tdist is None and opts.snd is False:
            print('WARNING: You must specify -s or -t to generate a z-value table, see -h for more details.')
    else:
        probability_table(opts)


if __name__ == '__main__':
    main()
