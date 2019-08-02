package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;

public class BrentRoot implements Reporter {
  /*
   * BrentRoot employs the Brent alorithm to find the root (zero) of
   * a function of one variable, x. The return value is the value of x 
   * at the root.
   * 
   * The arguments are:
   * fnctn - the function, passed as a NetLogo task variable.
   * a and b - the bounds between which the root is to be 
   * found. Note that if the range between the two bounds does not contain
   * a root of the function, that is if the function evaluated at x = a
   * does not have a different sign than the function evaluated at
   * x = b, then an exception is thrown. (It might be useful at some point
   * to add a bounds-finding routine that would allow a single initial 
   * guess or that would fix invalid bounds.)
   * rtol - the tolerance to which the solution is taken, as a proportion of
   * x at the root. If the root occurs very close to x = 0, the 
   * tolerance is set to a small, positive number.
   * atol - the absolute tolerance to which the solution is taken.
   */

  private double relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
  private double absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
  private static final int maxEvaluations = 1000;

  @Override
  public Syntax getSyntax() {
    return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
        Syntax.NumberType(), Syntax.NumberType(),
        Syntax.NumberType() | Syntax.RepeatableType()}, Syntax.NumberType(), 3);
  }

  @Override
  public Object report(Argument args[], Context context)
          throws ExtensionException, LogoException {

    final double SQRT_DBL_EPSILON = Math.sqrt(Double.MIN_NORMAL);

    AnonymousReporter fnctn = args[0].getReporter();
    double a = args[1].getDoubleValue();
    double b = args[2].getDoubleValue();
    if (args.length > 3) {
      double xx = args[3].getDoubleValue();
      relTolerance = (xx >= 0.0) ? xx : Bounds.RELATIVE_TOLERANCE_DEFAULT;
    }
    if (args.length > 4) {
      double xx = args[4].getDoubleValue();
      absTolerance = (xx >= 0.0) ? xx : Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
    }

    // Need to make sure both tols are not zero. (One is okay.)
    // If they both are, assume user wanted defaults for both.
    if (relTolerance == 0.0 && absTolerance == 0.0) {
      relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
      absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
    }

    // Evaluate the function at its bounds and make sure that it 
    // brackets zero.
    double fa = NumAnalUtils.getFofX(a, fnctn, context);
    double fb = NumAnalUtils.getFofX(b, fnctn, context);
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
      throw new ExtensionException("Brent-root: " + a
              + " and " + b + " do not bracket the root.");
    }

    double c = b;
    double fc = fb;
    double d = b - a;
    double e = d;

    for (int iter = 0; iter < maxEvaluations; iter++) {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
        // Need to rename a, b and c.
        c = a;
        fc = fa;
        d = b - a;
        e = d;
      }
      if (Math.abs(fc) < Math.abs(fb)) {
        a = b;
        fa = fb;
        b = c;
        fb = fc;
        c = a;
        fc = fa;
      }

      // Check for convergence to a root.
      double tol1 = Math.max(absTolerance, Math.abs(b) * relTolerance);
      tol1 = Math.max(tol1, SQRT_DBL_EPSILON);
      double xmid = 0.5 * (c - b);
      if (Math.abs(xmid) <= tol1 || fb == 0.0) {
        return (b);
      }

      // Otherwise try inverse quadratic interpolation.
      double p, q;
      if (Math.abs(e) >= tol1 && Math.abs(fa) >= Math.abs(fa)) {
        double s = fb / fa;
        // test for a == c uses abs diff < upl().
        if (Math.abs(a - c) < Math.ulp(a)) {
          p = 2.0 * xmid * s;
          q = 1.0 - s;
        } else {
          q = fa / fc;
          double r = fb / fc;
          p = s * (2.0 * xmid * q * (q - r) - (b - a) * (r - 1.0));
          q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        }
        if (p > 0.0) {
          q = -q;
        }
        p = Math.abs(p);
        double min1 = 3.0 * xmid * q - Math.abs(tol1 * q);
        double min2 = Math.abs(e * q);
        if (2.0 * p < Math.min(min1, min2)) {
          // Interpolate
          e = d;
          d = p / q;
        } else {
          // No luck, try bisection instead.
          d = xmid;
          e = d;
        }
      } else {
        // Convergence is too slow, just do a bisection.
        d = xmid;
        e = d;
      }

      // Move best guess to a.
      a = b;
      fa = fb;
      b += (Math.abs(d) > tol1) ? d : tol1 * ((xmid < 0.0) ? -1 : 1);
      fb = NumAnalUtils.getFofX(b, fnctn, context);
    }

    // We're out of the loop. Too many iterations.
    throw new ExtensionException("Brent-root: Exceeded the maximum number"
            + " of steps: " + maxEvaluations);
  }
}
