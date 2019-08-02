package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;

public class BrentMinimize implements Reporter {

  private static final int maxEvaluations = 1000;
  private double relTolerance;
  private double absTolerance;
 

  public BrentMinimize() {
    this.absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
    this.relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
  }

  @Override
  public Syntax getSyntax() {
    return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
        Syntax.NumberType(), Syntax.NumberType(),
        Syntax.NumberType() | Syntax.RepeatableType()}, Syntax.NumberType(), 3);
  }

  @Override
  public Object report(Argument args[], Context context)
          throws ExtensionException, LogoException {
    // Takes at least three arguments: the task function itself and
    // the lower and upper bounds of the solution interval.
    // Optionally takes the relative tolerance.
    
    // Calculate the Golden Ratio (approximately 0.3819660) and the 
    // square root of the precision for doubles. Set the maximum 
    // number of steps allowed.
    final double GOLDR = 0.5 * (3.0 - Math.sqrt(5.0));
    final double SQRT_DBL_EPSILON = Math.sqrt(Double.MIN_NORMAL);

    AnonymousReporter fnctn = args[0].getReporter();
    double lowBound = args[1].getDoubleValue();
    double highBound = args[2].getDoubleValue();
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
    
    // Check for "real" bounds in case Bounds.UPPER_ or LOWER_BOUND_DEFAULT
    // were passed.
    if (highBound == Double.POSITIVE_INFINITY) {
      highBound = Double.MAX_VALUE;
    }
    if (lowBound == Double.NEGATIVE_INFINITY) {
      lowBound = -Double.MAX_VALUE;
    }

    // Make sure that the lower and upper bounds are in ascending order.
    if (lowBound > highBound) {
      double temp = lowBound;
      lowBound = highBound;
      highBound = temp;
    }
    /*
     * x is the initial point in any step, the best point found so far, 
     * or the most recent one in case of a tie. 
     * w is the second best point.
     * v is the previous value of w.
     * u is the point that was most recently evaluated. In general,
     * u = x + d, where d is the distance moved in the current step,
     * either through a golden section or a parabolic fit.
     */
    double x = lowBound + GOLDR * (highBound - lowBound);
    double v = x;
    double w = x;
    double fx = NumAnalUtils.getFofX(x, fnctn, context);
    double fw = fx;
    double fv = fx;
    double e = 0.0;
    double d = 0.0;

    for (int i = 0; i < maxEvaluations; i++) {
      double midPt = 0.5 * (lowBound + highBound);
      // the stopping tolerance is the larger of the absolute tolerance and 
      // the relative tolerance times the current value of x, but not smaller
      // than the "very small" SQRT_DBL_EPSILON (in case absTolerance = 0 and 
      // x is very small).
      double tol1 = Math.max(absTolerance, Math.abs(x) * relTolerance);
      tol1 = Math.max(tol1, SQRT_DBL_EPSILON);
      double tol2 = 2.0 * tol1;
      // Check stopping criterion. If it is satisfied, return the 
      // current best point.
      if (Math.abs(x - midPt) <= (tol2 - 0.5 * (highBound - lowBound))) {
        return x;
      }

      if (Math.abs(e) > tol1) {
        // Fit a parabola through x, w and v, and check to see if 
        // it fits well.
        double r = (x - w) * (fx - fv);
        double q = (x - v) * (fx - fw);
        double p = ((x - v) * q) - ((x - w) * r);
        q = 2.0 * (q - r);
        p = (q > 0.0) ? -p : p;
        q = Math.abs(q);
        double etemp = e;
        e = d;
        if ((Math.abs(p) >= Math.abs(0.5 * q * etemp))
                || (p <= q * (lowBound - x))
                || (p >= q * (highBound - x))) {
          // The parabola does not fit well. Use a golden section
          // step instead.
          e = (x >= midPt) ? (lowBound - x) : (highBound - x);
          d = GOLDR * e;
        } else {
          // The parabola fits well enough, use a parabolic 
          // interpolation step.
          d = p / q;
          double uu = x + d;
          // fnctn must not be evaluated too close to the 
          // boundaries. Shorten the step if necessary.
          if ((uu - lowBound) < tol2 || (highBound - uu) < tol2) {
            d = (midPt > x) ? tol1 : -tol1;
          }

        }
      } else {
        // Go straight to a golden section step.
        e = (x >= midPt) ? (lowBound - x) : (highBound - x);
        d = GOLDR * e;
      }
      // Now take the step of distance, d, from the current point, x, 
      // to the new point, u. Take care that d is at least
      // as great as tol1. If not, step the distance tol1 instead as
      // fnctn must not be evaluated too close to x.
      double u = x + ((Math.abs(d) > +tol1) ? d : ((d > 0) ? tol1 : -tol1));

      // Evaluate the function at the new point, u.
      double fu = NumAnalUtils.getFofX(u, fnctn, context);

      // Update lowBound, highBound, v, w, and x
      if (fu <= fx) {
        // The new point is at least as good as the old. Make 
        // the old best point the new high or low bound and set x
        // to the new point.
        if (u >= x) {
          lowBound = x;
        } else {
          highBound = x;
        }
        v = w;
        fv = fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
      } else {
        // The new point is worse than the old. Replace the 
        // high or low bound with the new point and enter the next
        // iteration with the same x.
        if (u < x) {
          lowBound = u;
        } else {
          highBound = u;
        }
        // changed <= comparison to !> , and == to an approximately equal
        // test to avoid testing doubles for equality.
        if (!(fu > fw) || Math.abs(w - x) < Math.ulp(w)) {
          v = w;
          fv = fw;
          w = u;
          fw = fu;
          // change <= to !>, and == to abs of diff < ulp, so as to avoid 
          // the equality comparison to doubles.
        } else if (!(fu > fv) || Math.abs(v - x) < Math.ulp(v) 
                || Math.abs(v - w) <Math.ulp(v)) {
          v = u;
          fv = fu;
        }
      }
    }

    // Too many steps. Throw an exception.
    throw new ExtensionException(
            "Brent-minimize exceeded the maximum number of steps: "
            + maxEvaluations);
  }

}