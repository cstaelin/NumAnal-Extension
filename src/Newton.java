package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.core.LogoList;

import Jama.Matrix;

public class Newton {

  /* These procedures implement the Newton algorithm for finding the root
   * (the "zero") of a set of n nonlinear equations in n variables. The 
   * central routine is NewtonFindRoot.  It takes an initial guess for 
   * the n-dimentsional vector of inputs, x, passed as a NetLogo list, 
   * and a NetLogo AnonymousReporter, fnctn, referencing the set of equations 
   * for which the root is to be found, and returns the vector x (again,
   * as a NetLogo list) which yields the root.
   * 
   * Specifically, the NetLogo reporter, fnctn, should take a list of x 
   * values and should then return in a list the results of each equation 
   * evaluated at x. 
   * 
   * The following parameters determine the accuracy with which the Newton
   * procedure finds the root. We are looking for the root of the series of
   * eqations (functions), and so one way of knowing when we've found it is
   * to look at the devation from zero of each the function values.  In 
   * particular, if F is the vector of results for any given input vector, X,
   * we look at f = 0.5*F*F', that is half of the sum of squared values of 
   * the function results. If this is small enough, i.e., less than tolf, 
   * we've found the root. (F' is the transpose of F.)
   * 
   * On the other hand, it is possible that we instead find a local or 
   * global minimum of f before we find the root.  In that case, the 
   * change in X that is required to get f to fall will get smaller and 
   * smaller as we approach that minimum.  The Newton procedure also 
   * checks for that by seeing if the maximum proportionate change in 
   * the elements of X falls below tolx, or if the largest gradient in 
   * any of the X directions falls below tolmin.  Of course, it is possible
   * that the minimum is actually at the root, so the calling program
   * may want to check on that.
   * 
   * Default values for tolf, tolx and tolmin are given below. Experimentation 
   * suggests that tolx and tolmin be a couple of magnitudes smaller 
   * than tolf.
   * 
   * Other parameters are:
   * 
   * max_its: the maximum number of iterations (steps) allowed.
   * PRECISION: the square root of the machine precision for doubles.
   * stpmx: a parameter in the calculation of the maximum  step size
   *   in performLineSearch.
   * epsilon: the proportional change in each x element that is used to 
   *   calculate the computeJacobian (with a minimum change PRECISION).
   * alpha: a parameter that ensures that the performLineSearch routine has
   *   been able to find a new x vector that reduces f by a suffient 
   *   amount.
   * 
   * 
   * NewtonFailed returns true if the prior 
   * call to NewtonFindRoot resulted in a "soft" error, i.e., if we seem
   * to have found a local or global minimum, or false if a true root seems
   * to have been found.
   * 
   * These routines derive from the description of the Newton method in Numerical
   * Recipes in C, but have been altered substatially to use use matrix
   * arithemetic.  This simplifies the code significantly and makes the
   * algorithm more transparent.  The matrices are defined and handled by
   * the Jama Matrix package.
   */
  static final int MAX_ITS_DEFAULT = 1000;
  static final double PRECISION = Math.sqrt(Double.MIN_NORMAL);
  static final double TOLF_DEFAULT = 1.0e-6;
  static final double TOLX_DEFAULT = 1.0e-8;
  static final double STPMX_DEFAULT = 100.0;
  static final double TOLMIN_DEFAULT = 1.0e-8;
  static final double EPSILON_DEFAULT = 1.0e-4;
  static final double ALPHA_DEFAULT = 1.0e-6;

  // Indicates whether the most recent call to NewtonFindRoot ended in 
  // failure. Even if true is returned, NewtonFindRoot may nevertheless 
  // have returned a point very close to the root, so it may be worth 
  // checking it.
  public static class NewtonFailed implements Reporter {

    static boolean failedToFindRoot;

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{}, Syntax.BooleanType());
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {
      return failedToFindRoot;
    }

    private static void setFlag(boolean flag) {
      failedToFindRoot = flag;
    }
  }

  // The main routine for finding the root.
  public static class NewtonFindRoot implements Reporter {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{Syntax.ListType(),
        Syntax.ReporterType(),
        Syntax.NumberType() | Syntax.RepeatableType()},
              Syntax.ListType(), 2);
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      int max_its = MAX_ITS_DEFAULT;
      double tolf = TOLF_DEFAULT;
      double tolx = TOLX_DEFAULT;
      double stpmx = STPMX_DEFAULT;
      double tolmin = TOLMIN_DEFAULT;
      double epsilon = EPSILON_DEFAULT;
      double alpha = ALPHA_DEFAULT;

      int nargs = args.length;
      // get the initial guess and put it in a 1xn Matrix, X
      LogoList xlist = args[0].getList();
      int n = xlist.size();
      double[] x = NumAnalUtils.convertSimpleLogoListToArray(xlist);
      Matrix X = new Matrix(x, 1);

      // Save the remainder of the arguments.
      AnonymousReporter fnctn = args[1].getReporter();

      if (nargs > 2) {
        double xx = args[2].getDoubleValue();
        stpmx = (xx > 0.0) ? xx : STPMX_DEFAULT;
      }
      if (nargs > 3) {
        double xx = args[3].getDoubleValue();
        tolf = (xx > 0.0) ? xx : TOLF_DEFAULT;
      }
      if (nargs > 4) {
        double xx = args[4].getDoubleValue();
        tolx = (xx > 0.0) ? xx : TOLX_DEFAULT;
      }
      if (nargs > 5) {
        double xx = args[5].getDoubleValue();
        tolmin = (xx > 0.0) ? xx : TOLMIN_DEFAULT;
      }
      if (nargs > 6) {
        int mm = args[6].getIntValue();
        max_its = (mm > 0) ? mm : MAX_ITS_DEFAULT;
      }
      if (nargs > 7) {
        double xx = args[7].getDoubleValue();
        epsilon = (xx > 0.0) ? xx : EPSILON_DEFAULT;
      }
      if (nargs > 8) {
        double xx = args[8].getDoubleValue();
        alpha = (xx > 0.0) ? xx : ALPHA_DEFAULT;
      }

      // Evaluate all the equations at the initial guess and put the
      // results in a 1xn matrix, F.
      // ftest is then 0.5 * F*F'.
      Matrix F = NumAnalUtils.getFofX(X, fnctn, context);
      double ftest = 0.5 * (F.times(F.transpose())).get(0, 0);

      // Test to see if the guess is a root, using a tougher test
      // than TOLF.  The test is based on the maximum deviation of any 
      // f(x) from zero.
      if (NAMatrix.findMaxAbsElement(F) < 0.01 * tolf) {
        NewtonFailed.setFlag(false);
        return NumAnalUtils.convertArrayToSimpleLogoList(x);
      }

      // Calculate the maximum step for line searches.
      double sum = (X.times(X.transpose())).get(0, 0);
      double stpmax = stpmx * Math.max(Math.sqrt(sum), (double) n);

      for (int its = 0; its < max_its; its++) {
        // Compute the computeJacobian. Put it into a Matrix object, then 
        // compute (nabla)F as F.J.
        Matrix J = NumAnalUtils.computeJacobian(X, F, fnctn, epsilon, context);
        Matrix nablaF = F.times(J);

        // Save old values of X and ftest.
        Matrix X_old = X.copy();
        double ftest_old = ftest;

        Matrix deltaX;
        deltaX = (J.solve(F.times(-1.0).transpose())).transpose();

        // Now set up the performLineSearch and, upon return, recalculate
        // F and ftest at the new point (unless the search failed!).
        // NOTE: performLineSearch actually returns three values.  See the 
        // comments associated with the performLineSearch, below.
        X = NumAnalUtils.performLineSearch(X_old, ftest_old, nablaF,
                deltaX, X, stpmax, fnctn, tolx, alpha, context);
        if (!NumAnalUtils.linesearchFailed) {
          F = NumAnalUtils.getFofX(X, fnctn, context);
          ftest = 0.5 * (F.times(F.transpose())).get(0, 0);
        }

        // Check for convergence and, if found, return the point.
        if (NAMatrix.findMaxAbsElement(F) < tolf) {
          NewtonFailed.setFlag(false);
          x = X.getRowPackedCopy();
          return NumAnalUtils.convertArrayToSimpleLogoList(x);
        }

        // We need this matrix in two places below, so form it now.
        // Each element of X1 is the larger of the absolute value of
        // the corresponding element of X or 1.
        Matrix X1 = new Matrix(1, n);
        for (int j = 0; j < n; j++) {
          X1.set(0, j, Math.max(Math.abs(X.get(0, j)), 1.0));
        }

        if (NumAnalUtils.linesearchFailed) {
          // Check for a gradient of zero, i.e., spurious 
          // convergence.
          double den = Math.max(ftest, 0.5 * n);
          double test
                  = NAMatrix.findMaxAbsElement(nablaF.arrayTimes(X1));
          if (test / den < tolmin) {
            NewtonFailed.setFlag(true);
            x = X.getRowPackedCopy();
            NumAnalUtils.writeToNetLogo(
                    "NewtonFindRoot: spurious convergence",
                    false, context);
            return NumAnalUtils.convertArrayToSimpleLogoList(x);
          }
        }

        // Test for convergence on deltaX, based on the maximum
        // proportionate change in an x value.
        double test
                = NAMatrix.findMaxAbsElement(X.minus(X_old).arrayRightDivide(X1));
        if (test < tolx) {
          NewtonFailed.setFlag(false);
          x = X.getRowPackedCopy();
          return NumAnalUtils.convertArrayToSimpleLogoList(x);
        }
      }

      // Out of the loop. Throw an error.
      throw new ExtensionException(
              "Newton error - maximum number of iterations exceeded.");
    }
  }
}
