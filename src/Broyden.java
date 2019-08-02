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

public class Broyden {

  /* These procedures implement the Broyden algorithm for finding the root
   * (the "zero") of a set of n nonlinear equations in n variables. The 
   * central method is BroydenFindRoot.  It takes an initial guess for 
   * the n-dimentsional vector of inputs, x, passed as a NetLogo list, 
   * and a NetLogo AnonymousReporter, fnctn, referencing the set of equations 
   * for which the root is to be found, and returns the vector x (again,
   * as a NetLogo list) which yields the root.
   * 
   * Specifically, the NetLogo reporter, fnctn, should take a list of x 
   * values and should then return a list of the value of each equation
   * evaluated at x.
   * 
   * The following parameters determine the accuracy with which the Broyden
   * procedure finds the root. We are looking for the root of the series of
   * eqations (functions), and so one way of knowing when we've found it is
   * to look at the devation from zero of each the function values.  In 
   * particular, if F is the vector of results for any given input vector, x,
   * we look at f = 0.5*F*F, that is half of the sum of squared values of 
   * the function results. If this is small enough, i.e., less than tolf, 
   * we've found the root.
   * 
   * On the other hand, it is possible that we instead find a local or 
   * global minimum of f before we find the root.  In that case, the 
   * change in x that is required to get f to fall will get smaller and 
   * smaller as we approach that minimum.  The Broyden procedure also 
   * checks for that by seeing if the maximum proportionate change in 
   * the elements of x falls below tolx, or if the largest gradient in 
   * any of the x directions falls below tolmin.  Of course, it is possible
   * that the minimum is actually at the root, so the calling program
   * may want to check on that.
   * 
   * Default values for tolf, tolx and tolmin are given below. Experimentation 
   * suggests that tolx and tolmin be a couple of order of magnitudes smaller 
   * than tolf.
   * 
   * Other parameters are:
   * 
   * max_its: the maximum number of iterations (steps) allowed.
   * PRECISION: the square root of the machine precision for doubles.
   * stpmx: a parameter in the calculation of the maximum  step size
   *   in doLineSearch.
   * epsilon: the proportional change in each x element that is used to 
   *   calculate the computeJacobian (with a minimum change of epsilon itself).
   * alpha: a parameter that ensures that the doLineSearch routine has
   *   been able to find a new x vector that reduces f by a suffient 
   *   amount.
   * 
   * BroydenFailed returns true if the prior 
   * call to BroydenFindRoot resulted in a "soft" error, i.e., if we seem
   * to have found a local or global minimum, or false if a true root seems
   * to have been found.
   * 
   * 
   * This procedure is modeled on Numerical Recipes in C, but with a number 
   * of significan changes.  In terms of form, the algorithm has been 
   * recast to use matrix arithmetic.  This considerably simplifies the 
   * programming and makes the code much more transparent.  There is 
   * presumably a cost to this in terms of computational efficiency, but
   * it is likely to be very small.  We use the Jama Matrix class for
   * manipulating the matrices.  In terms of method, we do not use the 
   * QR decomposition of the computeJacobian, rather we use the computeJacobian itself
   * to solve for deltaX using the standard Jama Matrix.solve method.  This 
   * is likely a bit slower than using QR decomposition, but it makes for
   * much more transparent code since the Broyden update can be taken 
   * straight from any matrix description of the Broyden algorithm.
   */
  static final int MAX_ITS_DEFAULT = 1000;
  static final double PRECISION = Math.sqrt(Double.MIN_NORMAL);
  static final double TOLF_DEFAULT = 1.0e-6;
  static final double TOLX_DEFAULT = 1.0e-8;
  static final double STPMX_DEFAULT = 100.0;
  static final double TOLMIN_DEFAULT = 1.0e-8;
  static final double EPSILON_DEFAULT = 1.0e-4;
  static final double ALPHA_DEFAULT = 1.0e-6;

  public static class BroydenFailed implements Reporter {

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
    
    private static void setFlag (boolean flag) {
      failedToFindRoot = flag;
    }
  }

  public static class BroydenFindRoot implements Reporter {

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

      int nargs = args.length;
      int max_its = MAX_ITS_DEFAULT;
      double tolf = TOLF_DEFAULT;
      double tolx = TOLX_DEFAULT;
      double stpmx = STPMX_DEFAULT;
      double tolmin = TOLMIN_DEFAULT;
      double epsilon = EPSILON_DEFAULT;
      double alpha = ALPHA_DEFAULT;
      
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
        BroydenFailed.setFlag(false);
        return NumAnalUtils.convertArrayToSimpleLogoList(x);
      }

      // Calculate the maximum step for line searches.
      double sum = (X.times(X.transpose())).get(0, 0);
      double stpmax = stpmx * Math.max(Math.sqrt(sum), (double) n);

      // The compiler doesn't like it if these are not itialized.
      Matrix J = new Matrix(n, n);
      Matrix X_old = new Matrix(1, n);
      Matrix F_old = new Matrix(1, n);

      boolean restart = true;

      for (int its = 0; its < max_its; its++) {
        if (restart) {
          // Initialize or reintitialize the computeJacobian matrix, J.
          // Then do a QR decomposition of it and extract the parts.  
          J = NumAnalUtils.computeJacobian(X, F, fnctn, epsilon, context);
        } else {
          // Carry out a Broyden update, calculating the updated
          // computeJacobian in steps.
          // Note, we never get here on the first step, so X_old &
          // F_old will have been set on the previous step.
          // deltaX, deltaF are 1xn.  J is nxn.  So temp is nx1.
          // So temp*deltaX is an nxn outerproduct.
          Matrix dX = X.minus(X_old);
          Matrix dF = F.minus(F_old);
          Matrix temp
                  = dF.transpose().minus(J.times(dX.transpose()));
          // Here we need to make sure that at least one element of
          // temp has the required magnitude as we don't want to 
          // update with "noisy" elements. Elements that are too
          // small are set to zero.  Note that if all are zero, we
          // skip the update of J and go into the next line search
          // with the old J, but a new value of X.
          //
          // First create tmp as the matrix sum of the absolute 
          // values of F and F_old, then 
          // multiply each resulting element by the machine 
          // precision. Finally, if any element of temp is less than
          // the corresponding element of tmp, set that element of
          // temp to zero.
          Matrix tmp
                  = (NAMatrix.formAbsMatrix(F)).plus(NAMatrix.formAbsMatrix(F_old));
          tmp = tmp.timesEquals(PRECISION);
          for (int i = 0; i < n; i++) {
            if (temp.get(i, 0) < tmp.get(0, i)) {
              temp.set(i, 0, 0.0);
            }
          }
          if (NAMatrix.findMaxAbsElement(temp) > 0.0) {
            // At least one element of temp is big enough. 
            // Continue with the update of J.
            double dFsquared
                    = dF.times(dF.transpose()).get(0, 0);
            temp = temp.timesEquals(1.0 / dFsquared);
            J = J.plusEquals(temp.times(dX));
          }
        }

        // Compute nablaF = F*J = (Q*R)transpose * F = F*(Q*R), 
        // for the line search.
        Matrix nablaF = F.times(J);

        // Save old values of X and ftest.
        X_old = X.copy();
        F_old = F.copy();
        double ftest_old = ftest;

        // Solve J*deltaX' = -F' for deltaX.  Remember that J = Q*R 
        // (or at least approximately so after a Broyden update).
        Matrix deltaX
                = (J.solve(F.times(-1.0).transpose())).transpose();

        // Now set up the doLineSearch and, upon return, recalculate
        // F and ftest at the new point (unless the search failed!).
        // NOTE: doLineSearch actually returns three values.  See the 
        // comments associated with the doLineSearch, below.
        X = NumAnalUtils.performLineSearch(X_old, ftest_old, nablaF,
                deltaX, X, stpmax, fnctn, tolx, alpha, context);
        if (!NumAnalUtils.linesearchFailed) {
          F = NumAnalUtils.getFofX(X, fnctn, context);
          ftest = 0.5 * (F.times(F.transpose())).get(0, 0);
        }

        // Check for convergence and, if found, return the point.
        if (NAMatrix.findMaxAbsElement(F) < tolf) {
          BroydenFailed.setFlag(false);
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
          /* The line search failed to find a new x.
           * If we have already tried reinitializing the computeJacobian,
           * there is nothing left to do and we return the current
           * x.  Otherwise check for a gradient of zero signalling
           * a local minimmum. If we find it, we return the
           * current x.  In both cases we set failedToFindRoot to
           * true so the calling program can check it. 
           * If we have not just calculated the Jocobian and if 
           * the gradient is not zero, we signal a reinitialization of
           * the computeJacobian on the next step to try again.
           */
          if (restart) {
            BroydenFailed.setFlag(true);
            x = X.getRowPackedCopy();
//            NumAnalExtension.writeToNetLogo(
//                   "doLineSearch failed after restart.", false, context);
            return NumAnalUtils.convertArrayToSimpleLogoList(x);
          } else {
            // Check for a gradient of zero, i.e., spurious 
            // convergence.
            double den = Math.max(ftest, 0.5 * n);
            double test
                    = NAMatrix.findMaxAbsElement(nablaF.arrayTimes(X1));
            if (test / den < tolmin) {
              BroydenFailed.setFlag(true);
              x = X.getRowPackedCopy();
//              NumAnalExtension.writeToNetLogo(
//                "NewtonFindRoot: spurious convergence", false, context);
              return NumAnalUtils.convertArrayToSimpleLogoList(x);
            }
            restart = true;
          }

        } else {
          // doLineSearch did find a new x. Get set for using a Broyden
          // update on the next step, testing first for convergence.
          double test
                  = NAMatrix.findMaxAbsElement(X.minus(X_old).arrayRightDivide(X1));
          if (test < tolx) {
            BroydenFailed.setFlag(false);
            x = X.getRowPackedCopy();
            return NumAnalUtils.convertArrayToSimpleLogoList(x);
          }
          restart = false;
        }
      }
      // Out of the loop. Throw an error.
      throw new ExtensionException(
              "Broyden error - maximum number of iterations exceeded.");
    }
  }
}

  