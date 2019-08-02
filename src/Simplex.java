package org.nlogo.extensions.numanal;

// import java.util.Arrays;
import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.core.LogoList;

/*
 * Simplex finds the minimum of a multivariate function passed to it in 
 * the task variable, fnctn, and returns list of the input values at
 * the minimum.  fnctn should be a NetLogo reporter that takes
 * a list of input values and reports the value of the function at that
 * point. Simplex looks to the Bounds class to see if variable-by-variable
 * lower and/or upper bounds on the solution have been set by the user, 
 * and incorporates those bounds if they have.
 * 
 * These routines derive from the description of the Simplex method in Numerical
 * Recipes in C, but have been altered substantially, with input from other 
 * sources.
 */
public class Simplex {

  /* The Simplex class contains the global parameters and the subclasses
   * that actually do the work.
   * 
   * simplex requires two inputs:
   * guess	is a NetLogo list containing the initial guess for the 
   *  point (input values) at which the function is minimized.  
   *  Simplex converges pretty rapidly to a minimum from any initial 
   *  guess, but if there are local minima, it may get caught at one.  
   *  Thus entering a guess that is closer to the global minimum is 
   *  more likely to avoid the problem of getting stuck at a local minimum.
   * fnctn	is a task variable that refers to a NetLogo reporter 
   * evaluating the function to be minimized.  The reporter itself should 
   * take a single argument, a NetLogo list of input values, that is the 
   * coordinates of the n-dimensional point at which the function is to 
   * be evaluated, and report the value of the function at that point.
   *
   * The following optional inputs require that the primitive and its 
   * inputs be enclosed in parentheses.  Not all need be included, but 
   * all those to the left of the last one to be set must be.
   * delta	is the initial the amount by which each element of the 
   * initial guess is perturbed to form the simplex.  A value that 
   * is roughly 10% of the absolute value of the smallest element of 
   * the guess is not a bad place to start, although larger values 
   * often work quite nicely.  Larger values may lead to faster 
   * convergence, but may also lead to overshooting.  If delta is zero, 
   * the default value (10.0) is used.
   * rtol	is the relative tolerance to which the solution is to be found.  
   * If any step reduces the value of fnctn by a proportion less than rtol, 
   * we assume that the minimum has been found.  If rtol is negative, the 
   * default value (10e-12) is used.  rtol may be zero as long as atol 
   * is not also zero.
   * atol	is the absolute tolerance to which the solution is to be found.  
   * If any step reduces the value of fnctn by an amount less than atol, we 
   * assume that the minimum has been found.  If atol is negative, the 
   * default value (10e-12) is used.  atol may be zero as long as rtol is 
   * not also zero.
   * maxevals	insures that the procedure will not continue in an 
   * infinite loop if there is no convergence.  This sets the maximum 
   * number of evaluations of the function to be minimized and throws an 
   * ExtensionException if that number is exceeded. In the case of a 
   * restart, the number of function evaluations is reset to zero.  If 
   * maxevals is zero, the default (10,000) is used.
   * nrestarts	specifies the desired number of restarts of the simplex
   * procedure.  Many sources suggest a restart after the initial solution 
   * is found as the initial solution may be a false minimum. This, of 
   * course, will require more iterations, but given that it begins at 
   * the putative minimum, it should not require too many.  Note that 
   * the user can specify more than one restart, although it is not clear 
   * that there is any benefit for doing so.  The default number is zero.
   * nevalsmod	when the number of evaluations reaches an approximate 
   * multiple of this number, both the relative and absolute tolerances 
   * are increased by a factor of tolfactor.  This allows the routine 
   * to relax the tolerance required for a solution if the number of function 
   * evaluations grows too large.  By default, nevalsmod is set to 
   * (maxevals + 1) so that the tolerance is not changed.
   * tolfactor	is the factor by which the tolerances are multiplied 
   * each nevalsmod evaluations.  tolfactor must be >= 1.0.  Its default 
   * value is 2.0.
   */
  static final double SIDE_LENGTH_DEFAULT = 10.0;
  static final int NRESTARTS_MAX_DEFAULT = 0;
  static final int NEVALS_MAX_DEFAULT = 10000;
  static final int NEVALS_MOD_DEFAULT = NEVALS_MAX_DEFAULT + 1;
  static final double NEVALS_TOLFACTOR_DEFAULT = 2.0;
  static final double PRANGE_DEFAULT = 0.50;


  /*
   * We took out the random element in creating the initial simplex as
   * with it in, NetLogo will not give reproducable results when the 
   * NetLogo random-seed is set. It appears that this is so even when 
   * the Java Random class is constructed with a fixed seed.
   */
  // This method sets up the call to performSimplex.
  public static class SimplexSolve implements Reporter {

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

      double[] lwrBounds = null;
      double[] uprBounds = null;
      double sideLength = SIDE_LENGTH_DEFAULT;
      double relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
      double absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
      int nrestarts_max = NRESTARTS_MAX_DEFAULT;
      int nevals_max = NEVALS_MAX_DEFAULT;
      int nevals_mod = NEVALS_MOD_DEFAULT;
      double nevals_tolfactor = NEVALS_TOLFACTOR_DEFAULT;

      // Get the initial guess and turn it from a LogoList to a vector.
      LogoList xlist = args[0].getList();
      int nvar = xlist.size();
      double[] x = NumAnalUtils.convertSimpleLogoListToArray(xlist);

      // If this is a constrained simplex, check the bounds.
      if (Bounds.lowerBounds != null) {
        // There are bounds. Check that they are of the right
        // dimension and that the initial point is within them.
        // If all is okay, save them.
        if (Bounds.lowerBounds.length == nvar) {
          for (int i = 0; i < nvar; i++) {
            if (x[i] < Bounds.lowerBounds[i] || x[i] > Bounds.upperBounds[i]) {
              throw new ExtensionException("The initial guess "
                      + "must be within the specified bounds.");
            }
          }
          lwrBounds = Bounds.lowerBounds.clone();
          uprBounds = Bounds.upperBounds.clone();
        } else {
          // The bound dimentions do not match the problem.
          throw new ExtensionException(
                  "The dimension of the problem " + nvar
                  + " does not match the dimension of the bounds "
                  + Bounds.lowerBounds.length + ".");
        }
      }

      // Get reporter task.
      AnonymousReporter fnctn = args[1].getReporter();

      // Save the remainder of the arguments. If any is <= 0, use the
      // default for that parameter, except for the tolerances, one of
      // which can be zero. In that case, < 0 signifies the default.
      int nargs = args.length;
      if (nargs > 2) {
        double xx = args[2].getDoubleValue();
        sideLength = (xx > 0.0) ? xx : SIDE_LENGTH_DEFAULT;
      }
      if (nargs > 3) {
        double xx = args[3].getDoubleValue();
        relTolerance = (xx >= 0.0) ? xx : -1;
        // negative value means use absTolerance only.
      }
      if (nargs > 4) {
        double xx = args[4].getDoubleValue();
        absTolerance = (xx >= 0.0) ? xx : -1;
        // negative value means use relTolerance only.
      }
      if (nargs > 5) {
        int mm = args[5].getIntValue();
        nevals_max = (mm > 0) ? mm : NEVALS_MAX_DEFAULT;
        nevals_mod = nevals_max + 1;
      }
      if (nargs > 6) {
        nrestarts_max = args[6].getIntValue();  // the default is zero
      }
      if (nargs > 7) {
        int mm = args[7].getIntValue();
        nevals_mod = (mm > 0)
                ? Math.max(Math.min((nevals_max + 1), mm), 1)
                : nevals_max + 1;
      }
      if (nargs > 8) {
        double xx = args[8].getDoubleValue();
        nevals_tolfactor = (xx > 0.0) ? Math.max(1.0, xx)
                : NEVALS_TOLFACTOR_DEFAULT;
      }

      // Need to make sure both tols are not negative. (One is okay.)
      // If they are, assume user wanted defaults for both.
      if (relTolerance < 0.0 && absTolerance < 0.0) {
        relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
        absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
      }

      // Perform the simplex. Loop through the specified number of 
      // restarts, keeping track of the best solution and returning 
      // that solution at the end.
      double bestResult = Double.MAX_VALUE;
      double[] bestSolution = new double[nvar];
      for (int i = 0; i <= nrestarts_max; i++) {
        x = performSimplex(x, fnctn, context, relTolerance,
                absTolerance, sideLength, lwrBounds, uprBounds,
                nevals_max, nevals_mod, nevals_tolfactor);
        double thisResult = NumAnalUtils.getFofXvec(x, fnctn,
                context);
        if (thisResult < bestResult) {
          bestResult = thisResult;
          bestSolution = x.clone();
        }
      }
      return NumAnalUtils.convertArrayToSimpleLogoList(bestSolution);
    }
  }

  // This is the guts of the simplex algorithm.  It minimizes the 
  // function, fntn, to a given tolerance.  x is the initial guess of
  // a solution and sideLength is the amount by which each dimension of 
  // the initial guess is perturbed to form the initial simplex.
//  @SuppressWarnings("UnusedAssignment")
  private static double[] performSimplex(double[] x, AnonymousReporter fnctn,
          Context context, double relTol, double absTol, double sideLength,
          double[] lowerBounds, double[] upperBounds, int nevals_max,
          int nevals_mod, double nevals_tolfactor)
          throws ExtensionException, LogoException {

    int nvar = x.length;
    int nevals = 0;
    int nevalsPrior = 0;

    // Now construct the simplex matrix. Set each row to the intial vertex
    // and then perturb a different element of the second through last
    // rows by sideLength times a random factor in the range 1.0 plus or minus
    // PRANGE_DEFAULT/2, giving us nvar + 1 different vertices.
    // We have taken out the random element.
    int nvert = nvar + 1;
    double[][] s = new double[nvert][nvar];
    s[0] = (double[]) x.clone();
    for (int i = 1; i < nvert; i++) {
      s[i] = (double[]) x.clone();
//    s[i][i - 1] += ((1.0 - PRANGE_DEFAULT/2.0) + 
//            randgen.nextDouble()*PRANGE_DEFAULT) * sideLength;
      s[i][i - 1] += sideLength;
    }

    // y is a vector of results.  Fill it with the values of the
    // objective function for each vertex.
    double y[] = new double[nvert];
    for (int i = 0; i < nvert; i++) {
      y[i] = NumAnalUtils.getFofXvec(s[i], fnctn, context);
    }
    nevals += nvert;

    // sums is a vector of column sums of s.
    double sums[] = getColumnSums(s);

    // Now we enter into the solution loop.  Note that we break out of
    // the loop manually when a solution is found or when the maximum
    // number of function evaluations has been reached or exceeded.
    while (true) {
      // Find the vertices with the worst (highest), next worst
      // (next highest) and best (lowest) values of the function.
      int iwrst = 0;
      int i2wst = 1;
      if (y[1] > y[0]) {
        iwrst = 1;
        i2wst = 0;
      }
      int ibest = 0;
      for (int i = 0; i < nvert; i++) {
        if (y[i] <= y[ibest]) {
          ibest = i;
        }
        if (y[i] > y[iwrst]) {
          i2wst = iwrst;
          iwrst = i;
        } else if (y[i] > y[i2wst] && i != iwrst) {
          i2wst = i;
        }
      }

      /*
       * Now check to see if we have achieved the desired tolerance,
       * either relative or absolute, by comparing the values at the 
       * best and worst points.  NOTE, relTol or absTol may be negative
       * (but not both). Since diff is always positive, it will never be 
       * <= to the negative tolerance and thus that tolerance is ignored.
       * Other possibilities include looking at the size of
       * the simplex itself and stopping when it gets very small.
       */
      double diff = Math.abs(y[iwrst] - y[ibest]);
      double size = Math.max(Math.abs(y[iwrst]), Math.abs(y[ibest]));
      if ((diff <= size * relTol) || (diff <= absTol)) {
        return s[ibest];
      }

      // Check for an "infinite" loop.
      if (nevals > nevals_max) {
        throw new ExtensionException(
                "Simplex exceeded " + nevals_max
                + " iterations, the maximum number specified.");
      }

      // Check to see if the tolerance should be increased.  NOTE, this will
      // preserve a negative tolerance as nevals_tolfactor is always positive.
      if ((nevals / nevals_mod) > (nevalsPrior / nevals_mod)) {
        relTol *= nevals_tolfactor;
        absTol *= nevals_tolfactor;
        // and report this to the command center.
        NumAnalUtils.writeToNetLogo("At step " + nevals
                + " the relative and absolute tolerances were "
                + "increased by a factor of "
                + nevals_tolfactor + " to " + relTol + " and "
                + absTol, false, context);
      }
      nevalsPrior = nevals;

      // Time to actually begin the next iteration.
      // With the factor of -1, tryNewVertex reflects the current worst 
      // vertex through the opposite face of the simplex and, if the 
      // reflected vertex is better than the original, replaces the 
      // original vertex with the reflected one.  
      // However, if tryNewVertex finds that the reflected vertex
      // is not better than the original, the original vertex is not 
      // replaced, but the (even worse) value of the function at 
      // that reflected vertex is returned in order to signal the need 
      // for a contraction of from the original vertex, in the test 
      // below.
      double ynew = tryNewVertex(s, y, sums, fnctn, context, iwrst,
              -1.0, lowerBounds, upperBounds);
      nevals += 1;

      if (ynew <= y[ibest]) {
        // The reflected vertex is better than the current best vertex.
        // Expand the simplex in the direction of the new vertex to 
        // hopefully achieve an even better point.  Note that iwrst
        // is pointing to the new, reflected vertex, not to the worst
        // vertex in the new simplex. Note too that the return value of
        // ynew is not used in this case.
        tryNewVertex(s, y, sums, fnctn, context, iwrst,
                2.0, lowerBounds, upperBounds);
        nevals += 1;

      } else if (ynew >= y[i2wst]) {
        // If ynew is as bad or worse than the 2nd worst point in the
        // simplex, iwrst must still point to the worst vertex. 
        // That vertex may either be the original worst vertex, or 
        // the reflection of it (which apparently is only marginally
        // better).  In either case, contract the simplex from that 
        // vertex to hopefully find a better point.
        double ysave = y[iwrst];
        ynew = tryNewVertex(s, y, sums, fnctn, context, iwrst,
                0.5, lowerBounds, upperBounds);
        nevals += 1;

        if (ynew >= ysave) {
          // Well, even with the contraction away from the worst
          // vertex we did not find a better point. So, contract 
          // the entire simplex around the best point, moving 
          // each vertex halfway to the best vertex. (Note that
          // the if() statement is not strictly necessary as the 
          // best vertex would not move anyway. But the if() saves 
          // nvar calculations for a new vertex and one function
          // evaluation at the expense of nvert logical 
          // comparisons.)
          for (int i = 0; i < nvert; i++) {
            if (i != ibest) {
              for (int j = 0; j < nvar; j++) {
                s[i][j] = 0.5 * (s[i][j] + s[ibest][j]);
              }
              y[i] = NumAnalUtils.getFofXvec(s[i], fnctn,
                      context);
            }
          }
          nevals += (nvert - 1);
          // Update sums.
          sums = getColumnSums(s);
        }
      }  // End of the if, else if block.
    }  // End of the while loop.
  }  // End of performSimplex

  private static double[] getColumnSums(double[][] matrx) {
    // Returns an array containing the column sums of array, matrx.
    int nrows = matrx.length;
    int ncols = matrx[0].length;
    double[] colsums = new double[ncols];
    for (int j = 0; j < ncols; j++) {
      double sum = 0.0;
      for (int i = 0; i < nrows; i++) {
        sum += matrx[i][j];
      }
      colsums[j] = sum;
    }
    return colsums;
  }

  private static double tryNewVertex(double[][] s, double[] y, double[] sums,
          AnonymousReporter fnctn, Context context, int iwrst, double factor,
          double[] lowerBounds, double[] upperBounds)
          throws ExtensionException, LogoException {
    /*
     * If factor is positive and greater than one, the simplex is expanded
     * by pushing out the the vertex iwrst (which may no longer be the
     * worst!). If factor is positive and less than one, the simplex 
     * is contracted by pulling in the vertex iwrst. Finally, if factor 
     * negative, the vertex iwrst is reflected through the opposite 
     * face of the simplex with abs(factor) being the degree of 
     * reflection.
     */

    int nvar = sums.length;
    double factor1 = (1.0 - factor) / nvar;
    double factor2 = factor1 - factor;
    double ptry[] = new double[nvar];
    // First, find the new expanded, contracted or reflected vertex.
    for (int j = 0; j < nvar; j++) {
      ptry[j] = sums[j] * factor1 - s[iwrst][j] * factor2;
    }
    // Apply constraints, if any.
    // NumAnalUtils.printString("testing bounds." + " " + Arrays.toString(lowerBounds) + 
    //    " " + Arrays.toString(upperBounds) + " " + Arrays.toString(ptry), context);
    if (lowerBounds != null) {
      for (int j = 0; j < nvar; j++) {
        ptry[j] = Math.max(lowerBounds[j], ptry[j]);
        ptry[j] = Math.min(upperBounds[j], ptry[j]);
      }
    //  NumAnalUtils.printString("applied bounds." + " " + Arrays.toString(lowerBounds) + 
    //    " " + Arrays.toString(upperBounds) + " " + Arrays.toString(ptry), context);
    }

    double ynew = NumAnalUtils.getFofXvec(ptry, fnctn, context);
    // If the new vertex is better than the original one, replace the 
    // original one with the new one, updating column sums as well. If 
    // it is not better, the original vertex is not replaced.
    if (ynew < y[iwrst]) {
      y[iwrst] = ynew;
      for (int j = 0; j < nvar; j++) {
        sums[j] += ptry[j] - s[iwrst][j];
        s[iwrst][j] = ptry[j];
      }
    }

    // Return the value at the new vertex, whether or not it has 
    // replaced the original one.
    return ynew;
  }
}
