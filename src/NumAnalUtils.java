package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Context;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.core.LogoList;
import org.nlogo.api.LogoListBuilder;

import Jama.Matrix;
import java.util.ArrayList;
import java.util.Arrays;

/**
 This class is package-private by default.
 */
class NumAnalUtils {

  static boolean linesearchFailed;

  protected static Matrix computeJacobian(Matrix X, Matrix F, AnonymousReporter fnctn,
          double epsilon, Context context)
          throws ExtensionException, LogoException {
    // Computes the forward-difference approximation to the computeJacobian at
    // point X of the set of functions contained in fnctn. F is a
    // vector of the function values at point X and the Jacobean is
    // returned as a Matrix.  Both X and F are 1xn Matrix objects.

    int n = X.getColumnDimension();
    Matrix J = new Matrix(n, n);

    // Compute the computeJacobian one column of partials at a time.
    for (int j = 0; j < n; j++) {
      // save the current value of X[j] and replace it with a value
      // a small distance, h, away.  Then compute the partial
      // derivatives and set X[j] back to its original value.
      double temp = X.get(0, j);
      double h = Math.max(epsilon * Math.abs(temp), epsilon);
      X.set(0, j, temp + h);
//          h = X.get(0, j) - temp;  // Trick to reduce finite precision error.
      Matrix F_new = getFofX(X, fnctn, context);
      Matrix Partial_j = (F_new.minus(F)).times(1.0 / h);

      // Partial_j is a 1xn Matrix, that needs to be inserted as the j'th
      // column of the computeJacobian.  So, we transpose it to an nx1 Matrix.
      J.setMatrix(0, n - 1, j, j, Partial_j.transpose());

      X.set(0, j, temp);
    }
    return J;
  }

  protected static Matrix performLineSearch(Matrix X_old, double ftest_old,
          Matrix nablaF, Matrix deltaX, Matrix X, double stpmax,
          AnonymousReporter fnctn, double tolx, double alpha, Context context)
          throws ExtensionException, LogoException {
    /*
     * Given an n-dimensional point X_old and, at that point, 
     * nablaF and deltaX, find and return a the new point in the 
     * direction given by deltaX where the function has decreased
     * "sufficiently". 
     * stpmax limits the length of the steps so that you do not try to
     * evaluate the function in regions where it is not defined or
     * subject to overflow. If the step succeeds, the global 
     * linesearchFailed is set to true. If the step fails, i.e., if the 
     * new point is too close to the old, linesearchFailed is set to 
     * false. In a minimization this usually signals
     * convergence and can be ignored. In root finding, the calling
     * program should check whether the convergence is spurious.
     * 
     * NOTE: deltaX is also changed in this procedure.  It is not 
     * explicitly returned, but rather modified "in place" using the 
     * .timesEquals() method.  This is a Java no-no, but we use it anyway.
     */

    int n = X.getColumnDimension();

    // Check step size and then compute slope = nablaF * (transpose)deltaX.
    double test = Math.sqrt((deltaX.times(deltaX.transpose())).get(0, 0));
    if (test > stpmax) {
      // Attempted step is too large. Scale down.
      // NOTE that deltaX is altered "in place", so its new values
      // will be seen in the calling program.
      deltaX = deltaX.timesEquals(stpmax / test);
    }
    double slope = (nablaF.times(deltaX.transpose())).get(0, 0);

    // Compute the minimum value for lambda and set the initial value of
    // lambda to one, a full Newton step.
    Matrix X_old1 = new Matrix(1, n);
    for (int j = 0; j < n; j++) {
      X_old1.set(0, j, Math.max(Math.abs(X_old.get(0, j)), 1.0));
    }
    test = NAMatrix.findMaxAbsElement(deltaX.arrayRightDivide(X_old1));
    double lambdaMin = tolx / test;
    double lambda = 1.0;

    // Don't really need to be initialized as they are not used in the
    // first step, but the compiler likes it.
    double tmplam = 0.0, lambda2 = 0.0, ftest2 = 0.0, ftest_old2 = 0.0;
    while (true) {

      if (lambda < lambdaMin) {
        // The change in the x vector required to get f to fall
        // significantly during backtracking has gotten very small.
        // We may have hit a local minimum, either at a zero root
        // or at some othe point. The calling program will need 
        // to check. We use the most recent X.
        linesearchFailed = true;
        return X;
      }

      // calculate the new point to be tried and try it.
      X = X_old.plus(deltaX.times(lambda));
      Matrix F = getFofX(X, fnctn, context);
      double ftest = 0.5 * (F.times(F.transpose())).get(0, 0);

      if (ftest <= ftest_old + alpha * lambda * slope) {
        // We've made a significant enough step toward the root
        // as measured by the decrease in ftest. Return the new X.
        linesearchFailed = false;
        return X;
      } else {
        // Backtrack.
        if (lambda == 1.0) {
          // This is the first backtrack.  Calculate a new trial
          // lambda value using a quadratic model for g(lambda).
          tmplam = -slope / (2.0 * (ftest - ftest_old - slope));
        } else {
          // The first backtrack did not work. Calculate a new
          // trial lambda using a cubic model for g(lambda).
          // NOTE, this could be set up as a matrix multiplication.
          // It might make it more transparent, but not necessarily
          // faster!
          double lambda_squared = lambda * lambda;
          double lambda2_squared = lambda2 * lambda2;
          double lambdaMlambda2 = lambda - lambda2;
          double rhs1 = ftest - ftest_old - lambda * slope;
          double rhs2 = ftest2 - ftest_old2 - lambda2 * slope;
          double a = (rhs1 / lambda_squared
                  - rhs2 / lambda2_squared) / lambdaMlambda2;
          double b = (-lambda2 * rhs1 / lambda_squared
                  + lambda * rhs2 / lambda2_squared)
                  / lambdaMlambda2;
          if (a == 0.0) {
            // coefficient on the cubic term is zero.
            tmplam = -slope / (2.0 * b);
          } else {
            double disc = b * b - 3.0 * a * slope;
            if (disc < 0.0) {
              throw new ExtensionException("Roundoff problem in LineSearch");
            } else {
              tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
            }
          }
          // constrain lambda <= 0.5*lambda1.
          tmplam = Math.min(0.5 * lambda, tmplam);
        }
      }
      lambda2 = lambda;
      ftest2 = ftest;
      ftest_old2 = ftest_old;

      // constrain lambda >= 0.1*lambda1.
      lambda = Math.max(tmplam, 0.1 * lambda);
    }
    // go back and try again.
  }

  public static double getFofX(double x, AnonymousReporter fnctn, Context context) {
    Object[] Array = {x};
    return (Double) fnctn.report(context, Array);
  }

  public static double getFofXvec(double[] x, AnonymousReporter fnctn, Context context) {
    LogoList xlist = convertArrayToSimpleLogoList(x);
    Object[] Array = {xlist};
    return (Double) fnctn.report(context, Array);
  }

   public static double [][] getMofXvec(double[] x, AnonymousReporter fnctn, 
           Context context) throws ExtensionException {
    LogoList xlist = convertArrayToSimpleLogoList(x);
    Object[] Array = {xlist};
    LogoList rslts = (LogoList) fnctn.report(context, Array);
    return NumAnalUtils.convertNestedLogoListToArray(rslts);
  }
   
public static Matrix getFofX(Matrix X, AnonymousReporter fnctn, Context context) {
    double[] x = X.getRowPackedCopy();
    LogoList xlist = NumAnalUtils.convertArrayToSimpleLogoList(x);
    Object[] Array = {xlist};
    LogoList rslts = (LogoList) fnctn.report(context, Array);
    double[] r = NumAnalUtils.convertSimpleLogoListToArray(rslts);
    return new Matrix(r, 1);
  }

  /*
   public static double GetFiOfXvec(double[] x,
   AnonymousReporter fnctn, int i, Context context) {
   // Calls the NetLogo reporter indicated by the task variable fntcn.
   // The reporter should take a list of x values and an index i, and
   // return the result of the ith equation evaluated at x.
   LogoList xlist = convertArrayToSimpleLogoList(x);
   Object[] Array = {xlist, (double) i};
   return (Double) fnctn.report(context, Array);
   }
   public static double[] GetFvecOfXvec(double[] x,
   AnonymousReporter fnctn, Context context) {
   // Calls for each value of i, 0 ... n, the NetLogo reporter
   // indicated by the task variable fntcn, and returns a
   // one-dimensional array of the function values evaluated at point x.
   int n = x.length;
   LogoList xlist = convertArrayToSimpleLogoList(x);
   Object[] Array = {xlist};
   LogoList rslts = (LogoList) fnctn.report(context, Array);
   return convertSimpleLogoListToArray(rslts);
   }
   */

  public static void writeToNetLogo(String mssg, Boolean toOutputArea, 
          Context context) throws ExtensionException {
    try {
      context.workspace().outputObject(mssg, null, true, true, 
              (toOutputArea) ? 
                      org.nlogo.api.OutputDestinationJ.OUTPUT_AREA() : 
                      org.nlogo.api.OutputDestinationJ.NORMAL());
    } catch (LogoException e) {
      throw new ExtensionException(e);
    }
  }

  public static void printMatrix(String lbl, Matrix M, Context context) 
          throws ExtensionException, LogoException {
    int m = M.getRowDimension();
    int n = M.getColumnDimension();
    NumAnalUtils.writeToNetLogo(lbl + " " + m + " " + n, false, context);
    for (int i = 0; i < m; i++) {
      StringBuilder buf = new StringBuilder();
      for (int j = 0; j < n; j++) {
        buf.append(M.get(i, j));
        buf.append(" ");
      }
      NumAnalUtils.writeToNetLogo(buf.toString(), false, context);
    }
  }

  public static double[][] convertNestedLogoListToArray(LogoList nestedLogoList) 
          throws ExtensionException {
    int numRows = nestedLogoList.size();
    if (numRows == 0) {
      throw new ExtensionException("input list was empty");
    }
    int numCols = -1;
    // find out the maximum column size of any of the rows,
    // in case we have a "ragged" right edge, where some rows
    // have more columns than others.
    for (Object obj : nestedLogoList.toJava()) {
      if (obj instanceof LogoList) {
        LogoList rowList = (LogoList) obj;
        if (numCols == -1) {
          numCols = rowList.size();
        } else if (numCols != rowList.size()) {
          throw new ExtensionException("To convert a nested list "
                  + "into a matrix, all nested lists must be the "
                  + "same length -- e.g. [[1 2 3 4] [1 2 3]] is "
                  + "invalid, because row 1 has one more entry.");
        }
      } else {
        throw new ExtensionException("To convert a nested list into "
                + "a matrix, there must be exactly two levels of "
                + "nesting -- e.g. [[1 2 3] [4 5 6]] creates a good "
                + "2x3 matrix.");
      }
    }
    if (numCols == 0) {
      throw new ExtensionException("input list contained only empty lists");
    }
    double[][] array = new double[numRows][numCols];
    int row = 0;
    for (Object obj : nestedLogoList.toJava()) {
      int col = 0;
      LogoList rowList = (LogoList) obj;
      for (Object obj2 : rowList.toJava()) {
        if (obj2 instanceof Number) {
          array[row][col] = ((Number) obj2).doubleValue();
          col++;
        } else {
          throw new ExtensionException("input list contains a non-number");
        }
      }
//      This should be unnecessary since we've checked to see that all 
//      columns are of equal length and only contain numbers.
//      // pad with zeros if we have a "ragged" right edge
//      for (; col < numCols; col++) {
//        array[row][col] = 0.0;
//      }
      row++;
    }
    return array;
  }
   
   public static LogoList convertArrayToNestedLogoList(double[][] dArray) {
    LogoListBuilder lst = new LogoListBuilder();
    for (double[] row : dArray) {
      LogoListBuilder rowLst = new LogoListBuilder();
      for (double elem : row) {
        rowLst.add(elem);
      }
      lst.add(rowLst.toLogoList());
    }
    return lst.toLogoList();
  }
  public static double[] convertSimpleLogoListToArray(LogoList xlist) {
    int n = xlist.size();
    double[] x = new double[n];
    for (int j = 0; j < n; j++) {
      x[j] = ((Number) xlist.get(j)).doubleValue();
    }
    return x;
  }
  
public static int[] convertSimpleLogoListToIntArray(LogoList xlist) {
    int n = xlist.size();
    int[] x = new int[n];
    for (int j = 0; j < n; j++) {
      x[j] = ((Number) xlist.get(j)).intValue();
    }
    return x;
  }

public static boolean[] convertSimpleLogoListToBooleanArray(LogoList xlist) {
    int n = xlist.size();
    boolean[] x = new boolean[n];
    for (int j = 0; j < n; j++) {
      x[j] = (Boolean) xlist.get(j);
    }
    return x;
  }

  public static LogoList convertArrayToSimpleLogoList(double[] x) {
    LogoListBuilder xlist = new LogoListBuilder();
    for (double elem : x) {
      xlist.add(elem);
    }
    return xlist.toLogoList();
  }

  public static void printValue(String lbl, double val, Context context) 
          throws ExtensionException, LogoException {
    NumAnalUtils.writeToNetLogo(lbl + " " + val, false, context);
  }
  
  public static void printBoolean(String lbl, boolean val, Context context) 
          throws ExtensionException, LogoException {
    NumAnalUtils.writeToNetLogo(lbl + " " + val, false, context);
  }
 
  public static void printInt(String lbl, int val, Context context) 
          throws ExtensionException, LogoException {
    NumAnalUtils.writeToNetLogo(lbl + " " + val, false, context);
  }
  
    public static void printString(String strng, Context context) {
        try {
            NumAnalUtils.writeToNetLogo(strng, false, context);
        } catch (ExtensionException e) {}
    }
  
  public static void printArray(String lbl, double[] vec, Context context)
          throws ExtensionException, LogoException {
      int n = vec.length;
      StringBuilder buf = new StringBuilder();
      for (int j = 0; j < n; j++) {
          buf.append(vec[j]);
          buf.append(" ");
      }
      NumAnalUtils.writeToNetLogo(lbl + " " + buf.toString(), false, context);
  }

  public static int findPattern(ArrayList<double[]> history, Context context) 
          throws ExtensionException, LogoException {
    // check for a pattern in the history array and return the length of
    // the pattern if one is found, otherwise zero.
    // first, see if the most recent entry is a repeat of an earlier one.
    // since the history lines are added at the "end" of the arraylist, we 
    // need to work backwards from the last line.
      
    int histSize = history.size();
//NumAnalUtils.writeToNetLogo("history "+histSize, false, context);    
    // there must be at least four lines to form a pattern as we don't allow
    // two identical contiguous lines to qualify.
    if (histSize < 2) { return 0; }
    
    double[] mostRecent = history.get(0);
    // see if the last entry is duplicated in the list.  we want the last
    // (most recent) duplicate that is not the last entry itself nor the one
    // immediately prior to the last entry.
    int dupIndex =0;
    boolean found = false;
    for (int i = 1; i < history.size() && !found; i++) {
      if (Arrays.equals(mostRecent, history.get(i))) {
        dupIndex = i;
        found = true;
      }
    }
    if (dupIndex == 0) { return 0; }
NumAnalUtils.writeToNetLogo("Duplicate found: "+dupIndex, false, context);
    // dupIndex points to the most recent duplicate of the last entry.
    // the length of the potential pattern is thus  also dupIndex.
    // to have a pattern, the arraylist must be at least twice the size of
    // potential pattern.
    int patternLength = dupIndex;
    if (patternLength * 2 > history.size()) { return 0; }
    // finally, check to see if we really have a pattern by comparing the 
    // two sublists that might potentially be duplicates. 
    for (int i = 0, j = dupIndex; i < patternLength; i++, j++) {
      if (!Arrays.equals(history.get(i), history.get(j))) {
        return 0;
      }
    }
    // the sublists match.
NumAnalUtils.writeToNetLogo("Pattern found: "+dupIndex, false, context);
    return dupIndex;
  }
}


