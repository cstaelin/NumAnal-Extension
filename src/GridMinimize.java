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

import java.util.Arrays;
import java.util.ArrayList;

public class GridMinimize implements Reporter {

  /* Ranges across a hyper-grid of points in n dimentions and evaluates
   the function at each point on the hyper-grid, returning the cordinates
   of the point yeilding the lowest value of the function. Once a minimum is 
   found, the routine constructs a new, finer grid around that minimum and 
   finds the minimum point within that new grid. The process continues until
   the maximum interval between grid points in every dimension is less than
   or equal to tol, the tolerance.
   Takes five arguments:
   1. the function as an anonymous reporter.
   2. a list of length n containing the lower bound of the range of values in
   each dimension.
   3. a list of length n containing the upper bound of the range of values in 
   each dimension.
   4. the maximum number of evaluatins of the function. This, along with the 
   range of values in each dimension (upper - lower) will determine the 
   fineness of the grid by setting the number of intervals in each dimension.
   5. the tolerance 
   */
  private class CurrentMin {

    // holds the current minimum value and the point at which it occurs.

    double value;
    double[] point;

    CurrentMin(int n) {
      value = Bounds.UPPER_BOUND_DEFAULT;
      point = new double[n];
    }
  }

  private final double dmv = 0.50;

  @Override
  public Syntax getSyntax() {
    return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
      Syntax.ListType(), Syntax.ListType(),
      Syntax.NumberType(), Syntax.NumberType()}, Syntax.ListType(), 5);
  }

  @Override
  public Object report(Argument args[], Context context)
          throws ExtensionException, LogoException {

    // get the four arguments.
    AnonymousReporter fnctn = args[0].getReporter();
    LogoList lowerList = args[1].getList();
    LogoList upperList = args[2].getList();
    int maxEvals = args[3].getIntValue();
    double tol = args[4].getDoubleValue();

    double[] lower = NumAnalUtils.convertSimpleLogoListToArray(lowerList);
    double[] upper = NumAnalUtils.convertSimpleLogoListToArray(upperList);

    int n = lower.length;
    if (n != upper.length) {
      throw new ExtensionException("GridMinimize: the lower and upper "
              + "lists are not of the same length.");
    }
    for (int i = 0; i < n; i++) {
      // make sure lower < upper.
      if (lower[i] > upper[i]) {
        double tmp = lower[i];
        lower[i] = upper[i];
        upper[i] = tmp;
      }
    }

    // call findGridMin to set up the grid and find the point at which the 
    // function is at its minimum.
    // findGridMin returns the minimum value and the point at which it was
    // found in currentMin. GridMinimize returns only the point as a list.
    CurrentMin currentMin = new CurrentMin(n);
    findGridMin(lower, upper, maxEvals, tol, currentMin, fnctn, context);

    return NumAnalUtils.convertArrayToSimpleLogoList(currentMin.point);
  }

  private void findGridMin(double[] lower, double[] upper, int maxEvals,
          double tol, CurrentMin currentMin, AnonymousReporter fnctn,
          Context context) throws ExtensionException, LogoException {

    int n = lower.length;
    // Construct the grid array using the range and the number of intervals
    // in each dimension.
    // First, allocate the maximum number of function evaluations that we 
    // have set among the markets in proportion to the ranges of each to 
    // keep the increment in each range roughly equal.  
    // But, insure that no interval falls below tol.
    double rng[] = new double[n];
    for (int i = 0; i < n; i++) {
      rng[i] = upper[i] - lower[i];
    }
    double r0 = rng[0];
    double x = Arrays.stream(rng).reduce(1, (a, b) -> a * (b / r0));
    double nInt1 = Math.pow((maxEvals / x), (1.0 / n));
    int[] numIntervals = new int[n];
    // By truncating the double to an int, rathr than rounding, we 
    // may end up with fewer grid points than specified by maxEvals.
    // This seemed safer than having more as we might if we rounded.
    // We also make sure that there are at least three invervals in each
    // dimension (see below), but that no interval falls below tol.
    for (int i = 0; i < n; i++) {
      numIntervals[i] = (int) ((rng[i] / rng[0]) * nInt1);
      if (rng[i] / numIntervals[i] < tol) {
        numIntervals[i] = (int) Math.ceil(rng[i] / tol);
      } else {
        numIntervals[i] = Math.max(3, numIntervals[i]);
      }
    }

    // Next construct a 2D array of grid points. Eash line is a dimension
    // of the grid and each column a value to evaluate the function at
    // in that dimension. Since different dimensions might have different
    // numbers of increments and thus values, this could be a ragged array.
    double ranges[][] = new double[n][];
    for (int i = 0; i < n; i++) {
      // Construct the row for this dimension. Fill in ranges[i][] with
      // values from lower[i] to upper[i], inclusive.
      ranges[i] = new double[numIntervals[i] + 1];
      double incrmnt = rng[i] / numIntervals[i];
      for (int j = 0; j < (numIntervals[i] + 1); j++) {
        ranges[i][j] = lower[i] + j * incrmnt;
      }
    }

    // We will keep a history of solutions to make sure that we are not
    // cycling.
    ArrayList<double[]> history = new ArrayList<>();
    
    // Now call the findMin method which will work recusively through 
    // the grid, keeping the curren minimum and its point in the 
    // CurrentMin structure.
    findMin(0, n, ranges, lower, currentMin, fnctn, context);
    history.add(0, currentMin.point.clone());
//NumAnalUtils.writeToNetLogo("Added-1 "+Arrays.toString(currentMin.point)+Arrays.toString(lower)+Arrays.toString(upper)+Arrays.toString(numIntervals)+currentMin.value, false, context);

    // So, we've found the minimum of the points in this grid.  Look 
    // first to see if it is at the boundary of any dimension.  If so, 
    // we move the grid in that direction to see if a point outside the 
    // current grid is better. The grid is moved by dmv*the size of the 
    // current range of values in that dimension. A dmv of 0.5, for
    // instance, would move the range from 2 -> 6 to 4 -> 8. dmv must be 
    // less than 1.0 or else there will be a gap between the new range 
    // the old that could conceivably contain the point we are looking for.
    // dmv = 0.5 seems to work well.

    boolean doAgain = true;
    while (doAgain) {
      // check ech dimension to see if the solution is at a boundary.
      // if so, move that range in the appropriate direction.
      boolean movedOne = false;
      for (int i = 0; i < n; i++) {
        double moveBy = 0.0;
        if (currentMin.point[i] <= lower[i]) {
          moveBy = -dmv * rng[i];
        }
        else if (currentMin.point[i] >= upper[i]) {
          moveBy = dmv * rng[i];
        }
        if (moveBy != 0.0) {
          // it is at a boundary. move the bounds and the range for this
          // dimension up or down depending on the sign of moveBy.
          lower[i] += moveBy;
          upper[i] += moveBy;
          for (int j = 0; j < ranges[i].length; j++) {
            ranges[i][j] += moveBy;
          }
          movedOne = true;
//NumAnalUtils.writeToNetLogo("moved "+i+" by "+moveBy, false, context);
        } 
      }
      if (movedOne) {
        // at least one dimension's bounds have been moved. Before solving 
        // with the new range(s), 
        // Check to see whether in the process of moving bounds we are 
        // cycling by seeing if there is a repeating
        // pattern in the history. If so, average the points in the pattern
        // and use that as the solution.
        // (We could also use this point to begin a new search on this grid,
        // but would that simply lead to a new cycling? It might if the averaged
        // point is the same as one of the points in the pattern, or if using it
        // simply resulted in a new pattern of cycling. As it is, we just
        // accept the average.  Think about this.)
        int p = NumAnalUtils.findPattern(history, context);
        if (p != 0) {
//NumAnalUtils.writeToNetLogo("pattern found of length "+p, false, context);
          // we've found a pattern of length p. average the points in the pattern.
          double[] avgPoint = new double[n];
          Arrays.fill(avgPoint, 0);
          for (int j = 0; j < n; j++) {
            for (int i = 0, lineNum = history.size() - 1; i < p; i++, lineNum--) {
              avgPoint[j] += history.get(lineNum)[j];
            }
          }
          avgPoint = Arrays.stream(avgPoint).map(a -> a / n).toArray();
          // now set bounds around this new point to get new values of 
          // lower and upper, using the same ranges in each dimension.
          for (int i = 0; i < n; i++) {
            lower[i] = avgPoint[i] - rng[i] / 2;
            upper[i] = avgPoint[i] + rng[i] / 2;
          }
          // put this point and its value in currentMin.
          System.arraycopy(avgPoint, 0, currentMin.point, 0, n);
          currentMin.value = NumAnalUtils.getFofXvec(avgPoint, fnctn, context);
          history.add(0, currentMin.point.clone());
//NumAnalUtils.writeToNetLogo("Added-2 "+Arrays.toString(currentMin.point)+Arrays.toString(rng)+currentMin.value, false, context);
          doAgain = false;
        }
        else {
          // no pattern (yet).  keep going with new ranges.
          findMin(0, n, ranges, lower, currentMin, fnctn, context);
          history.add(0, currentMin.point.clone());
//NumAnalUtils.writeToNetLogo("Added-3 "+Arrays.toString(currentMin.point)+Arrays.toString(lower)+Arrays.toString(upper)+Arrays.toString(numIntervals)+currentMin.value, false, context);
//          for (double[] entry : history) {
//            NumAnalUtils.writeToNetLogo(Arrays.toString(entry), false, context);
//          }
        }
      }
      else {
        // we've found an internal solution.
        doAgain = false;
      }
    }

    // We've finished with this grid. Now see if we need to subdivide it to
    // achieve the desired precision as given by tol. Note that we do this by
    // shrinking the range for any dimension that has a currentInterval > tol.
    // When we call the procedure recursively, the next interation should 
    // result in smaller intervals given that maxEvals does not change. 
    // Indeed, since we reduce the range to two currentIntervals and make sure
    // above that there are at least three intervals in any range, the 
    // interval on the new range must be smaller by at least a sixth.
//NumAnalUtils.writeToNetLogo(Arrays.toString(lower)+Arrays.toString(upper), false, context);
//NumAnalUtils.writeToNetLogo(Arrays.toString(rng)+Arrays.toString(numIntervals), false, context);
//for (double[] entry : history) {
//  NumAnalUtils.writeToNetLogo(Arrays.toString(entry), false, context);
//}
    boolean done = true;
    for (int i = 0; i < n; i++) {
      double currentInterval = rng[i] / numIntervals[i];
//NumAnalUtils.writeToNetLogo("subdividing? "+currentInterval+" "+tol, false, context);
      if (currentInterval > tol) {
        done = false;
        lower[i] = currentMin.point[i] - currentInterval;
        upper[i] = currentMin.point[i] + currentInterval;
      }
    }
    if (!done) {
//NumAnalUtils.writeToNetLogo("subdividing "+Arrays.toString(lower)+Arrays.toString(upper), false, context);
      // At least one dimension has not reached its tolerance. Go back and 
      // create a new, more refined grid with the narrower boundaries.
      findGridMin(lower, upper, maxEvals, tol, currentMin, fnctn, context);
    }
  }

  private void findMin(int i, int n, double[][] ranges, double[] lower,
          CurrentMin currentMin, AnonymousReporter fnctn, Context context) {

    // calls itself recursively to work through the grid, evaluating the
    // function at each point and saving the current minimum in
    // currentMin.  Note that ties go to the first point with that value
    // encountered. Note too that the we start at the lower bound of each
    // range and work our way up though one demension at a time.
    double[] pt = lower.clone();
    for (int j = 0; j < ranges[i].length; j++) {
      pt[i] = ranges[i][j];
      if (i < (n - 1)) {
        findMin(i + 1, n, ranges, pt, currentMin, fnctn, context);
      } else {
        double val = NumAnalUtils.getFofXvec(pt, fnctn, context);
        if (val < currentMin.value) {
          currentMin.value = val;
          System.arraycopy(pt, 0, currentMin.point, 0, n);
        }
      }
    }
  }
}
