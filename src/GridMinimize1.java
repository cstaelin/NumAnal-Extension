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

public class GridMinimize1 implements Reporter {

    /* Ranges across a hyper-grid of points in n dimentions and evaluates
     the function at each point on the hyper-grid, returning the cordinates
     of the point yeilding the lowest value of the function.
     Takes four arguments:
     the function
     a list of length n containing the start of the range of values in
     each dimension.
     a list of length n containing the end of the range of values in each
     dimension.
     a list of length n containing the increment between points in each 
     direction.
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

    @Override
    public Syntax getSyntax() {
        return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
            Syntax.ListType(), Syntax.ListType(),
            Syntax.ListType()}, Syntax.NumberType(), 4);
    }

    @Override
    public Object report(Argument args[], Context context)
        throws ExtensionException, LogoException {

        // get the four arguments.
        AnonymousReporter fnctn = args[0].getReporter();
        LogoList startList = args[1].getList();
        LogoList endList = args[2].getList();
        LogoList incrementList = args[3].getList();

        double[] start = NumAnalUtils.convertSimpleLogoListToArray(startList);
        double[] end = NumAnalUtils.convertSimpleLogoListToArray(endList);
        double[] increment = NumAnalUtils.convertSimpleLogoListToArray(incrementList);

        int n = start.length;
        if (n != end.length || n != increment.length) {
            throw new ExtensionException("GridMinimize: start,end and increment "
                + "lists are not of the same length.");
        }

        // Construct the grid in a two-demensional, possibly ragged array.
        double ranges[][] = new double[n][];
        for (int i = 0; i < n; i++) {
            double s = start[i];
            double e = end[i];
            double incr = increment[i];
            // make sure start < end and incr is positive.
            if (s > e) {
                double tmp = s;
                s = e;
                e = tmp;
            }
            incr = Math.abs(incr);
            // fill in ranges[i][] with values from start[i] to end[i],
            // making sure that the last value is at end[i].  The last interval
            // may then be smaller than the others if the range divided by the 
            // interval has a fractional part.
            int m = (int) Math.ceil((e - s) / incr) + 1;
            ranges[i] = new double[m];
            for (int j = 0; j < (m - 1); j++, s += increment[i]) {
                ranges[i][j] = s;
            }
            ranges[i][m - 1] = e;
        }

        // Now call the findMin method which will work recusively through the grid.
        CurrentMin currentMin = new CurrentMin(n);
        findMin(0, n, ranges, start, currentMin, fnctn, context);

        return NumAnalUtils.convertArrayToSimpleLogoList(currentMin.point);
    }

    private void findMin(int i, int n, double[][] ranges, double[] pt,
        CurrentMin currentMin, AnonymousReporter fnctn, Context context) {

        // calls itself recursively to work through the grid, evaluating the
        // function at each point and saving the current minimum in
        // currentMin.  Note that ties go to the first point with that value
        // encountered.
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
