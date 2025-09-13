package org.nlogo.extensions.numanal;

import java.util.Arrays;
import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.core.LogoList;

import pal.math.DifferentialEvolution;
import pal.math.MultivariateFunction;
// import pal.math.MultivariateMinimum;

/*
 * Differential Evolution Search:
 * DESMinimize finds the minimum of a multivariate function passed to it in 
 * the anonymous reporter, fnctn, and returns list of the input values at
 * the minimum.  fnctn should be a NetLogo reporter that takes
 * a list of input values and reports the value of the function at that
 * point. The Bounds class is checked to see if the user has set lower and/or
 * upper bounds on the solution vector.
 * 
 * DESMinimize finds the minimum without using derivatives, employing a 
 * metaheuristic search.  It does not require gradients and
 * therefore is appropriate for poorly conditioned functions.
 * 
 * These routines use the pal.math library, v 1.4,
 * Copyright (c) 1999-2002 by the PAL Development Core Team
 * This package may be distributed under the terms of the
 * GNU Lesser General Public License (LGPL).
 * PAL Development Core Team:
 * Alexei Drummond, School of Biological Sciences, University of Auckland
 * Korbinian Strimmer, Department of Zoology, University of Oxford
 * Ed Buckler, Department of Genetics, North Carolina State University.
 *  http://iubio.bio.indiana.edu/soft/molbio/java/pal/
 */

public class DESMinimize
{

    // These are the default upper and lower bounds for the solution 
    // vector elements.
    static final double UPPER_BOUND_DEFAULT = 1.0E13;
    static final double LOWER_BOUND_DEFAULT = -UPPER_BOUND_DEFAULT;

    // Converts the NetLogo reporter into a MultivariateFunction type
    // that can be used by the pal.math methods. There are two constructors:
    // one that allows the solution to be unbounded, and one that allows
    // upper and lower bounds on the solution vector.  
    public static class MVFunction implements MultivariateFunction {

        private final AnonymousReporter NLReporter;
        private final int nvars;
        private final Context NLContext;
        private final double[] lowerBounds;
        private final double[] upperBounds;

        public MVFunction(AnonymousReporter func, int n, Context context) {
            NLReporter = func;
            nvars = n;
            NLContext = context;
            lowerBounds = new double[nvars];
            upperBounds = new double[nvars];
            Arrays.fill(lowerBounds, LOWER_BOUND_DEFAULT);
            Arrays.fill(upperBounds, UPPER_BOUND_DEFAULT);
        }
        
        public MVFunction(AnonymousReporter func, int n, double[] lowerBnds, 
                double[] upperBnds, Context context) {
            NLReporter = func;
            nvars = n;
            NLContext = context;
            lowerBounds = new double[nvars];
            upperBounds = new double[nvars];
            // In filling the upperBounds and lowerBounds arrays, check for
            // "real" values as the Bounds class uses + and - infinity as 
            // its unbounded values.
            for (int i = 0; i < nvars; i++) {
                lowerBounds[i] = (lowerBnds[i] != Bounds.LOWER_BOUND_DEFAULT) ?
                    lowerBnds[i] : LOWER_BOUND_DEFAULT;
                upperBounds[i] = (upperBnds[i] != Bounds.UPPER_BOUND_DEFAULT) ?
                    upperBnds[i] : UPPER_BOUND_DEFAULT;
            }
        }

        // MultivariateFunction requires that these methods be defined.
        @Override
        public double evaluate(double[] xvec) {
            return NumAnalUtils.getFofXvec(xvec, NLReporter, NLContext);
        }

        @Override
        public int getNumArguments() {
            return nvars;
        }

        @Override
        public double getLowerBound(int n) {
            return lowerBounds[n];
        }

        @Override
        public double getUpperBound(int n) {
            return upperBounds[n];
        }
    }

    public static class DESSolve implements Reporter {

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

            // Get the initial guess and turn it from a LogoList to a vector.
            LogoList xlist = args[0].getList();
            int nvar = xlist.size();
            double[] x = NumAnalUtils.convertSimpleLogoListToArray(xlist);

            if (Bounds.lowerBounds != null) {
                // There are bounds. Check that they are of the right
                // dimension and that the initial point is within them.
                if (Bounds.lowerBounds.length == nvar) {
                    for (int i = 0; i < nvar; i++) {
                        if (x[i] < Bounds.lowerBounds[i] || x[i] > Bounds.upperBounds[i]) {
                            throw new ExtensionException("The initial guess "
                                    + "must be within the specified bounds.");
                        }
                    }
                } else {
                    // The bound dimentions do not match the problem.
                    throw new ExtensionException(
                            "The dimension of the problem, " + nvar
                            + ", does not match the dimension of the bounds, "
                            + Bounds.lowerBounds.length + ".");
                }
            }

            // Get the reporter task.
            AnonymousReporter fnctn = args[1].getReporter();

            // Check for the optional arguments, first setting defaults to
            // zero to indicate use of the pal.math library defaults.
            double tolfx = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            double tolx =  Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            int popSize = 0;        // zero indicates use of the default.
            double CR = 0.0;        // zero indicates use of the default = 0.9.
            double F = 0.0;         // zero indicates use of the default = 0.7.
            int nFun = 0;           // the default is no limit.

            int nargs = args.length;
            if (nargs > 2) {
                tolfx = args[2].getDoubleValue();
            }
            if (nargs > 3) {
                tolx = args[3].getDoubleValue();
            }
            if (nargs > 4) {
                popSize = args[4].getIntValue();
            }
            if (nargs > 5) {
                nFun = args[5].getIntValue();
            }
            if (nargs > 6) {
                CR = args[6].getDoubleValue();
            }
            if (nargs > 7) {
                F = args[7].getDoubleValue();
            }

            // Construct the a new MVFunction from the NetLogo reporter.
            MVFunction mvReporter;
            if (Bounds.lowerBounds == null) {
                mvReporter =
                        new MVFunction(fnctn, nvar, context);
            } else {
                mvReporter = new MVFunction(fnctn, nvar, Bounds.lowerBounds, 
                        Bounds.upperBounds, context);
            }

            // Then construct a new DifferentialEvolutionn and use its 
            // optimize method to find the solution vector that minimizes 
            // the value of the reporter. NOTE: x goes into the method as
            // the initial guess and comes out as the solution vector.
            // Note that we use the optimize method rather than findMinimum
            // because of the strange way that the latter accepts the values
            // of the tolerances.
            DifferentialEvolution desearch;
            if (popSize <= 0) {
                // Let the method calculate an appropriate population size.
                desearch = new DifferentialEvolution(nvar);
            } else {
                desearch = new DifferentialEvolution(nvar, popSize);
            }

            if (nFun > 0) {
                desearch.maxFun = nFun;
            }
            if (CR > 0.0) {
                desearch.CR = CR;
            }
            if (F > 0.0) {
                desearch.F = F;
            }

            /* tolx and tolfx are the stopping tolerances for the elements of
             * the solution vector and the value of the objective function,
             * repectively.  One or the other may be zero, but not both.  If 
             * both are zero, Bounds.ABSOLUTE_TOLERANCE_DEFAULT is used for both.
             * tolx - the optimization stops if on two successive evaluations
             * every element of the solution vector is within tolx of its
             * prior value.
             * tolfx - the optimization stops if over a specified number of 
             * evaluations of the objective function, the value of the 
             * objective function changes by less than tolfx at each 
             * evaluation. This presumably means that the value could change 
             * in the same directionby an amount less than tolfx each time,
             * meaning that the lowest and highest value differ by n*tolfx,
             * where n is the specified number of successive evaluations.
             * The default n is given by DifferentialEvolution.numFuncStops,
             * which is set at 4. It may be changed by changing the value
             * of that field.
            */

            if (tolfx <= 0.0 && tolx <= 0.0) {
                tolfx = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
                tolx = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            }
            
//            desearch.numFuncStops = 4;  // uncomment if you want to change this.
            desearch.optimize(mvReporter, x, tolfx, tolx);

            return NumAnalUtils.convertArrayToSimpleLogoList(x);
        }
    }
}