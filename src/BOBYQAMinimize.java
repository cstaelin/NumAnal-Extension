/*
 * This extension primitive makes use of the the BOBYQAOptimizer in the 
 * Apache Commons Math3 package. That software is licensed under the 
 * Apache Software Foundation license which grants in part 
 * "a perpetual, worldwide, non-exclusive, no-charge, royalty-free, 
 * irrevocable copyright license to reproduce, prepare Derivative Works of, 
 * publicly display, publicly perform, sublicense, and distribute the Work 
 * and such Derivative Works in Source or Object form."
 * 
 * The statement included with all the Apache Commons software follows:
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
// import org.nlogo.api.Command;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.core.LogoList;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.exception.TooManyEvaluationsException;


/*
 * BOBYQAMinimize finds the minimum of a multivariate function passed to it as 
 * an anonymous reporter, fnctn, and returns list of the input values at
 * the minimum.  fnctn should be a NetLogo ananymous reporter that takes
 * a list of input values and reports the value of the function at that
 * point. BOBYQAMinimize looks to the Bounds class to see if variable-by-variable
 * lower and/or upper bounds on the solution have been set by the user, 
 * and incorporates those bounds if they have.
 * 
 * This routine uses the Apache Commons library BOBYQAOptimizer class to 
 * actually perform the optimization.
 */
public class BOBYQAMinimize {

    static final int MAX_FUNCTION_EVALUATIONS = 50000;

    // Converts the NetLogo reporter into a MultivariateFunction type
    // that can be used by the Apache.commons.math3 methods. 
    public static class MVFunction implements MultivariateFunction {

        private final AnonymousReporter NLReporter;
        private final Context NLContext;

        public MVFunction(AnonymousReporter func, int n, Context context) {
            NLReporter = func;
            NLContext = context;
        }

        // MultivariateFunction requires that this method be defined to
        // return the value of the function at any point.
        @Override
        public double value(double[] xvec) {
            return NumAnalUtils.getFofXvec(xvec, NLReporter, NLContext);
        }
    }

    // Sets up the actual optimization.
    public static class BOBYQASolve implements Reporter {

        private int maxEvaluations = MAX_FUNCTION_EVALUATIONS;
        // The next two defaults are 10 and 1.0e-8, respectively.
        private double initial_radius = BOBYQAOptimizer.DEFAULT_INITIAL_RADIUS;
        private double stopping_radius = BOBYQAOptimizer.DEFAULT_STOPPING_RADIUS;

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

            // Get the arguments. The first two are required, the rest
            // optional. For the optional arguments, 0 represents the use
            // of the default value.
            //
            // Get the initial guess and turn it from a LogoList to a vector.
            // Check first to make sure that this is not a one-dimensional
            // problem. (The minimum dimension is 2.)
            LogoList xlist = args[0].getList();
            int nvar = xlist.size();
            if (nvar < BOBYQAOptimizer.MINIMUM_PROBLEM_DIMENSION) {
                throw new ExtensionException("BOBYQA requires at least a "
                        + BOBYQAOptimizer.MINIMUM_PROBLEM_DIMENSION
                        + "-dimensional problem.");
            }
            double[] x = new double[nvar];
            for (int j = 0; j < nvar; j++) {
                x[j] = ((Number) xlist.get(j)).doubleValue();
            }

            // We need to reset the optional parameters to their defaults as 
            // they may otherwise carryover from previous calls if they are not
            // specified in this call.
            initial_radius = BOBYQAOptimizer.DEFAULT_INITIAL_RADIUS;
            stopping_radius = BOBYQAOptimizer.DEFAULT_STOPPING_RADIUS;
            maxEvaluations = MAX_FUNCTION_EVALUATIONS;

            int nargs = args.length;
            if (nargs > 2) {
                double xx = args[2].getDoubleValue();
                initial_radius = (xx > 0.0) ? xx : BOBYQAOptimizer.DEFAULT_INITIAL_RADIUS;
            }
            if (nargs > 3) {
                double xx = args[3].getDoubleValue();
                stopping_radius = (xx > 0.0) ? xx : BOBYQAOptimizer.DEFAULT_STOPPING_RADIUS;
            }
            if (nargs > 4) {
                int mm = args[4].getIntValue();
                maxEvaluations = (mm > 0) ? mm : MAX_FUNCTION_EVALUATIONS;
            }
            // The number of interpolation points suggested by the literature
            // as a good default value is 2*nvar + 1.
            int numInterpolationPoints = 2 * nvar + 1;
            if (nargs > 5) {
                int mm = args[5].getIntValue();
                numInterpolationPoints = (mm > 0) ? mm : numInterpolationPoints;
            }

            // Check to see if the user has specified bounds on the solution.
            SimpleBounds smplBounds = null;
            boolean isBounded = false;
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
                    smplBounds = new SimpleBounds(Bounds.lowerBounds,
                            Bounds.upperBounds);
                    isBounded = true;
                } else {
                    // The bound dimentions do not match the problem.
                    throw new ExtensionException(
                            "The dimension of the problem, " + nvar
                            + ", does not match the dimension of the bounds, "
                            + Bounds.lowerBounds.length + ".");
                }
            }

            // Get the reporter task and turn it into a MultivariateFunction.
            AnonymousReporter fnctn = args[1].getReporter();
            MVFunction mvReporter = new MVFunction(fnctn, nvar, context);

            // Set up the Optimizer and get the resulting solution point.
            BOBYQAOptimizer optim = new BOBYQAOptimizer(
                    numInterpolationPoints, initial_radius, stopping_radius);

            // Find the solution.
            double rslt[];
            try {
                PointValuePair result = optim.optimize(
                        new MaxEval(maxEvaluations),
                        new ObjectiveFunction(mvReporter),
                        GoalType.MINIMIZE,
                        (isBounded) ? smplBounds : SimpleBounds.unbounded(nvar),
                        new InitialGuess(x));
                // put the solution in an array.
                rslt = result.getPoint();
            } catch (TooManyEvaluationsException e) {
                throw new ExtensionException(
                        "BOBYQA exceeded " + optim.getMaxEvaluations()
                        + " iterations, the maximum number specified.");
            }

            // Convert the result from an array to a NetLogo list and return it.
            return NumAnalUtils.convertArrayToSimpleLogoList(rslt);
        }
    }
}
