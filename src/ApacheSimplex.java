/*
 * This extension primitive makes use of the the SimplexOptimizer in the 
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
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
// import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionPenaltyAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.exception.TooManyEvaluationsException;


/*
 * ApacheSimplex finds the minimum of a multivariate function passed to it as  
 * an anonymous reporter, fnctn, and returns list of the input values at
 * the minimum.  fnctn should be a NetLogo reporter that takes
 * a list of input values and reports the value of the function at that
 * point.  ApacheSimplex looks to the Bounds class to see if variable-by-variable
 * lower and/or upper bounds on the solution have been set by the user, 
 * and incorporates those bounds if they have.
 * 
 * This routine uses the Apache Commons library SimplexOptimizer class to 
 * actually perform the optimization.
 */
public class ApacheSimplex {

    // These specify wich simplex method is to be used.
    static final int MULTI_DIRECTIONAL_SIMPLEX = 0;
    static final int NELDER_MEAD_SIMPLEX = 1;
    // These are the defaults listed in the Apache API.
    static final double KHI_DEFAULT = 2.0;    // Expansion coefficient
    static final double GAMMA_DEFAULT = 0.5;  // Contraction coefficient
    static final double RHO_DEFAULT = 1.0;    // Refection coefficient
    static final double SIGMA_DEFAULT = 0.5;  // Shrinkage coefficient
    static final int MAX_FUNCTION_EVALUATIONS = 5000;
    static final double SIDE_LENGTH_DEFAULT = 10.0;

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

    // Sets up the actual optimization
    public static class SimplexSolve implements Reporter {

        // These are the various default values.
        private int solutionType = NELDER_MEAD_SIMPLEX;
        private int maxEvaluations = MAX_FUNCTION_EVALUATIONS;
        private double sideLength = SIDE_LENGTH_DEFAULT;
        private double relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
        private double absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;

        // The constructor determines which Simplex method is used.
        public SimplexSolve(int slnType) {
            solutionType = slnType;
        }

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
            // of the default value, except for the tolerances, one of which
            // can be zero. In those cases < 0 signifies the default.
            // Get the initial guess and turn it from a LogoList to a vector.
            LogoList xlist = args[0].getList();
            int nvar = xlist.size();
            double[] x = NumAnalUtils.convertSimpleLogoListToArray(xlist);

            // We need to reset the optional parameters to their defaults as 
            // they may otherwise carryover from previous calls if they are not
            // specified in this call.
            sideLength = SIDE_LENGTH_DEFAULT;
            relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
            absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            maxEvaluations = MAX_FUNCTION_EVALUATIONS;
            
            int nargs = args.length;
            if (nargs > 2) {
                double xx = args[2].getDoubleValue();
                sideLength = (xx > 0.0) ? xx : SIDE_LENGTH_DEFAULT;
            }
            if (nargs > 3) {
                double xx = args[3].getDoubleValue();
                relTolerance = (xx > 0.0) ? xx : -1.0;
                // negative value means use absTolerance only.
            }
            if (nargs > 4) {
                double xx = args[4].getDoubleValue();
                absTolerance = (xx > 0.0) ? xx : -1.0;
                // negative value means use relTolerance only.
            }
            if (nargs > 5) {
                int mm = args[5].getIntValue();
                maxEvaluations = (mm > 0) ? mm : MAX_FUNCTION_EVALUATIONS;
            }

            // Need to make sure both tols are not negative. (One is okay.)
            // If they are, assume user wanted defaults for both.
            if (relTolerance < 0.0 && absTolerance < 0.0) {
                relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
                absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            }

            // Get the reporter task and turn it into a MultivariateFunction.
            AnonymousReporter fnctn = args[1].getReporter();
            MVFunction mvReporter = new MVFunction(fnctn, nvar, context);

            // The SimplexOptimizer does not directly support boundaries on 
            // the solution. Therefore, if there are, we wrap the function in 
            // a MultivariateFunctionMappingAdapter. This returns function 
            // values that indicate to the Optimizer that it is approaching 
            // a boundary and that it is not optimal to push this variable
            // any further.

            // Check to see if the user has specified bounds on the solution.
            boolean isBounded = false;
            if (Bounds.lowerBounds != null) {
                // There are bounds. Check that they are of the right
                // dimension and that the initial point is strictly within them.
                if (Bounds.lowerBounds.length == nvar) {
                    for (int i = 0; i < nvar; i++) {
                        if (x[i] <= Bounds.lowerBounds[i] || x[i] >= Bounds.upperBounds[i]) {
                            throw new ExtensionException("The initial guess "
                                    + "must be strictly within the specified bounds.");
                        }
                    }
                    isBounded = true;
                } else {
                    // The bound dimentions do not match the problem.
                    throw new ExtensionException(
                            "The dimension of the problem " + nvar
                            + " does not match the dimension of the bounds "
                            + Bounds.lowerBounds.length + ".");
                }
            }

            // Constructor arguments are the relative and absolute threshholds.
            SimplexOptimizer optimizer = new SimplexOptimizer(relTolerance, absTolerance);

            // These are default values.  We could have a method to allow
            // them to be set by the user, or add them to the list of calling
            // arguments.
            double khi = KHI_DEFAULT;   // Expansion coefficient
            double gamma = GAMMA_DEFAULT; // Contraction coefficient
            double rho = RHO_DEFAULT;   // Reflection coefficient
            double sigma = SIGMA_DEFAULT; // Shrinkage coefficient

            // Perform the appropriate optimization.
            double[] solution;
            if (!isBounded) {
                try {
                    PointValuePair result = optimizer.optimize(
                            new MaxEval(maxEvaluations),
                            new ObjectiveFunction(mvReporter),
                            GoalType.MINIMIZE,
                            (solutionType == MULTI_DIRECTIONAL_SIMPLEX)
                            ? new MultiDirectionalSimplex(nvar, sideLength, khi, gamma)
                            : new NelderMeadSimplex(nvar, rho, khi, gamma, sigma),
                            new InitialGuess(x));
                    // Put the solution in an array.
                    solution = result.getPoint();
                } catch (TooManyEvaluationsException e) {
                    throw new ExtensionException(
                            "Simplex exceeded " + optimizer.getMaxEvaluations()
                            + " iterations, the maximum number specified.");
                }
            } else {
                // Create the bounds wrapper.
                MultivariateFunctionMappingAdapter wrappedMVReporter =
                        new MultivariateFunctionMappingAdapter(mvReporter,
                        Bounds.lowerBounds, Bounds.upperBounds);
                PointValuePair result;
                // Find the solution.
                try {
                    result = optimizer.optimize(
                            new MaxEval(maxEvaluations),
                            new ObjectiveFunction(wrappedMVReporter),
                            GoalType.MINIMIZE,
                            (solutionType == MULTI_DIRECTIONAL_SIMPLEX)
                            ? new MultiDirectionalSimplex(nvar, sideLength, khi, gamma)
                            : new NelderMeadSimplex(nvar, rho, khi, gamma, sigma),
                            new InitialGuess(wrappedMVReporter.boundedToUnbounded(x)));
                } catch (TooManyEvaluationsException e) {
                    throw new ExtensionException(
                            "Simplex exceeded " + optimizer.getMaxEvaluations()
                            + " iterations, the maximum number specified.");
                }
                // Unwrap the solution and put it in an array.
                solution = wrappedMVReporter.unboundedToBounded(result.getPoint());
            }

            // Convert the solution to a LogoList and return it.
            return NumAnalUtils.convertArrayToSimpleLogoList(solution);
        }
    }
}