/*
 * This extension primitive makes use of the the CMAESOptimizer in the 
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
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.exception.TooManyEvaluationsException;


/*
 * CMAESMinimize finds the minimum of a multivariate function passed to it as 
 * an anonymous reporter, fnctn, and returns list of the input values at
 * the minimum.  fnctn should be a NetLogo anonymous reporter that takes
 * a list of input values and reports the value of the function at that
 * point. CMAESMinimize looks to the Bounds class to see if variable-by-variable
 * lower and/or upper bounds on the solution have been set by the user, 
 * and incorporates those bounds if they have.
 * 
 * This routine uses the Apache Commons library CMAESOptimizer class to 
 * actually perform the optimization.
 */
public class CMAESMinimize {

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
    public static class CMAESSolve implements Reporter {

        private int maxEvaluations = MAX_FUNCTION_EVALUATIONS;
        private int maxIterations = MAX_FUNCTION_EVALUATIONS;
        private double relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
        private double absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
        // stopFitness: stop if the range of objective function values
        // over the past and current generations is less than this.
        private static final double stopFitness = 1.0e-12;
        // suggested defaults for isActiveCMA and diagonalOnly.
        // isActiveCMA = true turns on "active CMA" with a negative update 
        // of the covariance matrix and checks for positive definiteness.
        // diagonalOnly defines the number of initial iterations where the 
        // covariance matrix remains diagonal and the algorithm has internally
        // linear time complexity. diagonalOnly = 1 means keeping the 
        // covariance matrix always diagonal and this setting also exhibits 
        // linear space complexity. This can be particularly useful for 
        // dimension > 100. Default = 0, always use full covariance matrix.
        private boolean isActiveCMA = true;
        private int diagonalOnly = 0;
        
        // Determines how often new random objective variables are generated
        // in case they are out of bounds or infeasible. Default is 0.
        private int checkFeasibleCount = 0;
        private final MersenneTwister random = new MersenneTwister();
        boolean generateStatistics = false;

        @Override
        public Syntax getSyntax() {
            return SyntaxJ.reporterSyntax(new int[]{Syntax.ListType(),
                        Syntax.ReporterType(), Syntax.ListType(),
                        Syntax.WildcardType() | Syntax.RepeatableType()},
                    Syntax.ListType(), 3);
        }

        @Override
        public Object report(Argument args[], Context context)
                throws ExtensionException, LogoException {

            // Get the initial guess and turn it from a LogoList to a vector.
            LogoList xlist = args[0].getList();
            int nvar = xlist.size();
            double[] x = NumAnalUtils.convertSimpleLogoListToArray(xlist);
            
            // Get the reporter task and turn it into a MultivariateFunction.
            AnonymousReporter fnctn = args[1].getReporter();
            MVFunction mvReporter = new MVFunction(fnctn, nvar, context);
            
            // Get the list of sigma valuess. Sources differ on the best
            // values, but suggest values ranging from one-third to one times 
            // the distance from the guess for a given variable to its likely
            // value at the minimum of the function.
            LogoList slist = args[2].getList();
            double[] inSigma = NumAnalUtils.convertSimpleLogoListToArray(slist);
            if (inSigma.length != nvar) {
                throw new ExtensionException("The length of the list of sigma values "
                        + "must be be the same as the length of the guess.");
            }

            // The default population size, lambda, is given by this commonly
            // accepted expression. Larger values increase the likelihood
            // of avoiding a local minimum, but slow down execution.
            int lambda = 4 + (int) (3.0 * Math.log(nvar));

            // Get the optional arguments
            int nargs = args.length;
            if (nargs > 3) {
                int mm = args[3].getIntValue();
                lambda = (mm > 0) ? mm : lambda;
            }
            if (nargs > 4) {
                double xx = args[4].getDoubleValue();
                relTolerance = (xx > 0.0) ? xx : -1.0;
                // negative value means use absTolerance only.
            }
            if (nargs > 5) {
                double xx = args[5].getDoubleValue();
                absTolerance = (xx > 0.0) ? xx : -1.0;
                // negative value means use relTolerance only.
            }
            if (nargs > 6) {
                int mm = args[6].getIntValue();
                maxEvaluations = (mm > 0) ? mm : MAX_FUNCTION_EVALUATIONS;
                maxIterations = 10 * maxEvaluations;
            }
            if (nargs > 7) {
                int mm = args[7].getIntValue();
                checkFeasibleCount = (mm > 0) ? mm : checkFeasibleCount;
            }
            if (nargs > 8) {
                isActiveCMA = args[8].getBooleanValue();                
            }
            if (nargs > 9) {
                int mm = args[9].getIntValue();
                diagonalOnly = (mm >= 0) ? mm : diagonalOnly;
            }

            // Need to make sure both tols are not negative. (One is okay.)
            // If they are, assume user wanted defaults for both.
            if (relTolerance < 0.0 && absTolerance < 0.0) {
                relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
                absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
            }

            // Set the stopping criteria for the objective function.
            SimpleValueChecker checker = new SimpleValueChecker(relTolerance,
                    absTolerance);

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
                            "The dimension of the problem " + nvar
                            + " does not match the dimension of the bounds "
                            + Bounds.lowerBounds.length + ".");
                }
            }

            CMAESOptimizer optim = new CMAESOptimizer(
                    maxIterations, stopFitness, isActiveCMA, diagonalOnly,
                    checkFeasibleCount, random, generateStatistics, checker);

            double[] rslt;
            try {
                PointValuePair result = optim.optimize(
                        new MaxEval(maxEvaluations),
                        new ObjectiveFunction(mvReporter),
                        GoalType.MINIMIZE,
                        (isBounded) ? smplBounds : SimpleBounds.unbounded(nvar),
                        new CMAESOptimizer.Sigma(inSigma),
                        new CMAESOptimizer.PopulationSize(lambda),
                        new InitialGuess(x));

                // Get the solution.
                rslt = result.getPoint();
            } catch (TooManyEvaluationsException e) {
                throw new ExtensionException(
                        "CMAES-minimize exceeded " + optim.getMaxEvaluations()
                        + " iterations, the maximum number specified.");
            }

            // Convert the solution to a LogoList and return it.
            return NumAnalUtils.convertArrayToSimpleLogoList(rslt);
        }
    }
}
