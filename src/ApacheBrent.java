/*
 * This extension primitive makes use of the the BrentOptimizer in the 
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

// Could add an optional initial guess and max evaluations.
import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.MaxEval;

/*
 * ApacheBrent employs the Brent alorithm to minimize a function 
 * of one variable, x. The return value is the value of x between an upper
 * and lower bound where the function takes on its minimum value. 
 * The arguments are:
 * fnctn - the function to be minimized, passed as a NetLogo task variable.
 * lowBound and highBound - the bounds between which the minimum is to be 
 * found. Note that if the range between the two bounds does not contain
 * the "true" minimum of the function, the algorithm returns the bound 
 * closest to the minimim.
 * 
 * The relTolerance and absTolerance arguments are optional.
 * Apache's Brent Optimizer uses a solution tolerance of 
 *        tol = relTolerance*abs(x) + absTolerance.
 * relTolerance must be > 2 times the machine precision and absTolerance 
 * must be > 0.  If either is not specified or passed as zero, the default 
 * values Bounds.RELATIVE_TOLERANCE_DEFAULT and Bounds.ABSOLUTE_TOLERANCE_DEFAULT
 * are used.
 */
public class ApacheBrent {

  // Default for the maximum number of iterations
  static final int MAX_EVALUATIONS_DEFAULT = 1000;

    // Converts the NetLogo reporter into a UnivariateFunction type
  // that can be used by the Apache.commons.math3 methods.
  public static class UVFunction implements UnivariateFunction {

    private final AnonymousReporter NLReporter;
    private final Context NLContext;

    public UVFunction(AnonymousReporter func, Context context) {
      NLReporter = func;
      NLContext = context;
    }

        // UnivariateFunction requires that this method be defined to return 
    // the value of the function at any x.
    @Override
    public double value(double x) {
      return NumAnalUtils.getFofX(x, NLReporter, NLContext);
    }
  }

  // Performs the actual solution.
  public static class BrentSolve implements Reporter {

    private double relTolerance = Bounds.RELATIVE_TOLERANCE_DEFAULT;
    private double absTolerance = Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
    private final int maxEvaluations;

    public BrentSolve() {
      // This could be made an argument to the primitive.
      this.maxEvaluations = MAX_EVALUATIONS_DEFAULT;
    }

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
        Syntax.NumberType(), Syntax.NumberType(),
        Syntax.NumberType() | Syntax.RepeatableType()}, Syntax.NumberType(), 3);
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      // Takes at least three arguments: the task function itself and
      // the lower and upper bounds of the solution interval.
      // Optionally takes the absolute and relative tolerances with 
      // a negative number for either signalling the default value.
      // One, but not both, may be zero.
      int nargs = args.length;
      AnonymousReporter fnctn = args[0].getReporter();
      double lowBound = args[1].getDoubleValue();
      double highBound = args[2].getDoubleValue();
      if (nargs > 3) {
        double xx = args[3].getDoubleValue();
        relTolerance = (xx > 0.0) ? xx : Bounds.RELATIVE_TOLERANCE_DEFAULT;
      }
      if (nargs > 4) {
        double xx = args[4].getDoubleValue();
        absTolerance = (xx > 0.0) ? xx : Bounds.ABSOLUTE_TOLERANCE_DEFAULT;
      }

      // Check for "real" bounds in case Bounds.UPPER_BOUND_DEFAULT
      // or Bounds.LOWER_BOUND_DEFAULT were passed.
      if (highBound == Double.POSITIVE_INFINITY) {
        highBound = Double.MAX_VALUE;
      }
      if (lowBound == Double.NEGATIVE_INFINITY) {
        lowBound = -Double.MAX_VALUE;
      }

      // Make sure that the lower and upper bounds are in ascending order.
      if (lowBound > highBound) {
        double temp = lowBound;
        lowBound = highBound;
        highBound = temp;
      }

      // Turn the task function into a UnivariateFunction.
      UVFunction uvReporter = new UVFunction(fnctn, context);

      // Arguments are the relative and absolute tolerances.
      UnivariateOptimizer optimizer = new BrentOptimizer(relTolerance,
              absTolerance);

      // find the minimum
      final UnivariatePointValuePair result = optimizer.optimize(
              new MaxEval(maxEvaluations),
              new UnivariateObjectiveFunction(uvReporter),
              GoalType.MINIMIZE,
              new SearchInterval(lowBound, highBound));

      // Return the solution.
      return result.getPoint();
    }
  }
}
