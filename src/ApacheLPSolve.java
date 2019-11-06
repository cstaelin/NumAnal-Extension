/*
 * ApacheLPSolve uses the Apache math3 library to solve relatively simple
 * linear programing pooblems. It will handle mixed constraints and 
 * unconstrained variables, but not integer problems. It provides solutions
 * for both the primal and the dual of the problem and so does yield
 * shadow prices. The Apache library does not provide a solution for the dual,
 * so we construct the dual explicitly, and perhaps less efficiently than we
 * might, but in a way that makes the method clear. I have tried to make 
 * setting up the problem in NetLogo as simple as possible. See the 
 * documentation for details.
 *
 * Now, some boilerplate.
 * This extension primitive makes use of the the
 * Apache Commons Math3 package to solve simple linear programming problems.
 * That software is licensed under the 
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

import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.core.LogoList;
import org.nlogo.api.LogoListBuilder;

import org.apache.commons.math3.optim.linear.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.exception.TooManyIterationsException;

import java.util.ArrayList;
import java.util.Arrays;


/*
 * ApacheLPSove uses the Apache math3 library to solve simple linear 
 * programming problems and their associated duals.  It can handle mixed 
 * constraints and unconstrained primal variables.
 */
public class ApacheLPSolve {

    static final int SOLVE_PRIMAL = 0;
    static final int SOLVE_DUAL = 1;

    // Parses each constraint, passed as a LogoList.
    private static LinearConstraint ParseConstraint(LogoList constraintList,
            int nvars, GoalType goal, Context context) throws ExtensionException {
        // the first item in the list is a list of the constraint coefficients.
        Object item = constraintList.get(0);
        double[] coeffs = new double[nvars];
        if (item instanceof LogoList) {
            coeffs = NumAnalUtils.convertSimpleLogoListToArray((LogoList) item);
            if (coeffs.length != nvars) {
                throw new ExtensionException("This constraint has different number of  "
                        + "coefficients, " + coeffs.length + ", than there are "
                        + "variables, " + nvars + ".\n" + constraintList.toString());
            }
        } else {
            throw new ExtensionException("the first item in the constraint "
                    + "is not a list of coefficients." + constraintList.toString());
        }
        // the second item in the list is the constraint relationship or the 
        // RHS value of the constraint. If the latter, the relationship is 
        // assumed, based on the goal, <= for maximizations and >= for 
        // minimizations.
        item = constraintList.get(1);
        Relationship relationship = (goal == GoalType.MAXIMIZE)
                ? Relationship.LEQ : Relationship.GEQ;
        boolean haveValue = false;
        double constant = 0.0;
        if (item instanceof Double) {
            // we've found the RHS. Save it and note that we don't need 
            // to look for a third item in the constraint list.
            constant = (double) item;
            haveValue = true;
        } else if (item instanceof String) {
            // we've found the relationship.
            String reln = ((String) item).toLowerCase();
            if (reln.contains("ge") || reln.contains(">=")) {
                relationship = Relationship.GEQ;
            } else if (reln.contains("le") || reln.contains("<=")) {
                relationship = Relationship.LEQ;
            } else if (reln.contains("eq") || reln.contains("=")) {
                relationship = Relationship.EQ;
            } else {
                throw new ExtensionException("Unknown relation in constraint:"
                        + constraintList.toString());
            }
        }
        // if the second item was not the constraint value, look for it here.
        if (!haveValue) {
            item = constraintList.get(2);
            if (item instanceof Double) {
                constant = (double) item;
            } else {
                throw new ExtensionException("Illegal RHS in the constraint: "
                        + constraintList.toString());
            }
        }

        return new LinearConstraint(coeffs, relationship, constant);
    }
    
    private static int RemoveNonNegConstraints(
            ArrayList<LinearConstraint> constraints, boolean[] isFree,
            Context context) {
        // run through the constraints to find those with >= 0, and only 
        // one non-zero variable coefficient. Work from the back so that
        // the for loop will behave properly even if a constraint is 
        // removed. Note variables whose nonNeg constraint has been removed
        // in the isFree array. Return the new number of constraints.
        int ncon = constraints.size();
        int nvar = constraints.get(0).getCoefficients().toArray().length;
        for (int i = ncon - 1; i >= 0; i--) {
            if (constraints.get(i).getValue() == 0.0
                    && constraints.get(i).getRelationship() == Relationship.GEQ) {
                // we've found a candidate. check the coefficients.
                double[] coeffs = constraints.get(i).getCoefficients().toArray();
                // check the number and placement of zero coefficients in
                // the constraint.
                int numZeros = 0;
                boolean isZero[] = new boolean[nvar];
                Arrays.fill(isZero, false);
                for (int j = 0; j < nvar; j++) {
                    if (coeffs[j] == 0.0) {
                        numZeros++;
                        isZero[j] = true;
                    }
                }
                if (numZeros == nvar - 1) {
                    // we've found a nonNegative constraint on a variable
                    // as there is only one non-zero coefficient.
                    // note variable, i.e., the one with the non-zero 
                    // coefficient, by putting a true in the removed array.
                    // Then remove the constraint.
                    for (int j = 0; j < nvar; j++) {
                        if (!isZero[j]) {
                            isFree[j] = false;
                        }
                    }
                    constraints.remove(i);
                    // NumAnalUtils.printString("Non-negative constraint " + i + " has been removed.", context);
                }
            }
        }
        return constraints.size();
    }

    // Sets up the actual linear optimization
    public static class LPSimplexSolve implements Reporter {

        int solutionType;

        // This constructor determines whether the primal or the dual is to 
        // be solved.
        public LPSimplexSolve(int primal_dual) {
            solutionType = primal_dual;
        }

        // Required arguments: ojective function, constraints.
        // Optional arguments: maximize?, bounded-by-zero?
        @Override
        public Syntax getSyntax() {
            return SyntaxJ.reporterSyntax(new int[]{
                Syntax.ListType(),
                Syntax.ListType(),
                Syntax.StringType() | Syntax.BooleanType() | 
                        Syntax.ListType() | Syntax.RepeatableType()},
                    Syntax.ListType(), 2);
        }

        @Override
        public Object report(Argument args[], Context context)
                throws ExtensionException {

            // Get the arguments. The first two are required, the rest are 
            // optional.
            int nargs = args.length;

            // 1. the objective function.
            LogoList objectiveList = args[0].getList();
            int nvar = objectiveList.size();
            double[] objective = NumAnalUtils.convertSimpleLogoListToArray(objectiveList);

            // 3. If there is a third argument, it is the goal. If not, the
            // goal defaults to a maximization.  We need the goal now in order 
            // to form the constraints with the proper defaults.
            // The goal argument must be either a boolean, true => maximize
            // and false => minimize, or a string. In the absence of this 
            // argument, a maximization is assumed.
            GoalType goal = GoalType.MAXIMIZE;
            if (nargs > 2) {
                if (args[2].get() instanceof Boolean) {
                    // it's a boolean
                    goal = (args[2].getBoolean()) ? GoalType.MAXIMIZE : GoalType.MINIMIZE;
                } else if (args[2].get() instanceof String) {
                    // it's a string.
                    String minMax = args[2].getString().toLowerCase();
                    if (minMax.contains("max")) {
                        goal = GoalType.MAXIMIZE;
                    } else if (minMax.contains("min")) {
                        goal = GoalType.MINIMIZE;
                    } else {
                        throw new ExtensionException(minMax
                                + " specifies neither 'maximize' nor 'minimize'.");
                    }
                } else {
                    // it must be a list, given the syntax specification.
                    // For the moment, throw an error.
                    throw new ExtensionException("The third arguemtn must be " +
                           "either true/false or a string indicating " +
                            "maximize (true/\"max\") or minimize (false/\"min\").");
                }
            }

            // 2. the second argumentcontains the constraints as a list of lists.
            LogoList nestedLogoList = args[1].getList();
            int ncon = nestedLogoList.size();
            // iterate through the list of lists to extract and check each constraint.
            ArrayList<LinearConstraint> constraints = new ArrayList<>();
            for (Object obj : nestedLogoList.toJava()) {
                if (obj instanceof LogoList) {
                    constraints.add(ParseConstraint((LogoList) obj, nvar,
                            goal, context));
                } else {
                    throw new ExtensionException("The constraint list "
                            + "is not formed correctly. " + "\n"
                            + "It should be a list of lists, with each sublist "
                            + "being one of the contraints.");
                }
            }

            // 4. the fourth, optional, parameter is a boolean or 
            // a string indicating whether ALL the variables are constrained
            // to be non-negative, or a list indicating which variables are
            // constrained and which are not.  The default is that they ALL
            // are constrained to be non-negative.
            // Set the optional parameter to its default.
            boolean nonNegative = true;
            if (nargs > 3) {
                if (args[3].get() instanceof Boolean) {
                    // the argument is a boolean, true for all variables
                    // constrained to non-negative, or false, indicating that
                    // one or more variables is "free". Any variables that are 
                    // constrained will presumably have their non-negative 
                    // constraint included in the constraint set.
                    nonNegative = args[3].getBoolean();
                } else if (args[3].get() instanceof String) {
                    // the argument is a string.
                    String nonNeg = args[3].getString().toLowerCase();
                    if (nonNeg.contains("nonneg")) {
                        nonNegative = true;
                    } else if (nonNeg.contains("free") || nonNeg.equals("")) {
                        nonNegative = false;
                    } else {
                        throw new ExtensionException(nonNeg
                                + " specifies neither 'nonNegative' nor 'free'.");
                    }
                }
                else {
                    // We must have a list of booleans or variable numbers.
                    // If booleans, there will be one for each primal variable,
                    // true indicating that the variable is constrained to be
                    // non-negative and false indicating that the variable is 
                    // not constrained, i.e., "free". If variable numbers, it 
                    // will contain the variables whihc should be set to non-
                    // negative.
                    // Get and check the list.
                    nonNegative = false;
                    LogoList isConstrainedList = args[3].getList();
                    boolean isConstrained[] = new boolean[nvar];
                    Arrays.fill(isConstrained, false);
                    // check that the list is not empty. If it is, we stick 
                    // with the default that all variables are free since
                    // there is a 4th argument and nonNegative was not spedified.
                    if (!isConstrainedList.isEmpty()) {
                        if (isConstrainedList.get(0) instanceof Boolean) {
                            // At least the first element of the list is a boolean.
                            // Assume (perhaps incorrectly?) that they all are.
                            if (isConstrainedList.length() != nvar) {
                                throw new ExtensionException("The true/false list "
                                        + "specifying non-negative variables must "
                                        + "have length " + nvar);
                            }
                            isConstrained = NumAnalUtils.convertSimpleLogoListToBooleanArray(isConstrainedList);                      
                        } else if (isConstrainedList.get(0) instanceof Double) {
                            // The first element of the list is a number.
                            // Assume then that the list contains the numbers of 
                            // the constrained-to-non-negative variables.
                            double x[] = NumAnalUtils.convertSimpleLogoListToArray(isConstrainedList);
                            for (double k : x) {
                                int j = (int)Math.round(k);
                                if (j >= 0 && j < nvar) {
                                    isConstrained[j] = true;
                                } else {
                                    throw new ExtensionException("The fouth " +
                                            "argument contains a variable reference " +
                                            "ourside the range of variables, " +
                                            "zero to " + nvar);
                                }
                            }
                        } else {
                            // It's an illegal fourth argument.
                            throw new ExtensionException("fourth argument contains "
                                    + "a list with neither true/false nor "
                                    + "variable numbers: " + 
                                    isConstrainedList.toString());
                        }
                    }
                    // NumAnalUtils.printString("non-negative list " + Arrays.toString(isConstrained), context);
                    
                    // We assume that the isConstrained list dominates any
                    // explicit non-negative constraints in the constraint set.
                    // Thus, we somewhat inefficiently delete any explicit 
                    // constraints that might be there, and then create new
                    // ones accrding to the isConstrained list. In this case,
                    // the isFree argument is superfluous. We don't use it.
                    boolean isFree[] = new boolean[nvar];
                    ncon = RemoveNonNegConstraints(constraints, isFree, context);
                    // check to make sure that not all the variables are 
                    // constrained.
                    boolean b = true;
                    for (boolean k : isConstrained) {
                        b = b & k;
                    }
                    if (b) {
                        // all the variables were constrained!  Simply set the 
                        // non-negative flag to true.
                        nonNegative = true;
                    } else {
                        // we do need to go through the list and add the 
                        // appropriate constraints.
                        for (int j = 0; j < nvar; j++) {
                            if (isConstrained[j]) {
                                // add a constraint for this variable.
                                double coeffs[] = new double[nvar];
                                Arrays.fill(coeffs, 0.0);
                                coeffs[j] = 1.0;
                                LinearConstraint constr = new LinearConstraint(coeffs, Relationship.GEQ, 0.0);
                                constraints.add(constr);
                                ncon++;
                                // NumAnalUtils.printString("Non-negative constraint " + j + " has been added.", context);
                            }
                        }
                    }
                }
            }

            // We've finally built the full primal!
            // We're ready to move to the solution of either the primal or 
            // the dual. The solution is in the form of a LogoList with the 
            // first element being a list of the solution values of the 
            // primal variables (or shadow prices), and the second element 
            // being the value of the primal's objective function at the
            // solution point.
            
            // NumAnalUtils.printString("prial nonNegative " + nonNegative, context);
            LogoList soln;
            if (solutionType == SOLVE_PRIMAL) {
                soln = SimplexSolver(objective, constraints, goal, nonNegative, context);
            } else {
                soln = SolveDual(objective, constraints, goal, nonNegative, context);
            }
            return soln;
        }
    }

    private static LogoList SimplexSolver(double[] objective, ArrayList<LinearConstraint> constraints,
            GoalType goal, boolean nonNegative, Context context)
            throws ExtensionException {
        // find the solution.  Note that we assume the constant in the
        // objective function is zero, since a constant does not affect
        // the solution, only the final value of the ojective funtion at
        // the solution point.
        PointValuePair solution;
        SimplexSolver solver = new SimplexSolver();
        try {
            solution = solver.optimize(
                    new LinearObjectiveFunction(objective, 0.0),
                    new LinearConstraintSet(constraints),
                    goal,
                    new NonNegativeConstraint(nonNegative));
            // we have a solution. calculate the value of the objective function
            // and return a LogoList with the first element being the solution
            // list and the second the value of the objective function.
            // NOTE that PointValuePair has no method for extracting the value!
            double result = 0.0;
            double point[] = solution.getPoint();
            for (int i = 0; i < objective.length; i++) {
                result += objective[i] * point[i];
            }
            LogoListBuilder lst = new LogoListBuilder();
            lst.add(NumAnalUtils.convertArrayToSimpleLogoList(point));
            lst.add(result);
            return lst.toLogoList();
        } catch (TooManyIterationsException | NoFeasibleSolutionException | UnboundedSolutionException e) {
            throw new ExtensionException(e.toString());
        }
    }

    private static LogoList SolveDual(double[] primalObjective,
            ArrayList<LinearConstraint> primalConstraints,
            GoalType primalGoal, boolean primalNonNegative, Context context)
            throws ExtensionException {

        int nvar = primalObjective.length;
        int ncon = primalConstraints.size();

        boolean primalIsFree[] = new boolean[nvar];
        // now if !primalNonNegative, check for and  
        if (primalNonNegative) {
            // All the primal variables are constrained to be non-negative, 
            // and thus put false entries in the primalIsFree array.
            Arrays.fill(primalIsFree, false);
        } else {
            // Some primal variables may be free.
            // Assume initially that all the primal variables are free and
            // then find those that are explicitly constrained.
            // Remove those constraints on individual variables, keeping 
            // track of those by RemoveNonNegConstraints() putting 
            // corresponding false values in the primalIsFree array. 
            Arrays.fill(primalIsFree, true);
            ncon = RemoveNonNegConstraints(primalConstraints, primalIsFree, context);
            // NumAnalUtils.printString("free primals " + Arrays.toString(primalIsFree), context);
            // on the off chance that all variables were explicitly constrained,
            // i.e., none are free, we can simply set primalNonNegative to true.
            boolean allConstrained = true;
            for (boolean k : primalIsFree) {
                allConstrained = allConstrained & !k;
            }
            if (allConstrained) {
                primalNonNegative = true;
            }
        }

        // clean up the constraints to have the proper relationship to the goal.
        // check too for any equality constraints, and keep track of them.
        boolean equalities[] = new boolean[ncon];
        Arrays.fill(equalities, false);       // this is redundant?
        if (primalGoal == GoalType.MAXIMIZE) {
            // rewrite any GEQ constraints as LEQ. Leave equality constraints alone.
            for (int i = 0; i < ncon; i++) {
                if ((primalConstraints.get(i)).getRelationship() == Relationship.GEQ) {
                    double[] primalCoeffs = primalConstraints.get(i).getCoefficients().toArray();
                    for (int j = 0; j < nvar; j++) {
                        primalCoeffs[j] = primalCoeffs[j] * -1.0;
                    }
                    Relationship primalReln = Relationship.LEQ;
                    double primalValue = primalConstraints.get(i).getValue() * -1.0;
                    LinearConstraint newConstraint = new LinearConstraint(primalCoeffs, primalReln, primalValue);
                    primalConstraints.set(i, newConstraint);
                }
                if ((primalConstraints.get(i)).getRelationship() == Relationship.EQ) {
                    equalities[i] = true;
                }
            }
        } else {
            // rewrite any LEQ constraints as GEQ. Leave equality constraints alone.
            for (int i = 0; i < ncon; i++) {
                if ((primalConstraints.get(i)).getRelationship() == Relationship.LEQ) {
                    double[] primalCoeffs = primalConstraints.get(i).getCoefficients().toArray();
                    for (int j = 0; j < nvar; j++) {
                        primalCoeffs[j] = primalCoeffs[j] * -1.0;
                    }
                    Relationship primalReln = Relationship.GEQ;
                    double primalValue = primalConstraints.get(i).getValue() * -1.0;
                    LinearConstraint newConstraint = new LinearConstraint(primalCoeffs, primalReln, primalValue);
                    primalConstraints.set(i, newConstraint);
                }
                if ((primalConstraints.get(i)).getRelationship() == Relationship.EQ) {
                    equalities[i] = true;
                }
            }
        }

        // now form the matrix of coefficients and transpose it.
        double[][] coeffMatrix = new double[ncon][nvar];
        for (int i = 0; i < ncon; i++) {
            coeffMatrix[i] = primalConstraints.get(i).getCoefficients().toArray();
        }
        double[][] dualCoeffMatrix = new double[nvar][ncon];
        for (int i = 0; i < ncon; i++) {
            for (int j = 0; j < nvar; j++) {
                dualCoeffMatrix[j][i] = coeffMatrix[i][j];
            }
        }

        // now form dual objective.
        double[] dualObjective = new double[ncon];
        // min becomes max, and vice versa.
        GoalType dualGoal = (primalGoal == GoalType.MINIMIZE)
                ? GoalType.MAXIMIZE : GoalType.MINIMIZE;
        for (int j = 0; j < ncon; j++) {
            dualObjective[j] = primalConstraints.get(j).getValue();
        }
        // NumAnalUtils.printArray("dual obj ", dualObjective, context);
        // NumAnalUtils.printString("dual goal " + dualGoal.toString(), context);

        // form the constraints with the relationship depending on the goal.
        ArrayList<LinearConstraint> dualConstraints = new ArrayList<>();
        Relationship dualRelationship = (dualGoal == GoalType.MINIMIZE)
                ? Relationship.GEQ : Relationship.LEQ;
        for (int i = 0; i < nvar; i++) {
            LinearConstraint constr = new LinearConstraint(dualCoeffMatrix[i],
                    dualRelationship, primalObjective[i]);
            dualConstraints.add(constr);
            // NumAnalUtils.printArray("dual constr ", constr.getCoefficients().toArray(), context);
            // NumAnalUtils.printString("dual reln " + constr.getRelationship().toString(), context);
            // NumAnalUtils.printString("dual value " + constr.getValue(), context);
        }

        // now check to see if any of the lamdas should be unconstrained,
        // based on an equality in the primal constraints.
        boolean anyEqualities = false;
        for (boolean eq : equalities) {
            anyEqualities = anyEqualities | eq;
        }
        boolean dualNonNegative = true;
        if (anyEqualities) {
            // we need to constrain only the lamdas associated with inequalities.
            dualNonNegative = false;
            for (int i = 0; i < ncon; i++) {
                if (!equalities[i]) {
                    // this lamda is associated with an inequality constraint.
                    double[] lamdaCoeffs = new double[ncon];
                    lamdaCoeffs[i] = 1.0;
                    LinearConstraint lamdaConstr
                            = new LinearConstraint(lamdaCoeffs, Relationship.GEQ, 0.0);
                    dualConstraints.add(lamdaConstr);
                    // NumAnalUtils.printArray("lamda constr", lamdaConstr.getCoefficients().toArray(), context);
                    // NumAnalUtils.printString("lamda reln " + lamdaConstr.getRelationship().toString(), context);
                    // NumAnalUtils.printString("lamda value " + lamdaConstr.getValue(), context);
                }
            }
        }

        // now check to see if any of the constraints should be made equalities
        // because the corresponding primal variable was unconstrained.
        if (!primalNonNegative) {
            // One or more of the elements in primalIsFree must be true.
            for (int i = 0; i < nvar; i++) {
                if (primalIsFree[i]) {
                    // Get the corresponding constraint, change the relationship,
                    // and reinsert it.
                    double coeffs[] = dualConstraints.get(i).getCoefficients().toArray();
                    double value = dualConstraints.get(i).getValue();
                    Relationship newReln = Relationship.EQ;
                    LinearConstraint newConstr = new LinearConstraint(coeffs, newReln, value);
                    dualConstraints.set(i, newConstr);
                    // NumAnalUtils.printString("dual Constraint " + i + " was changed to an equality.", context);
                }
            }
        }

        // Finally, we can find the solution to the constructed dual.
        // NumAnalUtils.printString("dual nonNegative " + dualNonNegative, context);
        LogoList soln = SimplexSolver(dualObjective, dualConstraints, dualGoal,
                dualNonNegative, context);

        return soln;
    }
}
