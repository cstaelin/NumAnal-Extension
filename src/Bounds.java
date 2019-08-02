package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Command;
import org.nlogo.api.Reporter;
import org.nlogo.core.LogoList;
import org.nlogo.api.LogoListBuilder;

/*
 * This class allows users to set solution bounds for a number of the 
 * multivariate optimizers. Bounds may be set independently on each of 
 * the variables. 
 * 
 * SetBounds takes two lists, the first a list of lower bounds, one for each
 * variable, and the second a list of the upper bounds.
 * ClearBounds clears out any bounds that have been set and sets the bounds
 * arrays to null. This is the signal used by calling programs within the 
 * numanal extension that no bounds have been set.
 * GetBounds returns in a nested list the current lower and upper bounds 
 * as sublists. If no bounds are set, returns two empty sublists: [[][]].
 * GetBoundsDefaults returns in a list the default lower bound (item 0) and the 
 * default upper bound (item 1).
 */
public class Bounds {

  // These are the default upper and lower bounds for the solution 
  // vector elements.
  public static final double UPPER_BOUND_DEFAULT = Double.POSITIVE_INFINITY;
  public static final double LOWER_BOUND_DEFAULT = Double.NEGATIVE_INFINITY;

  // We hold these defaults here for convenience of other classes.
  public static final double RELATIVE_TOLERANCE_DEFAULT = 1.0E-12;
  public static final double ABSOLUTE_TOLERANCE_DEFAULT = 1.0E-12;

  static double[] lowerBounds = null;
  static double[] upperBounds = null;

  private static void clear() {
    lowerBounds = null;
    upperBounds = null;
  }
  
  private static void set(double[] lower, double[] upper) {
    lowerBounds = lower.clone();
    upperBounds = upper.clone();
  }
          

  public static class SetBounds implements Command {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.commandSyntax(new int[]{Syntax.ListType(),
        Syntax.ListType()});
    }

    @Override
    public void perform(Argument args[], Context context)
            throws ExtensionException, LogoException {

      LogoList lowerBoundsList = args[0].getList();
      LogoList upperBoundsList = args[1].getList();
      if (lowerBoundsList.size() != upperBoundsList.size()) {
        throw new ExtensionException(
                "The upper and lower bounds lists must be of the "
                + "same length.");
      }
      
      int nvar = lowerBoundsList.size();
      double[] low = new double[nvar];
      double[] high = new double[nvar];
      for (int j = 0; j < nvar; j++) {
        low[j] = ((Number) lowerBoundsList.get(j)).doubleValue();
        high[j] = ((Number) upperBoundsList.get(j)).doubleValue();
        if (high[j] <= low[j]) {
          throw new ExtensionException(
                  "Upper bound " + j + " is less than or equal to the "
                  + "corresponding lower bound: " + high[j]
                  + " vs. " + low[j]);
        }
      }
      Bounds.set(low, high);
    }
  }

  public static class ClearBounds implements Command {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.commandSyntax(new int[]{});
    }

    @Override
    public void perform(Argument args[], Context context)
            throws ExtensionException, LogoException {

      Bounds.clear();
    }
  }

  public static class GetBounds implements Reporter {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{}, Syntax.ListType());
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      LogoListBuilder lst = new LogoListBuilder();
      if (Bounds.lowerBounds != null) {
        lst.add(NumAnalUtils.convertArrayToSimpleLogoList(Bounds.lowerBounds));
        lst.add(NumAnalUtils.convertArrayToSimpleLogoList(Bounds.upperBounds));
      } else {
        lst.add(NumAnalUtils.convertArrayToSimpleLogoList(new double[0]));
        lst.add(NumAnalUtils.convertArrayToSimpleLogoList(new double[0]));
      }
      return lst.toLogoList();
    }
  }

  public static class GetBoundsDefaults implements Reporter {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{}, Syntax.ListType());
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      LogoListBuilder lst = new LogoListBuilder();
      lst.add(LOWER_BOUND_DEFAULT);
      lst.add(UPPER_BOUND_DEFAULT);
      return lst.toLogoList();
    }
  }
}
