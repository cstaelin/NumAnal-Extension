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
import org.nlogo.api.LogoListBuilder;

/**
 * @author cstaelin
 */
 
public class DIRECTMinimize {
  

public static class DIRECTSetup implements Reporter {
  
    static final int MAXINT_DEFAULT = 20000;
    static final double MAXAREA_DEFAULT = 10.0e-4;


    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(),
        Syntax.ListType(), Syntax.ListType(), 
        Syntax.NumberType() | Syntax.RepeatableType()},
              Syntax.ListType(), 3);
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      
      // Get the bounds and turn them from LogoLists to arrays.
      LogoList xlist = args[1].getList();
      int nvar = xlist.size();
      double[] lb = NumAnalUtils.convertSimpleLogoListToArray(xlist);
      xlist = args[2].getList();
      double[] ub = NumAnalUtils.convertSimpleLogoListToArray(xlist);
      // Get the estimated minimum value.

      // Get reporter task.
      AnonymousReporter fnctn = args[0].getReporter();
      
      // Get the optional arguments.
      int nargs = args.length;
      double maxArea = MAXAREA_DEFAULT;
      if (nargs > 3) {
        maxArea = args[3].getDoubleValue();
      }
      int maxint = MAXINT_DEFAULT;
      if (nargs > 4) {
        maxint = args[4].getIntValue();
      }

      // call the direct routine for a solution.
      DIRECTSolve ds = new DIRECTSolve();
      DIRECTSolve.BoolArrayReturn ba = ds.direct(nvar, lb, ub, maxArea, 
              maxint, fnctn, context);
      boolean solutionFound = ba.trovato;
      double[] xbest = ba.xbest;

      LogoListBuilder lst = new LogoListBuilder();
      lst.add(NumAnalUtils.convertArrayToSimpleLogoList(xbest));
      lst.add(solutionFound);
      return lst.toLogoList();
    }
  }


}
