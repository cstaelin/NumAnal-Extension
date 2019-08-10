package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;
import java.util.ArrayList;
import org.nlogo.api.LogoListBuilder;

/**
 * *****************************************************************************
 * GESOLVER
 *
 * This package computes general equilibrium prices using Scarf's algorithm. See
 * the included PDF for more details.
 *
 * Richard Moshe Katzwer Princeton University 2012
 * ****************************************************************************
 */
/**
 * This class implements Scarf's Simplical Subdivision Algorithm for computing
 * general equilibrium prices
 */
public class ScarfsMethod {

  /* Default paramaters */
  private static final double TOLERANCE = 1.0E-8;
  private static final int MAX_ITERS = 250000;
  private static final long MAX_GRID_SIZE = (long) 1.0E17;
  private static final int GRID_REFINE = 2;
  private static final boolean USE_RELATIVE_ERROR = true;
  private static final boolean KEEP_ITERATION_LOG = false;
  
  public static class getScarfsFxdptInfo implements Reporter {
    
    static ScarfReturnType scfrt;
    
    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{}, Syntax.ListType());
    }
    
    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {
      
      if (scfrt.keepLog) {
        String iterLog = scfrt.getLog(0);
      }
      LogoListBuilder lst = new LogoListBuilder();
      lst.add((double)scfrt.exitFlag);
      lst.add(scfrt.isRelError);
      lst.add((double)scfrt.totalIters);
      lst.add((double)scfrt.totalTime);
      lst.add(scfrt.getExitMessage());
      return lst.toLogoList();
    }
    
    private static void setScrt(ScarfReturnType scrt) {
      scfrt = scrt;
    }
  }
  
  public static class getEquilibrium implements Reporter {
    
    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{Syntax.NumberType(),
        Syntax.ReporterType(),
        Syntax.WildcardType() | Syntax.RepeatableType()},
              Syntax.ListType(), 2);
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {

      /* numerical parameters */
      double tolerance = TOLERANCE;		// maximum tolerance of error in each market
      int maxIters = MAX_ITERS;			// maximum total iterations 
      int gridRefine = GRID_REFINE;			// factor by which we increase fineness of grid
      long maxGridSize = MAX_GRID_SIZE;		// maximum "fineness" of grid over simplex
      // in theory, bounded by maximum size of long data type
      boolean useRelativeError = USE_RELATIVE_ERROR;	// If true, use market clearing error by percentage
      // If false, use magnitude of market clearing error
      boolean keepIterationLog = KEEP_ITERATION_LOG;

      int nargs = args.length;
      // Get the number of variables.
      int nvar = args[0].getIntValue();
      // Get reporter task.
      AnonymousReporter mapping = args[1].getReporter();
      
      // Get optional arguments.
      if (nargs > 2) { useRelativeError = args[2].getBooleanValue(); }
      if (nargs > 3) { keepIterationLog = args[3].getBooleanValue(); }
        
      
      ScarfReturnType rslts = solve(mapping, nvar, tolerance,
              maxIters, maxGridSize, gridRefine, useRelativeError, 
              keepIterationLog, context);
      getScarfsFxdptInfo.setScrt(rslts);
      return NumAnalUtils.convertArrayToSimpleLogoList(rslts.prices);
    }
  }
  
  private static ScarfReturnType solve(AnonymousReporter mapping,
          int numMarkets, double tolerance, int maxIters, long maxGridSize,
          int gridRefine, boolean useRelativeError, boolean keepIterationLog, Context context)
          throws ExtensionException {
    
    ScarfReturnType srt = new ScarfReturnType(mapping);
    
    long startTime = System.currentTimeMillis();

    // the grid on the simplex
    long gridSize = numMarkets;
    
    long[] W = new long[numMarkets + 1];
    W[numMarkets] = 1;
    for (int i = 0; i < numMarkets; i++) {
      W[i] = gridSize / numMarkets;
    }

    // unit vectors
    long[][] unitVecs = getUnitVecs(numMarkets + 1);
    
    long[][] X = new long[numMarkets + 1][numMarkets + 1];
    for (int i = 0; i < numMarkets + 1; i++) {
      X[i] = pivot(W, unitVecs[i], unitVecs[numMarkets]);
    }

    //initialize vertex labels
    double[] label = new double[numMarkets + 1];
    
    int pivotIndex;
    int iter;
    int totalIter = 0;
    int numGridRefines = 0;

    // we will store prices here
    double[] prices = null;
    
    long time_cum = 0;
    long time_step;

    // temp dummy variable
    int y;
    
    int exitFlag = 0;
    
    while (exitFlag == 0) {
      
      iter = 0;
      
      time_step = System.currentTimeMillis();
      
      for (int i = 0; i < numMarkets + 1; i++) {
        X[i] = pivot(W, unitVecs[i], unitVecs[numMarkets]);
      }
      for (int i = 0; i < numMarkets; i++) {
        label[i] = i + 1;
      }
      
      pivotIndex = numMarkets;
      
      while (X[pivotIndex][numMarkets] < 2 && iter < maxIters) {
        
        if (X[pivotIndex][numMarkets] == 0) {
          y = 0;
          for (int j = 0; j < numMarkets + 1; j++) {
            if (X[pivotIndex][j] > 0) {
              y = j;
              break;
            }
          }
          for (int j = 0; j < numMarkets + 1; j++) {
            if (X[pivotIndex][j] > 0 && X[pivotIndex][j] - W[j] > X[pivotIndex][y] - W[y]) {
              y = j;
            }
          }
          label[pivotIndex] = y + 1;
        } else if (X[pivotIndex][numMarkets] == 1) {
          prices = new double[numMarkets];
          for (int j = 0; j < numMarkets; j++) {
            prices[j] = X[pivotIndex][j];
          }
          
          double[][] DandS = NumAnalUtils.getMofXvec(prices, mapping, context);
          
          y = 0;
          for (int j = 0; j < numMarkets; j++) {
            if (X[pivotIndex][j] > 0) {
              y = j;
              break;
            }
          }
          for (int j = 0; j < numMarkets; j++) {
            if (X[pivotIndex][j] > 0 && (DandS[0][j] / DandS[1][j] - 1d) * X[pivotIndex][j] < (DandS[0][y] / DandS[1][y] - 1d) * X[pivotIndex][y]) {
              y = j;
            }
          }
          label[pivotIndex] = y + 1;
        } else {
          label[pivotIndex] = numMarkets + 1;
        }

        //Update E
        for (int i = 0; i < numMarkets + 1; i++) {
          if (i != pivotIndex && label[i] == label[pivotIndex]) {
            pivotIndex = i;
            break;
          }
        }

        //pivot step
        if (pivotIndex == 0) {
          X[pivotIndex] = pivot(X[numMarkets], X[1], X[0]);
        } else if (pivotIndex == numMarkets) {
          X[pivotIndex] = pivot(X[numMarkets - 1], X[0], X[numMarkets]);
        } else {
          X[pivotIndex] = pivot(X[pivotIndex - 1], X[pivotIndex + 1], X[pivotIndex]);
        }
        
        iter++;
        
      }
      
      totalIter += iter;
      
      double[] p = new double[numMarkets];
      for (int i = 0; i < numMarkets + 1; i++) {
        X[pivotIndex][i] = 0;
      }
      for (int i = 0; i < numMarkets; i++) {
        for (int g = 0; g < numMarkets + 1; g++) {
          p[i] += X[g][i];
        }
        p[i] /= (gridSize * numMarkets);
      }
      
      prices = p;
      
      double[][] DandS = NumAnalUtils.getMofXvec(prices, mapping, context);
      
      time_step = System.currentTimeMillis() - time_step;
      time_cum += time_step;
      
      double currentError = useRelativeError ? maxRelativeError(DandS) : maxError(DandS);
      
      if (keepIterationLog) {
        srt.addIter(numGridRefines, prices, gridSize, time_step, iter, time_cum, totalIter);
      }

      if (currentError > tolerance) {
        gridSize *= gridRefine;
        int s = (pivotIndex == 0) ? 1 : 0;
        for (int i = 0; i < numMarkets; i++) {
          W[i] = X[s][i] * gridRefine;
        }
      } else {
        exitFlag = 1;
      }
      
      if (gridSize > maxGridSize) {
        exitFlag = 2;
      }
      
      numGridRefines++;
      
    }
    
    srt.prices = prices;
    srt.tolerance = tolerance;
    srt.isRelError = useRelativeError;
    srt.keepLog = keepIterationLog;
    srt.totalIters = totalIter;
    srt.totalTime = System.currentTimeMillis() - startTime;
    srt.exitFlag = exitFlag;
    return srt;
    
  }

  /**
   * Returns an n x n matrix of unit vectors, i.e. the n x n identity matrix.
   *
   * @param n dimension
   * @return n x n identity matrix
   */
  private static long[][] getUnitVecs(int n) {
    long[][] z = new long[n][n];
    for (int i = 0; i < n; i++) {
      z[i][i] = 1;
    }
    return z;
    
  }

  /**
   * A simple pivoting operation
   *
   * @param x1 a vector
   * @param x2 a vector
   * @param x3 a vector
   * @return x1 + x2 - x3
   */
  private static long[] pivot(long[] x1, long[] x2, long[] x3) {
    long[] q = new long[x1.length];
    for (int i = 0; i < q.length; i++) {
      q[i] = x1[i] + x2[i] - x3[i];
    }
    return q;
  }

  /**
   * Returns the maximum relative error among all markets, where DandS[0] =
   * Demand, and DandS[1] = supply.
   *
   * @param DandS the quantities demanded and supplied
   * @return max over m, abs(Demand(m)/Supply(m)-1)
   */
  private static double maxRelativeError(double[][] DandS) {
    double max = (DandS[0][0] == DandS[0][1]) ? 0 : Math.abs(DandS[0][0] - DandS[1][0]) / DandS[1][0];
    for (int i = 1; i < DandS.length; i++) {
      max = Math.max(max, (DandS[0][i] == DandS[1][i]) ? 0 : Math.abs(DandS[0][i] - DandS[1][i]) / DandS[1][i]);
    }
    return max;
  }

  /**
   * Returns the maximum error among all markets, where DandS[0] = Demand, and
   * DandS[1] = supply.
   *
   * @param DandS the quantities demanded and supplied
   * @return max over m, abs(Demand(m) - Supply(m))
   */
  private static double maxError(double[][] DandS) {
    double max = Math.abs(DandS[0][0] - DandS[1][0]);
    for (int i = 1; i < DandS.length; i++) {
      max = Math.max(max, Math.abs(DandS[0][i] - DandS[1][i]));
    }
    return max;
  }

  /**
   * A wrapper class that contains the solved-for equilibrium prices as well as
   * the log of the algorithm's run.
   */
  private static class ScarfReturnType {
    
    public AnonymousReporter demandSystem;
    public double[] prices;
    public ArrayList<GridRefinementIter> iterList;
    public int exitFlag;
    public long totalTime;
    public int totalIters;
    public double tolerance;
    public boolean isRelError;
    public boolean keepLog;
    
    private ScarfReturnType(AnonymousReporter mapping) {
      iterList = new ArrayList<>();
      this.demandSystem = mapping;
    }
    
    private void addIter(int curStep, double[] prices, long gridSize,
            long curStepTime, int curStepIters, long totalTime,
            int totalIter) {
      
      iterList.add(new GridRefinementIter(demandSystem, curStep, prices,
              gridSize, curStepTime, curStepIters, totalTime, totalIter));
    }
    
//    private void printLog(int digits) {
//      for (GridRefinementIter iterList1 : iterList) {
//        System.out.println(iterList1.toString());
//      }
//    }
    
    private String getLog(int digits) {
      String ret = "";
      for (GridRefinementIter iterList1 : iterList) {
        ret += iterList1.toString();
      }
      return ret;
    }
    
    private String getExitMessage() {
      switch (exitFlag) {
        case 1:
          return (isRelError ? "relative " : "") + "differences reduced to within " + tolerance;
        case 2:
          return "No convergence, grid became too large";
        default:
          return "Unknown return flag";
      }
    }
    
//    private void setPrice(int m, double newPrice) {
//      double mult = newPrice / prices[m];
//      prices[m] = newPrice;
//      for (int i = 0; i < prices.length; i++) {
//        if (i != m) {
//          prices[i] *= mult;
//        }
//      }
//    }
//    
//    private void normalizePrices(double newSum) {
//      double sum = 0;
//      for (int j = 0; j < prices.length; j++) {
//        sum += prices[j];
//      }
//      for (int j = 0; j < prices.length; j++) {
//        prices[j] *= newSum / sum;
//      }
//    }
    
    private static class GridRefinementIter {

      public long gridSize;
      public long curStepTime;
      public int curStepIters;
      public long totalTime;
      public int totalIter;
      public double[] prices;
      public int curStep;
      AnonymousReporter mapping;
      
      public GridRefinementIter(AnonymousReporter mapping, int curStep,
              double[] prices, long gridSize, long curStepTime,
              int curStepIters, long totalTime, int totalIter) {
        this.gridSize = gridSize;
        this.curStep = curStep;
        this.prices = prices;
        this.curStepIters = curStepIters;
        this.curStepTime = curStepTime;
        this.totalIter = totalIter;
        this.totalTime = totalTime;
        this.mapping = mapping;
      }

//      private String currentStepToString(Context context) {
//
//        String curIter = "";
//
//        curIter += ("\n\n--------------------------------------------------------------------------------------");
//        curIter += ("\n  Grid Refinement: " + curStep);
//        curIter += ("\n--------------------------------------------------------------------------------------\n");
//
//        curIter += "\n" + String.format("%6s", "") + "GridSize = " + gridSize;
//        curIter += "\n" + String.format("%6s", "") + "Current Iterations = " + curStepIters;
//        curIter += "\n" + String.format("%6s", "") + "Total Iterations   = " + totalIter;
//        curIter += "\n" + String.format("%6s", "") + "Current Step Time  = " + curStepTime + " ms";
//        curIter += "\n" + String.format("%6s", "") + "Cumulative Time    = " + totalTime + " ms";
//
//        double[][] DandS = NumAnalUtils.getMofXvec(prices, mapping, context);
//        curIter += "\n\n" + demandSystem.getMarketQuantityString(prices, 5);
//
//        return curIter;
//      }
    }
    
  }
  
}
