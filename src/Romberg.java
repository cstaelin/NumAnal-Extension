package org.nlogo.extensions.numanal;

import org.nlogo.api.LogoException;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.Argument;
import org.nlogo.core.Syntax;
import org.nlogo.core.SyntaxJ;
import org.nlogo.api.Context;
import org.nlogo.api.Reporter;
import org.nlogo.api.AnonymousReporter;

public class Romberg {

  /* EPS is the fractional accuracy desired, as determined by the 
   * extrapolation error estimate.  JMAX limits the total number of steps.
   * K is the number of points used in the extrapolation.
   * 
   * THINK ABOUT WHAT WE WANT TO USE FOR ACCURACY CHECKING AND ITS 
   * RELATIONSHIP TO JMAX.
   */
  static final double EPS = 1.0e-8;
  static final int JMAX = 20;
  static final int JMAXP = JMAX + 1;
  static final int K = 5;

  static double currentEstimate;

  public static class RombergFindIntegral implements Reporter {

    @Override
    public Syntax getSyntax() {
      return SyntaxJ.reporterSyntax(new int[]{Syntax.ReporterType(), Syntax.NumberType(),
        Syntax.NumberType()}, Syntax.NumberType(), 3);
    }

    @Override
    public Object report(Argument args[], Context context)
            throws ExtensionException, LogoException {
      // These store the successive trapezoidal approximations and
      // their relative step sizes.
      double s[] = new double[JMAXP];
      double h[] = new double[JMAXP + 1];

      AnonymousReporter fnctn = args[0].getReporter();
      double a = args[1].getDoubleValue();
      double b = args[2].getDoubleValue();

      h[0] = 1.0;
      for (int j = 0; j < JMAX; j++) {
        s[j] = trapzd(fnctn, a, b, j, context);
        if (j >= K - 1) {
          // once we have K successive estimates, 0 through K-1,
          // we can try the polynomial interpolation. We always use
          // the latest K estimates, estimates j-K+1 through j.
          double[] ss = polint(h, s, (j - K + 1), K, 0.0, context);
          if (Math.abs(ss[1]) <= EPS * Math.abs(ss[0])) {
            return ss[0];
          }
        }
        h[j + 1] = 0.25 * h[j];
      }
      NumAnalUtils.writeToNetLogo(
              "Romberg Integration: Too many steps", false, context);
      return 0.0;
    }

    private static double[] polint(double[] xa, double[] ya,
            int firstPoint, int npts, double x, Context context)
            throws ExtensionException, LogoException {
      double[] ydy = new double[2];
      double[] c = new double[npts];
      double[] d = new double[npts];

      double dif = Math.abs(x - xa[firstPoint]);
      int nsi = 0;
      int nsj = 0;
      // find the index of the closest table entry.
      for (int i = 0, j = firstPoint; i < npts; i++, j++) {
        double dift = Math.abs(x - xa[j]);
        if (dift < dif) {
          nsj = j;
          nsi = i;
          dif = dift;
        }
        // and initialize the c and d arrays with the values from ya.
        c[i] = ya[j];
        d[i] = ya[j];
      }

      double y = ya[nsj];  // the initial approximation to y.
      nsi--;
      double dy = 0.0;
      // for each column in the tableau we loop over the current
      // c's and d's, and update them.
      for (int m = 1; m < npts; m++) {
        for (int i = 0, j = firstPoint; i < (npts - m); i++, j++) {
          double ho = xa[j] - x;
          double hp = xa[j + m] - x;
          double w = c[i + 1] - d[i];
          double den = ho - hp;
          if (den == 0.0) {
            throw new ExtensionException(
                    "Romberg Integration: Identical x values.");
          }
          den = w / den;
          d[i] = hp * den;
          c[i] = ho * den;
        }
        if (2 * nsi < (npts - m)) {
          dy = c[nsi + 1];
        } else {
          dy = d[nsi];
          nsi--;
        }
        y += dy;
      }

      ydy[0] = y;
      ydy[1] = dy;

      return ydy;
    }

    private static double trapzd(AnonymousReporter fnctn, double a,
            double b, int estNum, Context context)
            throws ExtensionException, LogoException {
      /* 
       * trapzd performs increasingly refined estimates of the intergral
       * as it is repeatedly called with increasing values of estNum,
       * 0, 1, 2, ..., n. Estimates from the previous steps are saved
       * in the static variable currentEstimate.  
       * Calling trapzd with estNum = 0 restarts the process. 
       * At each stage for estNum > 0, 2^(estNum-1) new points are 
       * added to the estimate ().
       * 
       */
      if (estNum == 0) {
        // the initial estimate, estNum = 0, is the trapezoid formed by
        // the entire range.
        double fa = NumAnalUtils.getFofX(a, fnctn, context);
        double fb = NumAnalUtils.getFofX(b, fnctn, context);
        currentEstimate = 0.5 * (b - a) * (fa + fb);
        // NumAnalExtension.writeToNetLogo(
        //   "estNum " + estNum + " currentEstimate " + 
        //   currentEstimate, false, context);
      } else {
        double sum = 0.0;
        int newPts = (int) Math.pow(2, (estNum - 1));
        double del = (b - a) / newPts;
        double x = a + 0.5 * del;
        for (int j = 0; j < newPts; j++, x += del) {
          sum += NumAnalUtils.getFofX(x, fnctn, context);
        }
        currentEstimate = 0.5 * (currentEstimate + sum * (b - a) / newPts);
        // NumAnalExtension.writeToNetLogo(
        //   "estNum " + estNum + " " + newPts + " currentEstimate " + 
        //   currentEstimate + " " + sum, false, context);
      }
      return currentEstimate;
    }
  }
}
