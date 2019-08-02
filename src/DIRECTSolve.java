package org.nlogo.extensions.numanal;

/*============================================================================================
 !    DIRECT - C implementation of the well-known DiRect Algorithm. A 
 !    Derivative-Free algorithm for bound constrained global optimization problems 
 !    proposed by Jones et al. (see Ref. below)
 !    Copyright (C) 2011  G.Liuzzi, S.Lucidi, V.Piccialli
 !
 !    This program is free software: you can redistribute it and/or modify
 !    it under the terms of the GNU General Public License as published by
 !    the Free Software Foundation, either version 3 of the License, or
 !    (at your option) any later version.
 !
 !    This program is distributed in the hope that it will be useful,
 !    but WITHOUT ANY WARRANTY; without even the implied warranty of
 !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !    GNU General Public License for more details.
 !
 !    You should have received a copy of the GNU General Public License
 !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 !
 !    Jones, D.R., Perttunen, C.D., Stuckman, B.E.: Lipschitzian optimization without 
 !    the Lipschitz constant. J. Optim. Theory Appl. 79(1), 157-181 (1993)
 !
 ==============================================================================================
 */

/*
Some parameters:
nf - the number of function evaluations
nint - the number of hyperrectangle subdivisions.  This is also equal to the number
       of function evaluations and it is used to stop the search for a solution
       if maxint is exceeded.
nconv - within each major iteration, incremented each time a new "vertice" is 
        created.
nelim -
mindiam -
maxdiam - 
*/
import java.util.Arrays;
import org.nlogo.api.Context;
import org.nlogo.api.AnonymousReporter;
import org.nlogo.api.ExtensionException;
import org.nlogo.api.LogoException;

class DIRECTSolve {

  /* ------------------------------------------------------------------------ */
  
  // this is the class defining a hyperrectangle.  two new ones are created
  // at the subdivision of each hyperrectangle into three parts.  they are
  // all chained to each other through "next" and the new ones insert themselves
  // in the chain after the one they are subdividing.  we have added a field
  // "area", which beings at unity for the original hyperrectangle and is 
  // divided by three every every time a hyperrectagle is subdivided.  we have
  // also added a field "depth" which keeps track of the depth of the 
  // subdivision that created this hyperrectangle.  In some sense it is 
  // redundant with area.
  private class Intervallo {

    double[] cent;        // coordinates of the center of the hyperrectangle
    double[] dimen;
    double fint, diam, maxdim, der, area;
    boolean flagdiv, flagloc, flagcon;
    int id, depth;
    Intervallo next;

    Intervallo(int n) {
      cent = new double[n];
      dimen = new double[n];
      next = null;
    }
  }

  private class Vertice {

    Intervallo inter;
    Vertice next;
  }
  
  private class FDir {
    double fdir;
    double[] xdir;
    int fdirID;
    Intervallo fdirHR;
    
    FDir(int n, double fdir, double[] xdir, int fdirID, 
            Intervallo fdirHR) {
      this.fdir = fdir;
      this.xdir = new double[n];
      this.xdir = xdir.clone();
      this.fdirID = fdirID;
      this.fdirHR = fdirHR;
    }
  }

//  private class Fpunt {
//
//    double f;
//    Intervallo punt;
//  }

  private class ModBox {

    int n;
    double[] lb = null;
    double[] ub = null;
    double[] xtemp;
    double[] xbar;

    ModBox(int nn, double[] lower, double[] upper) {
      n = nn;
      lb = new double[n];
      ub = new double[n];
      System.arraycopy(lower, 0, this.lb, 0, n);
      System.arraycopy(upper, 0, this.ub, 0, n);
      xtemp = new double[n];
      xbar = new double[n];
    }
  }

  private class ModSuddividi {

    double[] vetf1;
    double[] vetf2;
    double[] xsud;
    double[] ysud;
    boolean[] mask;

    ModSuddividi(int n) {
      vetf1 = new double[n];
      vetf2 = new double[n];
      xsud = new double[n];
      ysud = new double[n];
      mask = new boolean[n];
    }
  }

  private class IntIntReturn {

    int nf;
    int nint;

    IntIntReturn(int nf, int nint) {
      this.nf = nf;
      this.nint = nint;
    }
  }

  private class VrtIntReturn {

    Vertice vertice;
    int nconv;

    VrtIntReturn(Vertice vertice, int nconv) {
      this.vertice = vertice;
      this.nconv = nconv;
    }
  }

  public class BoolArrayReturn {

    boolean trovato;
    double[] xbest;

    BoolArrayReturn(boolean trovato, double[] xbest) {
      this.trovato = trovato;
      this.xbest = xbest;
    }
  }

  public BoolArrayReturn direct(int n, double[] lb, double[] ub, 
          double maxArea, int maxint, 
          AnonymousReporter funct, Context context) 
          throws ExtensionException, LogoException {

    double[] xdir = new double[n];
    ModBox modBox = new ModBox(n, lb, ub);

//    NumAnalUtils.printString("lb " + Arrays.toString(modBox.lb), context);
//    NumAnalUtils.printString("ub " + Arrays.toString(modBox.ub), context);

    ModSuddividi modSuddividi = new ModSuddividi(n);
    
    // xdir starts at the center of the hypercube.
    for (int i = 0; i < n; i++) {
      xdir[i] = (ub[i] + lb[i]) / 2.0;
    }
 
    Intervallo primo = new Intervallo(n);
    Intervallo start = primo;

    // set the center and dimension of the inital unit hypercube. The center
    // of any hyper-rectangle/hyper-cube is defined relative to the upper 
    // and lower bounds.
    for (int i = 0; i < n; i++) {
      primo.cent[i] = 0.5;
      primo.dimen[i] = 1.0;  // this is a unit hypercube.
    }
    primo.maxdim = 1.0;
    primo.der = 0.0;
    primo.diam = norma(primo.dimen) / 2.0;
    primo.flagloc = false;
    primo.flagdiv = true;
    primo.flagcon = false;
    primo.id = 1;
    primo.area = 1.0;
    primo.depth = 0;

    double[] xbest = unscalevars(modBox, primo.cent);
    double fdir = NumAnalUtils.getFofXvec(xbest, funct, context);
    primo.fint = fdir;
    primo.next = null;
    
    // fDir holds the current minimum function value, it's location relative
    // to the upper and lower bounds, and the hyper-rectangle in which it is
    // located.  It begins by referring to the initial hyper-cube.
    FDir fDir = new FDir(n, fdir, xdir, primo.id, primo);

    int nf = 1;
    int nint = 1;
    int nelim = 1;
    double toldiam = 0.0 * Math.sqrt((double) n) / 2.0; //!!!!!!!!!!!!!!!!!!
    double eps = 1.e-4;
    boolean halt = false;
    boolean trovato = false;

//    NumAnalUtils.printString("xbest " + Arrays.toString(xbest), context);
//    NumAnalUtils.printValue("fdir ", fdir, context);
//    NumAnalUtils.printInt("nf ", nf, context);

    while (!halt) {
      Intervallo curr = start;
      double maxdiam = curr.diam;
      double mindiam = curr.diam;
      double maxdimen = maxVal(curr.dimen);
      int nmaxdimen = 0;

      while (curr != null) {
        if (maxdiam < curr.diam) {
          maxdiam = curr.diam;
        }
        if (mindiam > curr.diam) {
          mindiam = curr.diam;
        }
        double mv = maxVal(curr.dimen);
        if (maxdimen < mv) {
          maxdimen = mv;
          nmaxdimen = 1;
        } else if (mv == maxdimen) {
          nmaxdimen = nmaxdimen + 1;
        }
        curr = curr.next;
      }
      // sets trovato to T/F given the satisfaction of the stopping condition.
      // if trovato is true or we've exceded maxint, break out of the 
      // "while (!halt)" loop. (Setting halt here is redundant.)
      // NOTE: fglob - the known global minimum - was used as the stopping
      // condition, but in our problems, fglob is not known, so is irrelevant
      // as a stopping condition.
      // We will borrow the stopping condition from another implementation.
      // trovato = ((fdir - fglob) / (1.0 > fglob ? 1.0 : fglob)) < tolglob;
      trovato = (fDir.fdirHR.area <= maxArea);
      if ((nint > maxint) || trovato) {
        break;
      }

      // it appears that convex hull refers to the first vertex in the 
      // chain of vertices. 
      VrtIntReturn vi = ricintervallo(start);
      Vertice convexhull = vi.vertice;
      int nconv = vi.nconv;

      nelim = riduciconvexhull(convexhull, nelim, eps, toldiam, fDir.fdir);

      Vertice currch = convexhull;

      // looks like this finds the last vertice in the linked list.
      for (int i = 1; i <= nelim; i++) {
        currch = currch.next;
      }

      while (currch != null) {
        if (currch.inter.flagdiv) {
          IntIntReturn di = suddividi(currch.inter, n, nf, nint, fDir, funct,
                  modBox, modSuddividi, context);
          nf = di.nf;
          nint = di.nint;
        }
        currch = currch.next;
      }

      xbest = unscalevars(modBox, fDir.xdir);

//      String str = "fdir " + fdir + " nf " + nf + " nconv " + nconv + " nelim " + nelim + " mindiam " + mindiam + " maxdiam " + maxdiam;
//      NumAnalUtils.printString(str, context);
    }

//    NumAnalUtils.printString("DIRECT has satisfied the stopping condition.", context);
//    NumAnalUtils.printBoolean("trovato", trovato, context);

    return new BoolArrayReturn(trovato, xbest);
  }

  /* ------------------------------------------------------------------------ */
  private IntIntReturn suddividi(Intervallo curr, int n, int nf, int nint, 
          FDir fDir, AnonymousReporter funct, ModBox modBox,
          ModSuddividi modSuddividi, Context context) {
    
    // Goes through the list of hyperrectangles to find and mark those 
    // eligible for subdivision.  Then calls triplica to subdivide those 
    // indentified into three parts, inserting the new new parts into
    // the list of hyperrectangles.
    
    int i, j;
    int numtrue, ind1, ind2;

    numtrue = 0;
    for (i = 0; i < n; i++) {
      if (curr.maxdim == curr.dimen[i]) {
        System.arraycopy(curr.cent, 0, modSuddividi.ysud, 0, n);
        modSuddividi.ysud[i] = curr.cent[i] + 1.0 * curr.dimen[i] / 3.0;
        modSuddividi.xsud = unscalevars(modBox, modSuddividi.ysud);
        modSuddividi.vetf1[i] = NumAnalUtils.getFofXvec(modSuddividi.xsud, funct, context);

        modSuddividi.ysud[i] = curr.cent[i] - 1.0 * curr.dimen[i] / 3.0;
        modSuddividi.xsud = unscalevars(modBox, modSuddividi.ysud);
        modSuddividi.vetf2[i] = NumAnalUtils.getFofXvec(modSuddividi.xsud, funct, context);
        modSuddividi.mask[i] = true;
        numtrue = numtrue + 1;
        nf = nf + 2;
      } else {
        modSuddividi.vetf1[i] = 1.e+30;
        modSuddividi.vetf2[i] = 1.e+30;
        modSuddividi.mask[i] = false;
      }
    }

    for (i = 1; i <= numtrue; i++) {
      ind1 = minloc(n, modSuddividi.vetf1, modSuddividi.mask);
      ind2 = minloc(n, modSuddividi.vetf2, modSuddividi.mask);
      if (modSuddividi.vetf1[ind1] < modSuddividi.vetf2[ind2]) {
        modSuddividi.mask[ind1] = false;
        triplica(curr, n, ind1, modSuddividi.vetf1[ind1],
                modSuddividi.vetf2[ind1], nint, fDir);
      } else {
        modSuddividi.mask[ind2] = false;
        triplica(curr, n, ind2, modSuddividi.vetf1[ind2],
                modSuddividi.vetf2[ind2], nint, fDir);
      }
      nint = nint + 2;
    }

    // need to return the new values of both nf and nint, and fdir.
    return new IntIntReturn(nf, nint);
  }

  /* ------------------------------------------------------------------------ */
  private void triplica(Intervallo primo, int n, int ind, 
          double f1, double f2, int nint, FDir fDir) {

    // subdivides the hyperrectangle "primo" into three parts. First    
    // create two new hyperrectangle structures and insert them into the chained
    // list of hyperrectagles right after primo.  Then do the subdivision
    // of primo and fill in the new parmeters of all three.
    
    Intervallo secondo = new Intervallo(n);
    Intervallo terzo = new Intervallo(n);
    terzo.next = primo.next;
    secondo.next = terzo;
    primo.next = secondo;

    System.arraycopy(primo.cent, 0, secondo.cent, 0, n);
    System.arraycopy(primo.cent, 0, terzo.cent, 0, n);
    secondo.cent[ind] = secondo.cent[ind] + 1.0 * primo.dimen[ind] / 3.0;
    terzo.cent[ind] = terzo.cent[ind] - 1.0 * primo.dimen[ind] / 3.0;

    System.arraycopy(primo.dimen, 0, secondo.dimen, 0, n);
    System.arraycopy(primo.dimen, 0, terzo.dimen, 0, n);
    primo.dimen[ind] = primo.dimen[ind] / 3.0;
    secondo.dimen[ind] = secondo.dimen[ind] / 3.0;
    terzo.dimen[ind] = terzo.dimen[ind] / 3.0;
    primo.maxdim = maxVal(primo.dimen);
    primo.diam = norma(primo.dimen) / 2.0;
    secondo.maxdim = maxVal(secondo.dimen);
    secondo.diam = norma(secondo.dimen) / 2.0;
    terzo.maxdim = maxVal(terzo.dimen);
    terzo.diam = norma(terzo.dimen) / 2.0;
    
    secondo.flagloc = false;
    terzo.flagloc = false;
    secondo.flagdiv = true;
    terzo.flagdiv = true;
    secondo.flagcon = false;
    terzo.flagcon = false;
    secondo.id = nint + 1;
    terzo.id = nint + 2;
    
    // area and depth keep track of how far down in the subdivisions this set of 
    // hyperrectagles now is.  area it will be (1/3) ^ depth.  area is also the 
    // size of the hyperrectangle relative to the size of the original one.
    primo.area /= 3.0;
    primo.depth += 1;
    secondo.area = primo.area;
    secondo.depth = primo.depth;
    terzo.area = primo.area;
    terzo.depth = primo.depth;
    
    // test the function values of the new hyperrectangles. If either is a 
    // new minimum update fDir.
    secondo.fint = f1;
    if (f1 < fDir.fdir) {
      fDir.fdir = f1;
      System.arraycopy(secondo.cent, 0, fDir.xdir, 0, n);
      fDir.fdirID = secondo.id;
      fDir.fdirHR = secondo;
    }

    terzo.fint = f2;
    if (f2 < fDir.fdir) {
      fDir.fdir = f2;
      System.arraycopy(terzo.cent, 0, fDir.xdir, 0, n);
      fDir.fdirID = terzo.id;
      fDir.fdirHR = terzo;
    }

    // need to return new values of fdir and xdir.
    return;
  }

  /* ------------------------------------------------------------------------ */
  private int riduciconvexhull(Vertice convexhull, int nelim, double eps,
          double toldiam, double fmin) {

    nelim = 0;
    int nugua = 0;
    boolean halt = false;
    Vertice currch = convexhull;

    while (!halt) {
      if (currch != null) {
        if (currch.inter.diam < toldiam) {
          nelim = nelim + 1;
          currch = currch.next;
          continue;
        }
      } else {
        // setting halt here is redundant.
//        halt = true;
        break;
      }

      if (currch.next != null) {
        if ((currch.next.inter.diam - currch.inter.diam) > 0.0) {
          double L = ((((currch.next).inter).fint) - ((currch.inter).fint))
                  / ((((currch.next).inter).diam) - ((currch.inter).diam));
          if (currch.inter.fint - L * currch.inter.diam
                  > fmin - eps * (fmin > 0.0 ? fmin : -fmin)) {

            nelim = nelim + 1 + nugua;
            nugua = 0;
            currch = currch.next;
          } else {
            halt = true;
          }
        } else {
          nugua = nugua + 1;
          currch = currch.next;
        }
      } else {
        halt = true;
      }
    }
    return nelim;
  }

  /* ------------------------------------------------------------------------ */
  private VrtIntReturn ricintervallo(Intervallo root) {

    Vertice convexhull = new Vertice();
    Intervallo curr = root;
    Intervallo primo = root;

    double maxdiam = curr.diam;
    double minfunc = curr.fint;

    /*!----------------------------------------
     ! search for the interval with max. diam.
     ! and min. objective function
     !---------------------------------------- */
    while (curr.next != null) {
      if (((curr.next).fint < minfunc)
              || (((curr.next).diam > maxdiam)
              && ((curr.next).fint - minfunc <= 1.e-9))) {
        primo = curr.next;
        maxdiam = primo.diam;
        minfunc = primo.fint;
      }
      curr = curr.next;
    }
    /*!--------------------------------------
     ! record the first interval belonging
     ! to the convex hull so far identified
     !-------------------------------------- */
    convexhull.inter = primo;
    primo.flagcon = true;
    Vertice currch = convexhull;
    int nconv = 1;
    boolean halt = false;

    while (!halt) {
      /*!-------------------------------------
       ! among those intervals in the upper
       ! right region, find the one with
       ! maximum cosine with vector (1,0)
       !------------------------------------- */
      curr = root;
      double maxcos = -1.0;
      while (curr != null) {
        if ((curr.diam >= (currch.inter).diam) && (!(curr.flagcon))) {
          double norm2 = Math.sqrt(Math.pow(curr.diam - (currch.inter).diam, 2.0)
                  + Math.pow(curr.fint - (currch.inter).fint, 2.0));
          if (norm2 > 0.0) {
            double coseno = (curr.diam - (currch.inter).diam) / norm2;
            if (coseno > maxcos) {
              maxcos = coseno;
              primo = curr;
            }
          } else { // puo' aggiungere anche punti che coincidono sul piano f-d
            maxcos = 1.0;
            primo = curr;
            break;
          }
        }
        curr = curr.next;
      }
      if (maxcos > 0.0) {
        currch.next = new Vertice();
        (currch.next).inter = primo;
        currch = currch.next;
        nconv = nconv + 1;
        primo.flagcon = true;
      } else {
        halt = true;
      }
    }
    currch.next = null;

    //currch = convexhull;
    //if(*convexhull != NULL) printf("convexhull NON vuoto in ricintervallo\n\n");
    return new VrtIntReturn(convexhull, nconv);
  }

  /* ------------------------------------------------------------------------ */
  private static double maxVal(double[] ary) {
    double ris = ary[0];
    for (double d : ary) {
      if (ris < d) {
        ris = d;
      }
    }
    return ris;
  }

  /* ------------------------------------------------------------------------ */
  /*   
   private static double[] scalevars(ModBox modBox, double[] x) {

   double[] y = new double[modBox.n];
   for (int i = 0; i < modBox.n; i++) {
   y[i] = (x[i] - modBox.lb[i]) / (modBox.ub[i] - modBox.lb[i]);
   }
   return y;
   }
   */
  /* ------------------------------------------------------------------------ */
  private static double[] unscalevars(ModBox modBox, double[] y) {

    double[] x = new double[modBox.n];
    for (int i = 0; i < modBox.n; i++) {
      x[i] = modBox.lb[i] + y[i] * (modBox.ub[i] - modBox.lb[i]);
    }
    return x;
  }

  /* ------------------------------------------------------------------------ */
  private static double norma(double[] x) {
    double sum = 0.0;
    for (double d : x) {
      sum = sum + d * d;
    }
    return Math.sqrt(sum);
  }

  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  private int minloc(int n, double[] x, boolean[] mask) {

    int ind = -1;
    double xmin = 1.e+30;
    for (int i = 0; i < n; i++) {
      if ((x[i] < xmin) && mask[i]) {
        xmin = x[i];
        ind = i;
      }
    }
    return ind;
  }

  /* ------------------------------------------------------------------------ */
}
