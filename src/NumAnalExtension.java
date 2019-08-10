package org.nlogo.extensions.numanal;
/*
              VERSION 3.3.0 - TO BE USED WITH NETLOGO v6.1
*/ 

/*
  * This extension provides a number of primitives for finding the roots of 
  * both single variable equations and systems of n equations in n variables; 
  * for optimizing (finding the minima) of both univariate and multivariate 
  * equations; and for estimating integrals.  Some were more or less adapted 
  * from Numerical Recipes in C: The Art of Scientific Computing . Others are 
  * dependent upon the Apache Commons Math3 library  or the PAL Math  library.
  * 
  * This extension assumes that functions can be passed to it in the form 
  * of NetLogo anonymous reporters and therefore works only with NetLogo 
  * versions 6.0 and above.  
  * The JAMA Matrix class of java methods is used internally to the extension 
  * for performing linear algebra .  Jama-1.0.3.jar, commons-math3.3.6.1.jar 
  * and PalMathLibrary.jar, all distributed with this package, must be placed 
  * in the same directory as the numanal.jar file.  
  */

public class NumAnalExtension extends org.nlogo.api.DefaultClassManager {

  @Override
  public java.util.List<String> additionalJars() {
    java.util.List<String> list = new java.util.ArrayList<>();
    list.add("Jama-1.0.3.jar");
    list.add("commons-math3-3.2.jar");
    list.add("PalMathLibrary.jar");
    return list;
  }

  // Define the primitives.
  @Override
  public void load(org.nlogo.api.PrimitiveManager primManager) {
    primManager.addPrimitive("simplex",
            new Simplex.SimplexSolve());
    primManager.addPrimitive("Brent-minimize",
            new BrentMinimize());
    primManager.addPrimitive("BrentA-minimize",
            new ApacheBrent.BrentSolve());
    primManager.addPrimitive("Broyden-root",
            new Broyden.BroydenFindRoot());
    primManager.addPrimitive("Broyden-failed?",
            new Broyden.BroydenFailed());
    primManager.addPrimitive("Newton-root",
            new Newton.NewtonFindRoot());
    primManager.addPrimitive("Newton-failed?",
            new Newton.NewtonFailed());
    primManager.addPrimitive("Brent-root",
            new BrentRoot());
    primManager.addPrimitive("Romberg-integrate",
            new Romberg.RombergFindIntegral());
    primManager.addPrimitive("CMAES-minimize",
            new CMAESMinimize.CMAESSolve());
    primManager.addPrimitive("CGS-minimize",
            new CSMinimize.CGSSolve());
    primManager.addPrimitive("CDS-minimize",
            new CSMinimize.CDSSolve());
    primManager.addPrimitive("BOBYQA-minimize",
            new BOBYQAMinimize.BOBYQASolve());
    primManager.addPrimitive("simplex-MD",
            new ApacheSimplex.SimplexSolve(ApacheSimplex.MULTI_DIRECTIONAL_SIMPLEX));
    primManager.addPrimitive("simplex-NM",
            new ApacheSimplex.SimplexSolve(ApacheSimplex.NELDER_MEAD_SIMPLEX));
    primManager.addPrimitive("DES-minimize",
            new DESMinimize.DESSolve());
    primManager.addPrimitive("bounds-set",
            new Bounds.SetBounds());
    primManager.addPrimitive("bounds-get",
            new Bounds.GetBounds());
    primManager.addPrimitive("bounds-get-defaults",
            new Bounds.GetBoundsDefaults());
    primManager.addPrimitive("bounds-clear",
            new Bounds.ClearBounds());
    primManager.addPrimitive("scarfs-fxdpt",
            new ScarfsMethod.getEquilibrium());
    primManager.addPrimitive("scarfs-fxdpt-info", 
            new ScarfsMethod.getScarfsFxdptInfo());
//    primManager.addPrimitive("DIRECT-minimize", 
//            new DIRECTMinimize.DIRECTSetup());
    primManager.addPrimitive("grid-minimize", 
            new GridMinimize());
  }
}
