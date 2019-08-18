# NetLogo NumAnal extension

The NumAnal extension contains methods for finding the roots of single equations (Brent) and multivariable systems of equations (Newton and Broyden), for finding the minima of single equations (Brent) and multivariable functions (BOBYQA, CDS, CGS, CMAES, DES and Simplex), for finding fixed points (Scarf), and for evaluating definite integrals (Romberg).  

To install this extension, simply put `extensions [ numanal ]` as the first line of your model code (or add `numanal` to the current extensions line) and the Extension Manager of NetLogo 6.1 should install it.  

The documentation for the numanal extension is contained in the PDF file *NumAnal-3.3.0.pdf* and examples of its use are in the `Examples` directory.  Java source files, a manifest and makefile for those who want examine the code or modify the extension for their own purposes are available at https://github.com/cstaelin/NumAnal-Extension/releases.

## Feedback? Bugs? Feature Requests?

Please visit the [github issue tracker](https://github.com/cstaelin/NumAnal-Extension/issues) to submit comments, bug reports, or feature requests.  I'm also more than willing to accept pull requests.

## Credits

The NumAnal extension was written by Charles Staelin, but is based on several freely available numerical analysis libraries. See the documentation for full credits.

## Terms of Use

[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/)

The NetLogo numanal extension is in the public domain.  To the extent possible under law, Charles Staelin has waived all copyright and related or neighboring rights.  However, note that the external libraries may be subject to other license terms. For more information please refer to the `licence.md` file accompanying this release and to  <http://unlicense.org>.
