# Anisotropic Smoothing
This repo is code and methods for a general class of smoothing estimators for change-point regression functions. Articles associated with this work are:
   - Thompson, J.R.J. (2024) [Iterative Smoothing for Change-point Regression Function Estimation](https://www.tandfonline.com/doi/full/10.1080/02664763.2024.2352759), *Journal of Applied Statistics*, 1–25.

For ease of usage, the methods for this paper have been coerced into the <tt>R</tt> package [nonsmooth](https://cran.rstudio.com/web/packages/nonsmooth/) through the function <tt>alc()</tt>. 

## Experimental fire spread data
The experimental fire data used in the article is associated with the following papers:
 - Thompson, J.R.J., Wang, X.J., & Braun, W.J. (2020) [A mouse model for studying fire spread rates using experimental micro-fires](http://www.jenvstat.org/v09/i06). *Journal of Environmental Statistics*, 9(1), 1-19.
 - Wang, X.J., Thompson, J.R.J., Braun, W.J., & Woolford, D.G. (2019) [Fitting a stochastic fire spread model to data.](https://ascmo.copernicus.org/articles/5/57/2019/) *Advances in Statistical Climatology, Meteorology and Oceanography*, 5(1), 57-66.

The micro-fire imagery data is available through a Github <tt>R</tt> package [firedata](https://github.com/jrjthompson/R-package-firedata). For any questions about data access or otherwise, please contact me at john.thompson@ubc.ca.
