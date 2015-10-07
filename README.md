# Super-resolution microscopy analysis in R

__WARNING:__ This package is still under heavy development and documentation is very incomplete. Better documentation and a vignette will be forthcoming at some point.

Package of tools for doing various types analysis of super-resolution microscopy data. At the moment the package has very limited functionality, but the hope is to expand functionality substantially in the future. Currently the tools are primarily focused around provided Ripley's K type cluster analysis<sup>1</sup> and co-cluster<sup>2</sup> analysis of single molecule super-resolution microscopy data. The code heavily utilises [_spatstat_](http://spatstat.github.io/) and, at this early stage, any potential users are strongly encouraged to get familiar with _spatstat_ before attempting to use this package.

<sup>1</sup>Owen, D. M., Rentero, C., Rossy, J., Magenau, A., Williamson, D., Rodriguez, M., & Gaus, K. (2010). PALM imaging and cluster analysis of protein heterogeneity at the cell surface. _Journal of Biophotonics_, 3(7), 446–454. [http://doi.org/10.1002/jbio.200900089](http://doi.org/10.1002/jbio.200900089)

<sup>2</sup>Rossy, J., Cohen, E., Gaus, K., & Owen, D. M. (2014). Method for co-cluster analysis in multichannel single-molecule localisation data. _Histochemistry and Cell Biology_, 141(6), 605–612. [http://doi.org/10.1007/s00418-014-1208-z](http://doi.org/10.1007/s00418-014-1208-z)

## Installation
You will need to install the supr package from github using `devtools`:

```R
# install.packages("devtools")
devtools::install_github("keithschulze/supr")
```