# Differential-Correlation-Mining
Welcome to the Github repo for the DCM R and Matlab packages! These packages contain two implementations of Bodwin et al.'s *A testing-based approach to the discovery of differentially correlated variable sets* which can be found at https://arxiv.org/abs/1509.08124.

This document will provide some tips and tricks for using the DCM R code successfully.

## Getting Started
Currently this package is not available from the CRAN R package repository. To use it you should download it directly from Github. To do this, you can open up the repo in your favorite web browser and click the green **Clone or Download** button on the right.

Once you have the code saved in a local directory, open up your R IDE and source the directory
```
library(R.utils)
sourceDirectory("~/Differential-Correlation-Mining/DCM")
```
Now you are ready to use the DCM package! If you are using Rstudio and have the `devtools` package installed, you can enable our documentation by using the commands,
```
library(devtools)
document("~/Differential-Correlation-Mining/DCM")
```
Now if you have questions about how to use one of our scripts, feel free to type something like `?DCM` or `?init_DCM` and the documentation should be displayed.

## Running the DCM Algorithm on Data
The entire DCM algorithm can be run by using the `DCM` function found in the DCM.R script. All that is required are two numeric data matrices, say `dat.1` and `dat.2`, both with the same number of rows (variables) and any number of (columns). Calling 
```
DCM(dat.1, dat.2)
```
is all you need to do! Of course, you'll probably be wanting to tweak things to fit your needs and the characteristics of your data. There are many options that can be specified including the maximum number of groups to search for (`max.groups`), the maximum number of iterations to perform in each search (`max.iter`), and size of the initial set (`est.size`). In practice one would commonly want to specify `max.groups` and call 
```
DCM(dat.1, dat.2, max.groups, QN=TRUE, echo=TRUE)
```
Please see the documentation for DCM.R for a complete list of options.

## Validation Features
In our usage of the DCM algorithm, we found that at times we wanted a bit more information than is typically returned from each run to get a better idea of what the algorithm is doing behind the scenes. To this end, we have included a validation switch, `validation`, which when switched to `TRUE`, will create a pdf doc and txt doc that provides plots of p-values, test statistics, standard errors and other relevant quantities from each search. These will be saved in a folder called *DCM_Validation_Search{i}* where *i* is the search number that the plots correspond to, within the specified directory, `validation.dir`. If `validation==TRUE` and no `validation.dir` is specified, these documents will be saved in the current working directory.

## Running the DCM Algorithm on Simulated Data
You can run DCM on simulated data just as was done in the paper by using the `makeSim` function to generate data and then calling `DCM` on it. For example, try
```
sim.output <- makeSim(p=500, k=50, rho=0.6, n1=50, n2=40)
dat1 <- sim.output$dat1
dat2 <- sim.output$dat2
DCM(dat1, dat2)
```
Note that since we use the `mvrnorm` function to generate multivariate normal random variables, making `p` large significantly slows things down.

## Questions?
If you have questions about the details of the this method, please refer to Bodwin et al.'s *A testing-based approach to the discovery of differentially correlated variable sets*. If you have questions about how to use this code, please feel free to email Kevin O'Connor at koconn@live.unc.edu.
