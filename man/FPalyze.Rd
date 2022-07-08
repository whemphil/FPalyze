\name{FPalyze}
\alias{FPalyze}
\title{Analysis of FP-Based Kd and Off-Rate Experiments}
\description{
Function to analyze Kd and off-rate data from fluorescence polarization experiments.
}
\usage{
FPalyze(experiment.type,path.to.file='./',file.name=c('par.txt','perp.txt'),enzyme='Predator',prey.molecule='Prey',decoy.molecule='Decoy',parameters.file=NULL,save.data=T,save.name='Results',plot.pdf=F,plot.name='Graphs',save.console=T,save.id='Console',variant.concentrations=NULL,data.size='full',time.step=30,t.zero=1,incubation.time=NULL,coerce.timepoints=F,equilibrium.points=10,outliers='none',scale.data=F,default.mP.values=F,background.subtraction=T,G.factor=1,estimate.initials=T,reverse.timepoints=F,regression.approach='both',estimated.kd=NULL,regress.data=T,coerce.offrates=T,start.kn1=NULL,start.ktheta=NULL,start.Bhalf=NULL,start.N=NULL,match.optimizers=T,manual.fEP.adjust=F,FP.baseline=NULL,legend.location='bottomright')
}
\arguments{
\item{experiment.type}{A single character string indicating the type of data to be analyzed; can be 'Kd' or 'COMP'.}
\item{path.to.file}{A single character string indicating the path to the directory where input file(s) are located. Output files will also go to this directory. DEFAULT=(working directory)}
\item{file.name}{A character string vector of length 1 or 2, indicating the names of the raw input txt files. If the vector is of length 2, as expected when default.mP.values=F, the vector is expected to identify the parallel THEN perpindicular intensity files. DEFAULT=c('par.txt','perp.txt')}
\item{enzyme}{A single character string indicating the name of the protein tested. DEFAULT='Predator'}
\item{prey.molecule}{A single character string indicating the name of the prey molecule tested. DEFAULT='Prey'}
\item{decoy.molecule}{A single character string indicating the name of the decoy molecule tested. DEFAULT='Decoy'}
\item{parameters.file}{A single character string indicating the name of a file containing explicit analysis parameters. DEFAULT=NULL.}
\item{save.data}{A logical indicating if the generated R data should be written to an output file. DEFAULT=T}
\item{save.name}{A single character string indicating the name of the output file if save.data=T. DEFAULT='Results'}
\item{plot.pdf}{A logical indicating if the generated plots should be written to an output file. DEFAULT=F}
\item{plot.name}{A single character string indicating the name of the output file if plot.pdf=T. DEFAULT='Graphs'}
\item{save.console}{A logical indicating if the generated console information should be written to an output file. DEFAULT=T}
\item{save.id}{A single character string indicating the name of the output file if save.console=T. DEFAULT='Console'}

Input File Parameters:

\item{variant.concentrations}{what are the concentrations of the variant molecule, nM; DEFAULT -- 2-fold dilutions from 10µM if experiment.type='COMP', and from 1µM if experiment.type='Kd'}
\item{data.size}{is the experiment size half (2 replicates) or full (4 replicates), half/full; DEFAULT -- 'full'}
\item{time.step}{frequency of data collection, s; DEFAULT -- 30}
\item{t.zero}{estimated start-read delay time, min; DEFAULT -- 1}
\item{incubation.time}{window of data collection, min; DEFAULT -- 30 (experiment.type='Kd') or 120 (experiment.type='COMP')}
\item{coerce.timepoints}{should these time-points be forced if they don't match the raw data, T/F; DEFAULT -- F}
\item{equilibrium.points}{average the last __ data points for equilibrium calculations; DEFAULT -- 10}
\item{outliers}{what are the outlier samples in your experiment, none/A1:P24; DEFAULT -- 'none'}
\item{scale.data}{should the mP values be normalized to a 0-1 scale, T/F; DEFAULT -- F}
\item{default.mP.values}{are you using pre-generated mP values, T/F; DEFAULT -- F}
\item{background.subtraction}{should background subtraction be applied, T/F; DEFAULT -- T}
\item{G.factor}{value that the G-factor should be set to, domain=(x ε R | 0 < x ≤ Inf); DEFAULT -- 1}
\item{estimate.initials}{should initial regression values be estimated from data, T/F; DEFAULT -- T}
\item{reverse.timepoints}{should equilibrium points be sampled from the beginning instead, T/F; DEFAULT -- F}
\item{regression.approach}{what regression method should be used for equilibrium binding curve, std/hill/both/none; DEFAULT -- 'both'}
\item{estimated.kd}{estimated Kd for regression, nM; DEFAULT -- NULL}
\item{regress.data}{should the raw data be fit with regression and further modeled, T/F; DEFAULT -- T}
\item{coerce.offrates}{should negative off-rates be coerced to zero, T/F; DEFAULT -- T}
\item{start.kn1}{starting estimate of k-1 for null/displacement modeling, 1/s; DEFAULT -- NULL}
\item{start.ktheta}{starting estimate of kø for displacement modeling, 1/M/s; DEFAULT -- NULL}
\item{start.Bhalf}{starting estimate of normalization constant 1 for optimized model, nM; DEFAULT -- NULL}
\item{start.N}{starting estimate of normalization constant 2 for optimized model; DEFAULT -- NULL}
\item{match.optimizers}{should the optimization parameter values for null model be constrained to those solved for the displacement model, T/F; DEFAULT -- T}
\item{manual.fEP.adjust}{should the baseline FP value for fraction of initial EP calculations be explicitly provided, T/F; DEFAULT -- F}
\item{FP.baseline}{what is the baseline FP value associated with no EP formation, mP; DEFAULT -- NULL}
\item{legend.location}{location of legend for null/displacement models graph; DEFAULT -- 'bottomright'}
}

\value{
Returns a list containing all variables created as part of the function.

Files:

\item{Results.RData}{A compressed R data file containing all variables created as part of the function.}
\item{Graphs.pdf}{A pdf of all graphs produced as part of the function.}
\item{Console.txt}{A text file of all console output generated as part of the function.}
}

\examples{

EXAMPLE 'PARAMETERS.R' FILE:

### Explicit Analysis Definition

# Experiment conditions
variant.concentrations=signif(c(1e4/2^c(0:10),0),3) # what are the concentrations of the variant molecule, nM
data.size='full' # is the experiment size half (2 replicates) or full (4 replicates), half/full

# Data collection parameters
time.step=30 # frequency of data collection, s
t.zero=1 # estimated start-read delay time, min
incubation.time=120 # window of data collection, min
coerce.timepoints=F # should these time-points be forced if they don't match the raw data, T/F

# Analysis parameters
equilibrium.points=10 # average the last __ data points for equilibrium calculations
outliers=c('none') # what are the outlier samples in your experiment, none/A1:P24
scale.data=F # should the mP values be normalized to a 0-1 scale, T/F
default.mP.values=F # are you using pre-generated mP values, T/F
if (default.mP.values==F){
  background.subtraction=T # should background subtraction be applied, T/F
  G.factor=1 # value that the G-factor should be set to, domain={x ε R | 0 < x ≤ Inf}
}
estimate.initials=T # should initial regression values be estimated from data, T/F
# if experiment.type=='Kd'
  reverse.timepoints=F # should equilibrium points be sampled from the beginning instead, T/F
  regression.approach='both' # what regression method should be used for equilibrium binding curve, std/hill/both/none
  if (estimate.initials==F){
    estimated.kd=NULL # estimated Kd for regression, nM
  }
# if experiment.type=='COMP'
  regress.data=T # should the raw data be fit with regression and further modeled, T/F
  if (regress.data==T){
    coerce.offrates=T # should negative off-rates be coerced to zero, T/F
    if (estimate.initials==F){
      start.kn1=NULL # starting estimate of k-1 for null/displacement modeling, 1/s
      start.ktheta=NULL # starting estimate of kø for displacement modeling, 1/M/s
      start.Bhalf=NULL # starting estimate of normalization constant 1 for optimized model, nM
      start.N=NULL # starting estimate of normalization constant 2 for optimized model
    }
    match.optimizers=T # should the optimization parameter values for null model be constrained to those solved for the displacement model, T/F
    manual.fEP.adjust=F # should the baseline FP value for fraction of initial EP calculations be explicitly provided, T/F
    if (manual.fEP.adjust==T){
      FP.baseline=NULL # what is the baseline FP value associated with no EP formation, mP
    }
    legend.location='bottomright' # location of legend for null/displacement models graph
  }

}
