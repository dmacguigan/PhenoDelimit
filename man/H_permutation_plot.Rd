\name{H_permutation_plot}
\alias{H_permutation_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
H_permutation_plot
}
\description{
Plot permuted CLUMPP H' results. Returns plot of obsserved H' minus mean permuted H'. Number inside bar is the H' difference, number on top of bar is p-value from permutation test. 
Also returns plots for each model of permuted H' distribution and observed H' (red dashed line). Saves plot as file (svg or pdf).
}
\usage{
H_permutation_plot(wd, clumpp.data, clumpp.data.permuted, model.numbers, best.perc.var, plot.type, plot.prefix, plot.width, plot.height, sig.threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wd}{
  working directory to create plots, should contain "CLUMPP_permuted" subdirectory
  }
  \item{clumpp.data}{
  data frame from read_clumpp_results function
  }
  \item{clumpp.data.permuted}{
  data frame from read_clumpp_results_permuted function
  }
  \item{model.numbers}{
  vector containing delimitaiton model numbers
  }  
  \item{best.perc.var}{
  which percent retained variance to plot? numeric
  }  
  \item{plot.prefix}{
  prefix for plot name, character
  }
  \item{plot.type}{
  save plot as "pdf", "svg", or "png"
  }
  \item{plot.width}{
  width of plot in inches
  }
  \item{plot.height}{
  height of plot in inches
  }
  \item{sig.threshold}{
  NOT IN USE, will allow user to specify statistical significance threshold
  }
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
