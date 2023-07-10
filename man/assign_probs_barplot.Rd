\name{assign_probs_barplot}
\alias{assign_probs_barplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assign_probs_barplot
}
\description{
Make "structure-like" bar plots using discriminant analysis assignment probabilities for chosen model and proportion of retained variance. Returns plot and saves plot as file (svg or pdf). Plots will be named based on model selected.
}
\usage{
assign_probs_barplot(wd, clumpp.wd, sample.names,
                    sample.plot.groups, sample.plot.groups.order,
                    sample.order,
                    best.perc.var, best.model.number,
                    plot.type, plot.width, plot.height, colors, border.color, label.size,
                    apriori))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{wd}{
  working directory to create plots
  }
  \item{clumpp.wd}{
  directory containing CLUMPP results
  }
  \item{sample.names}{
  vector of sample names, same order as original data, numeric or character
  }
  \item{sample.plot.groups}{
  vector of sample groups, same order as original data, numeric or character
  }
  \item{sample.plot.groups.order}{
  vector of sample groups in desired order, numeric or character
  }
  \item{best.perc.var}{
  which percent retained variance to plot? numeric
  }
  \item{best.model.number}{
  which model number to plot? numeric
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
  \item{colors}{
  vector of colors for bar plot bars
  }
  \item{border.color}{
  color for borders on bar, character
  }
  \item{cex.axis}{
  font size for x axis sample names, default = 0.3
  }
  \item{apriori}{
  do you wish to plot results from apriori individual assignment to species/populations/clusters, or results from assignment using using k-means clustering? TRUE or FALSE
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
