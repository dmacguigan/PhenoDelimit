\name{assignProbs_barplot}
\alias{assignProbs_barplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assignProbs_barplot
}
\description{
Make loading plot of variables on selected discriminant axis. Saves plot as file (svg or pdf). Also returns file with variable contributions and loadings for that axis.
}
\usage{
assignProbs_barplot(wd, clumpp.wd,
                    best.perc.var, best.model.number,
                    plot.type, plot.width, plot.height,
                    axis))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{wd}{
  working directory to create plots
  }
  \item{clumpp.wd}{
  directory containing CLUMPP results
  }
  \item{best.perc.var}{
  which percent retained variance to plot? numeric
  }
  \item{best.model.number}{
  which model number to plot? numeric
  }
  \item{plot.type}{
  save plot as "pdf" or "svg"
  }
  \item{plot.width}{
  width of plot in inches
  }
  \item{plot.height}{
  height of plot in inches
  }
  \item{axis}{
  which discriminant axis to examine? numeric
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
