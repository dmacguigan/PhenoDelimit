\name{plot_discriminant_axes}
\alias{plot_discriminant_axes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_discriminant_axes
}
\description{
Make density plot (2 clusters) or scatter plot (>2 clusters) of discriminant axes, as specified by the user. Returns plot and saves plot as file (svg or pdf).
}
\usage{
plot_discriminant_axes(wd, clumpp.wd,
                       sample.plot.groups = NULL, sample.plot.groups.order = NULL,
                       best.perc.var, best.model.number,
                       plot.type, plot.width, plot.height, colors, shapes,
                       x.axis=NULL, y.axis=NULL))
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
  vector of colors for groups
  }
  \item{shapes}{
  vector of shapes for K-means clusters
  }
  \item{x.axis}{
  which discriminant axis to show on x axis? If only two clusters, discriminant axis 1 is automatically used. numeric
  }
  \item{y.axis}{
  which discriminant axis to show on y axis? If only two clusters, this argument is not used. numeric
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
