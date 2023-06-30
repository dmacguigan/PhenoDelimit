\name{read_clumpp_results}
\alias{read_clumpp_results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read_clumpp_results
}
\description{
Read in results from CLUMPP miscfiles, return data frame of H' values for each delimitation model
}
\usage{
read_clumpp_results(wd, perc.var, model.number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wd}{
  working directory to store CLUMPP files, should contain a CLUMPP directory created by the dapc_clumpp function
  }
  \item{model.number}{
  vector containing delimitaiton model numbers
  }
  \item{perc.var}{
  vector containing cumulative percentages of variance to retain for discriminant analyses
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
