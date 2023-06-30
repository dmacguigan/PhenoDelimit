\name{dapc_clumpp_permuted}
\alias{dapc_clumpp_permuted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dapc_clumpp_permuted
}
\description{
Permutes the data, then performs performs K-means clustering, runs discriminant analysis, and generates "indfiles" and "paramfiles" for CLUMPP based on user supplied info.
}
\usage{
dapc_clumpp_permuted(wd, data, n.groups, model.numbers, models, perc.var, scale, center, apriori)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wd}{
  working directory to store CLUMPP files, results will be written to new subdirectory "CLUMPP_permuted"
  }
  \item{data}{
  data frame containing phyotypic data, one row per individual, one column per trait, no missing values
  }
  \item{n.groups}{
  vector of the number of groups (populations, species, etc) in each delimitation model
  }
  \item{model.numbers}{
  vector containing delimitaiton model numbers
  }
  \item{models}{
  data frame containing species delimitation models, one row per individual, one column per model
  }
  \item{perc.var}{
  vector containing cumulative percentages of variance to retain for discriminant analyses
  }  
  \item{permutations}{
  number of permutations to perform
  }
  \item{scale}{
  scale PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
  }
  \item{center}{
  center PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
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
