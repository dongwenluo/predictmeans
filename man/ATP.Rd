\name{ATP}
\alias{ATP}
\docType{data}
\title{\code{ATP} containing data}

\description{\code{ATP} containing data from an experiment to study the effects 
of preserving liquids on the enzyme content of dog hearts. There were 23 hearts 
and two treatment factors, A and B, each at two levels. Measurements were made 
of ATP as a percentage of total enzyme in the heart, at one and two hourly
intervals during a twelve hour period following initial preservation.
}
\usage{
data(ATP)
}
%- maybe also `usage' for other objects documented here.

\format{\code{ATP} is a 230 row data frame with the following columns 
\describe{
\item{heart}{dog heart id.}
\item{time}{time in hour for \code{ATP} measurement.}
\item{A}{treatment with two levels.}
\item{B}{treatment with two levels.}
\item{ATP}{percentage of total enzyme in the heart.}
}
} 

