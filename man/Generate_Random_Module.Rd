% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_Random_Module.R
\name{Generate_Random_Module}
\alias{Generate_Random_Module}
\title{Generate a random module of size 'N'}
\usage{
Generate_Random_Module(All_KOs, N)
}
\arguments{
\item{All_KOs}{A vector of all annotated features}

\item{N}{The number of elements to be included in the randomly generated module}
}
\value{
a vector containing N random annoated features (without replacement).
}
\description{
This function generates a vector of length 'N' composed of randomly selected annotated features
}
\examples{
Generate_Random_Module(All_KOs,6)
}
