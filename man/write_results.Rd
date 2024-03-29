% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/9b-write_excel_model_outputs.R
\name{write_results}
\alias{write_results}
\title{Write excel file of Rceattle results}
\usage{
write_results(Rceattle, file = "Rceattle_results.xlsx")
}
\arguments{
\item{Rceattle}{Single or list of Rceattle model objects exported from \code{\link{Rceattle}}}

\item{file}{name of a file to identified the files exported by the
function.}
}
\description{
Write excel file of Rceattle results
}
\examples{

# Load package and data
library(Rceattle)
data(BS2017SS) # ?BS2017SS for more information on the data

# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- fit_mod(data_list = BS2017SS,
   inits = NULL, # Initial parameters = 0
   file = NULL, # Don't save
   debug = 0, # Estimate
   random_rec = FALSE, # No random recruitment
   msmMode = 0, # Single species mode
   avgnMode = 0,
   silent = TRUE)

write_results(Rceattle = BS2017SS, file = 'Rceattle_results.xlsx')
}
