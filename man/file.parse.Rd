\name{file.parse}
\alias{file.parse}
\title{Analysis File Generator}
\description{
Function to take master csv files for parallel and perpindicular FP experiment data and parse them into individual data files for analysis.
}
\usage{
file.parse(data.rows,background.rows,exp.type='full',par.file='par.csv',perp.file='perp.csv',path.to.file='./')
}
\arguments{
\item{data.rows}{A character string vector that sequentially indicates the plate rows containing experimental/non-blank reactions; ex. c('A','B','C','D') or LETTERS[1:4].}
\item{background.rows}{A character string vector that sequentially indicates the corresponding background/blanking rows for each reaction/row set; ex. if experimental set 1 is in rows A:B and 2 is in C:D, or if 1 is in A and 2 is in B, and experimental sets 1 and 2 both correspond to the blanks in row F, then background.rows=c('F','F').}
\item{exp.type}{A single character string indicating the size of each experimental set; can be 'single' for 1-replicate/0.5-row per set, 'half' for 2-replicates/1-row per set, or 'full' for 4-replicates/2-rows per set. DEFAULT='full'.}
\item{par.file}{A single character string indicating the name of the master csv file containing all parallel intensity data. DEFAULT='par.csv'.}
\item{perp.file}{A single character string indicating the name of the master csv file containing all perpendicular intensity data. DEFAULT='perp.csv'.}
\item{path.to.file}{A single character string indicating the directory the master csv input files are located in. Output files will also go here. DEFAULT=(working directory).}
}

\value{

Files:

\item{parN.txt}{A text file containing the perpendicular intensity values for each of N experimental sets, appropriate for direct input to the FPalyze function.}
\item{perpN.txt}{A text file containing the parallel intensity values for each of N experimental sets, appropriate for direct input to the FPalyze function.}
}

