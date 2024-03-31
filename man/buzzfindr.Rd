\name{buzzfindr}
\alias{buzzfindr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
buzzfindr: Automated Detection of Bat Feeding Buzzes
}
\description{
This function couples the signal detection algorithm from the "bioacoustics" package developped by Marchal et al. 2022 with sequential bandpass filtering to detect and classify feeding buzzes in full spectrum bat echolocation recordings
}
\usage{
buzzfindr(path =, passes = FALSE, pass_thresh = FALSE, exp = 1, buzzprob = 0.8, channel = "left", minfreq = 20, maxfreq = 50, out.file = NULL, num_cores = FALSE, vote = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
Path to the files to analyse. Specify as "path = path:/to/files".
}
  \item{passes}{
Takes a vector or dataframe column containing the names of files known to be true bat passes. FALSE if you don't want to specify specific files. Default is FALSE.  Providing buzzfindr with true bat passes reduces the likelihood of detecting false positives.
}
   \item{pass_thresh}{
integer. Argument to specify whether to analyse only recordings that contain potential bat call and the integer value of the threshold parameter that controls the sensitivity of the audio event detection algorithm (dB; 'threshold_detection' function in the 'Bioacoustics' package) (Marchal et al. 2022). Specify FALSE to analyse all recordings in the directory, TRUE to analyse only recordings containing potential bat calls with a default threshold value of 8 dB, or a different threshold value. Default is FALSE.
}
  \item{exp}{
integer. Time expansion factor of the recordings. Default is '1' for real time files. Argument from the 'threshold_detection' function in the 'Bioacoustics' package.
}
  \item{buzzprob}{
numeric value between 0 and 1. Probability threshold for deciding if the detection is a buzz (0 is low, 1 is high). Default is 0.8, which detects almost all buzzes. A high value will miss more buzzes but lead to greater classification accuracy. A lower value will detected more buzzes but will lead to detecting more false positives.
}
  \item{channel}{
character. The channel for the recording (default is 'left'). Argument from the 'threshold_detection' function in the 'Bioacoustics' package.
}
  \item{minfreq}{
integer. Minimum frequency for the lowest frequency band. Default is '20'. This can be adjusted to remove unwanted frequencies and to just include frequencies of interest.
}
  \item{maxfreq}{
integer. Maximum frequency for the highest frequency band. Default is '50'. This can be adjusted to remove unwanted frequencies and to just include frequencies of interest.
}
  \item{out.file}{
Character vector. Format for the results. Options are "csv", "wav", "png". "csv" creates a .csv file of the results with a row for each buzz detection and information about the recording file, path, location of buzz in the file, and the detection parameters and model used. "wav" creates a new .wav file for each detected buzz containing only the buzz. This is useful for easy manual vetting. "png" creates a .png image file of the spectrogram of each buzz This also facilitates manual vetting. Results are output in the specified path for the recordings. Multiple output formats can be obtained by specifying more than one format (e.g. out.file = c("csv","wav")). Default is NULL (i.e. no output to the path directory).
}
  \item{num_cores}{
Argument to specify if parallel processing is desired. Default is FALSE (i.e. no parallel processing). Parallel processing will be implemented if the user specifies the number of threads to use in parallel processing (e.g. num_cores = 8).
}
  \item{vote}{
Argument to apply a voting rule to the final classification decision. Default is FALSE (i.e. no voting). If vote = TRUE, only potential buzz sequences comprised of two or more consecutive buzz segments (a segment is a sequence of four consecutive detection events) will be classified as a true buzz.
}
}
\details{
Note about noisy files: When there is alot of noise, there is a greater probability that unwanted detected signals conform to a 'buzz' pattern. This can be minimized by doing one or more of the following: Specifying a higher classification threshold (buzzprob=), providing buzzfindr with the list of files that are true bat passes (passes=), specifying an initial detection/filtering of files containing a bat call (pass_thresh=), and by applying the voting procedure (vote=).
}
\value{
A dataframe or csv file containing: 1) Time: The time at which the buzz was detected within the file; 2) File: The name of the file containing the buzz; and 3) Probability: The probability, attributed by the model, that the detected signal is a buzz; 4) Path: Path of the recordings analysed; passes: Specified passes argument; pass_thresh: Specified pass_thresh argument; exp: Specified exp argument; buzzprob: Specified buzzprob argument; channel: Specified channel argument; minfreq: Specified minfreq argument; maxfreq: Specified maxfreq argument; model: Model used by buzzfindr for the classification.
}
\references{
\itemize{
\item Jameson, J. W. 2023. buzzfindr: Automating the detection of feeding buzzes in bat echolocation recordings. In Review.
\item Marchal J, Fabianek F, Scott C. 2022. Bioacoustics: Analyse Audio Recordings and Automatically Extract Animal Vocalizations. R package version 0.2.8.
}
}
\author{
Joel W. Jameson
}

