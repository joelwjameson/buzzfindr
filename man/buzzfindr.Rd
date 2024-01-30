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
buzzfindr(path =, passes = FALSE, exp = 1, buzzprob = 0.8, clickfiltr = FALSE, channel = "left", minfreq = 20, maxfreq = 45, threshold = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
Path to the files to analyse. Specify as "path = path:/to/files".
}
  \item{passes}{
Takes a vector or dataframe column containing the names of files known to be true bat passes. FALSE if you don't want to specify specific files. Default is FALSE.  Providing buzzfindr with true bat passes reduces the likelihood of detecting false positives.
}
  \item{exp}{
integer. Time expansion factor of the recordings. Default is '1' for real time files. Argument from the 'threshold_detection' function in the 'Bioacoustics' package.
}
  \item{buzzprob}{
numeric value between 0 and 1. Probability threshold for deciding if the detection is a buzz (0 is low, 1 is high). Default is 0.8, which detects almost all buzzes. A high value will miss more buzzes but lead to greater classification accuracy. A lower value will detected more buzzes but will lead to detecting more false positives.
}
  \item{clickfiltr}{
numeric value between 0 and 1. Threshold for removing broadband noise in the file (higher removes more noise but can lead to missed buzzes). Default is FALSE. For files with some broadband noise, a value of 0.3 is ususally sufficient. Note: This function still requires fine-tuning as it can lead to missed buzzes even at a low threshold setting.
}
  \item{channel}{
character. The channel for the recording (default is 'left'). Argument from the 'threshold_detection' function in the 'Bioacoustics' package.
}
  \item{minfreq}{
integer. Minimum frequency for the lowest frequency band. Default is '20'. This can be adjusted to remove unwanted frequencies and to just include frequencies of interest.
}
  \item{maxfreq}{
integer. Maximum frequency for the highest frequency band. Default is '45'. This can be adjusted to remove unwanted frequencies and to just include frequencies of interest.
}
  \item{threshold}{
integer. "Sensitivity of the audio event detection function (peak-picking algorithm) in dB" (Marchal et al. 2022). Default is 4. Higher values lead to fewer detections but fewer false positives and much faster analysis (this is useful for accurately quantifying the loudest buzzes). Argument from the 'threshold_detection' function in the 'Bioacoustics' package.
}
}
\details{
Note about noisy files: When there is alot of noise, there is a greater probability that unwanted detected signals conform to a 'buzz' pattern. This can be minimized by providing buzzfindr with the list of files that are true bat passes (passes=), and by adjusting the detection parameters of the function (i.e. adjusting buzzprob, clickfiltr, minfreq, maxfreq, threshold). It is recommended to try adjusting the threshold argument before the clickfiltr argument when analyzing noisy files as the clickfiltr is still being worked on.
}
\value{
A txt file containing: 1) Time: The time at which the buzz was detected within the file; 2) File: The name of the file containing the buzz; and 3) Probability: The probability, attributed by the model, that the detected signal is a buzz.
}
\references{
\itemize{
\item Jameson, J. W. 2023. Automated detection of feeding buzzes in echolocation recordings of bats. In Review.
\item Marchal J, Fabianek F, Scott C. 2022. Bioacoustics: Analyse Audio Recordings and Automatically Extract Animal Vocalizations. R package version 0.2.8.
}
}
\author{
Joel W. Jameson
}

