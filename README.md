# buzzfindr: Automated Detection of Bat Feeding Buzzes

## Description

`buzzfindr` is an R function that automates the detection and classification of feeding buzzes in full-spectrum bat echolocation recordings. It uses the signal detection algorithm from the "bioacoustics" R package developed by Marchal et al. 2022, combined with sequential bandpass filtering.

## Installation

You can install the `buzzfindr` function directly from GitHub using the `devtools` package:

```
devtools::install_github("joelwjameson/buzzfindr")
```


Usage
To use the buzzfindr function, simply call it with your echolocation recordings as input:

```
library(buzzfindr)

#Detect and classify feeding buzzes
path <- "path/to/your/recordings"
detected_buzzes <- buzzfindr(path=path)

#View the results
detected_buzzes
```

License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See the LICENSE file for details.

The buzzfindr function utilizes the "bioacoustics" R package, which is also released under the GNU General Public License v3.0 (GPL-3.0). Please ensure compliance with the terms of the "bioacoustics" package license when using this function.
