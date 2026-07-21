# APECS Shiny App

ALS Family History – Monogenic Probability Calculator

This repository contains a Shiny application for estimating the probability of monogenic disease in ALS based on family history patterns observed in simulated pedigrees under Mendelian and complex inheritance models.

The app allows users to specify the number of affected relatives by degree of relationship and, optionally, include frontotemporal dementia (FTD) history. It then returns an estimated posterior probability of monogenic disease, a 95% confidence interval, and a comparison between prior and posterior probabilities.

## Features

- Interactive Shiny web application.
- Supports two models:
  - ALS only
  - ALS + FTD
- Displays:
  - Estimated probability of monogenic disease
  - 95% confidence interval
  - Prior versus posterior probability plot
  - Matching pedigree summary table

## Repository structure

```text
.
├── app.R
├── manifest.json
├── output_ppv_als_grid.csv
├── output_ppv_alsftd_grid.csv
├── src/
│   └── PPV_grid.R
└── www/
    └── APECS_logo.png
    └── APECS_relative_count.svg
```

## File overview

- `app.R` – main Shiny application.
- `manifest.json` – deployment manifest generated with `rsconnect::writeManifest()`.
- `output_ppv_als_grid.csv` – lookup table for the ALS-only model.
- `output_ppv_alsftd_grid.csv` – lookup table for the ALS + FTD model.
- `src/PPV_grid.R` – supporting code used to generate or process lookup grid content.
- `www/APECS.png` – application logo served as a static web asset.

## Running locally

Open R in the project directory and run:

```r
shiny::runApp()
```

Alternatively, from an interactive R session:

```r
setwd("APECS_shiny_app")
shiny::runApp()
```

## Deployment

This repository is prepared for deployment using Posit Connect Cloud or another Shiny-compatible hosting platform.

The deployment manifest can be regenerated with:

```r
install.packages("rsconnect")
rsconnect::writeManifest()
```

## Requirements

The application depends on the following R packages:

- shiny
- readr
- dplyr
- bslib
- plotly
- wesanderson
- rsconnect (for deployment only)

## Purpose

This application was developed to support interpretation of ALS family history patterns using pedigree simulation output generated under simple and complex disease models.

## Author

Paul Beele  
UMC Utrecht

## License

MIT License