# Simple Pathway Visualization

## Intro

Given a *curated* input `csv` file of experimental data, does multiple t-tests to test for differential metabolite abundances between case and control, followed by high-level pathway visualization with `iPath3`.


### Requirements

```
pacman
shiny
shinydashboard
shinycssloaders
ggplot2
MetaboAnalystR
KEGGREST
tidyverse
ggfortify
httr
XML
DT
```

### Files

This repo should only contain:

* `README.md` - This readme file. 
* `sample-data.csv` - A `csv` file of example data.
* `app.R` - the Shiny App itself.
* `help.html` - `html` that's a "help" page, accessed from within the app.
* `about.html` - information about the app, accessed from within the app.
