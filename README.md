# Simple Pathway Visualization

## Intro

Given a *curated* input `csv` file of experimental data, does multiple t-tests to test for differential metabolite abundances between case and control, followed by high-level pathway visualization with `iPath3`.

### Files

`.gitignore` has been set up such that this repo will only contain:

* `README.md` - This readme file.
* `sample-data.csv` - A `csv` file of example data.
* `app.R` - the Shiny App itself.
* `help.html` - `html` that's a "help" page, accessed from within the app.
* `about.html` - information about the app, accessed from within the app.

Calls to `MetaboAnalystR` write out a bunch of intermediary files, so these get `.gitignored`.

### Methodological Details

lorem ipsum

### To Do

 - split main.R into ui.R and server.R
 - volcano plots, heatmaps, boxplots?
 - MetaboAnalystR pathway enrichment module?
 - Add species filter
 - Download all files button
 - sort significantly changing metabs by FC
 - Dropdown: which group is the denominator
 - ipath needs legend
 - add ipath:coverage and ipath:DE metabs
 - include ANOVA?
 - Check DE logic - is p-value input necessary, or only q-value?
 - Add command history
 - make test data and test suite
 - Allow name_map.csv download

### Technical Details

The app makes `POST` calls to `iPath3` to generate `svg` files.

### Requirements for Local Installation

Verified working versions:

* OSX - 3.6.0

Packages:
```
pacman
shiny
shinydashboard
shinycssloaders
ggplot2
MetaboAnalystR
KEGGREST
tidyverse
ggfortify # for pca
httr
XML
DT
```
### Contributors
Don Teng and Veronica Lui
