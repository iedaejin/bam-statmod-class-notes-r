# BAM Statistical Modelling — class notes R code

Runnable R code extracted from the **class notes** (bookdown) for Statistical Modelling (Bachelor in Applied Mathematics, IE University).

This is **not** the session/slide examples repo. For those (one folder per class session), see [iedaejin/bam-statmod-r-examples](https://github.com/iedaejin/bam-statmod-r-examples).

One `.R` file per chapter. Source: course `class-notes/*.Rmd`.

## Setup

```r
install.packages(c(
  "tidyverse", "MASS", "faraway", "car", "glmnet", "nlme",
  "lme4", "mgcv", "gamm4", "gratia", "survival",
  "pscl", "brglm2", "leaps", "lmerTest", "pbkrtest",
  "gridExtra"
))
```

From this folder in R:

```r
setwd("path/to/bam-statmod-class-notes-r")
source("_common.R")          # theme + knitr-equivalent options
source("02-statistical-models.R")
```

Or from the shell:

```bash
Rscript -e 'source("_common.R"); source("02-statistical-models.R")'
```

Each chapter script already starts with `source("_common.R")` when the notes do.

## Files

| File | Topic |
|------|--------|
| `_common.R` | Shared ggplot theme and options |
| `01-prerequisites.R` … `20-survival.R` | Chapter examples |
| `B-r-toolkit.R` | R toolkit appendix |

## Notes

- Scripts are meant to be read and run **section by section**, not always end-to-end without packages/data loaded earlier in the same file.
- Some chunks in the notes are `eval=FALSE` (install hints, long demos); they still appear here so you can run them deliberately.
- Regenerated with `knitr::purl` from the class-notes Rmds (`class-notes/extract-r-examples.R` in the course folder).

## Licence

Code: MIT (see `LICENSE`). Course text and slides are not included.
