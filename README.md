# AnalyteR
An R package for estimation of analytes reference intervals.

## Instalation

Install `AnalyteR` from github with `devtools`:

```{r, eval = F}
install.packages("devtools")
devtools::install_github("rdosreis/AnalyteR")
```

## Examples

### Load `AnalyteR`

```{r, eval = F}
library(AnalyteR)
```

### Analyte distribution plots

```{r, eval = F}
data(ca)
AnalyteDistPlot(ca$value)
AnalyteDistPlot(ca$value, analyte.name = "Calcium (mg/dL)")
```

```{r, eval = F}
data(alt)
AnalyteDistPlot(alt$value[alt$sex == "m"],
analyte.name = "Alanine Aminotrasferase (U/L)", main = "Male", probability = FALSE)
```

### Reference intervals

```{r, eval = F}
data(ca)
RefInterval(data = ca, x = "value", analyte.name = "Calcium (mg/dL)")
RefInterval(data = ca, x = "value",
type.ri = "non-parametric", analyte.name = "Calcium (mg/dL)",
partition = "sex")
```

