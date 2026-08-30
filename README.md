# Estimation of the mean function in sparse longitudinal functional data

R code for the Monte Carlo study of a master's thesis on mean function
estimation for sparse longitudinal functional data.

The data are of the form

$$
Y_{ij}(u) = \mu(u, T_{ij}) + X_i(u, T_{ij}) + \varepsilon_{ij}(u),
$$

where $Y_{ij}(u)$ is a curve observed for subject $i$ at visit time $T_{ij}$,
$\mu(u,t)$ is the mean function to be estimated, $X_i(u,t)$ a subject-specific
process and $\varepsilon_{ij}(u)$ an additive noise process. Each subject is
seen only a few times, at random visit times drawn from a density $g$ that is
not uniform, and each visit yields a whole curve rather than a scalar.

The estimators expand $\mu$ in the time direction on an orthonormal system,

$$
\mu(u,t) \approx \sum_{k=1}^{K} \beta_k(u)\,\phi_k(t),
$$

and estimate each coefficient function by a weighted average of the observed
curves. The study compares the weighting schemes, the choice of the system
$\{\phi_k\}$, and the way $K$ is selected.

## What is compared

Weighting schemes:

1. uniform weights;
2. subject-balanced weights;
3. control-neighbours weights (Monte Carlo, leave-one-out Voronoi);
4. spacing weights, built from order-statistic windows and using no design
   density;
5. plug-in variant of scheme 1 with a leave-one-out kernel estimate of $g$.

External benchmarks:

6. bivariate local-linear smoother, at oracle bandwidths;
7. penalised regression spline in time, at an oracle smoothing parameter.

Orthonormal systems, produced by Gram-Schmidt on the ordered dictionary
$(c_1, p_1, \dots, p_r, c_2, c_3, \dots)$ where $\{c_k\}$ is the half-cosine
system and $p_j(t) = t^j$:

- `B0`: plain half-cosine;
- `B1`: augmented with $t$;
- `B2`: augmented with $t$ and $t^2$.

Mean functions, one per boundary regime: `mu` (boundary derivatives non-zero
and different), `mu_B` (non-zero and equal), `mu_C` (vanishing).

## Repository structure

```text
.
├── README.md
├── script.R
└── outputs/
```

`script.R` runs the whole study and writes everything to the output directory.

## Requirements

R, plus the optional package `tikzDevice`. When it is installed together with a
working LaTeX distribution providing `pgf`, figures are exported as TikZ files;
otherwise the code falls back to PDF.

```r
install.packages("tikzDevice")
tinytex::tlmgr_install(c("pgf", "preview"))
```

## Running the code

Quick run, small Monte Carlo sizes, for checking that everything works:

```bash
SIM_SETTINGS=QUICK SIM_OUTDIR=./outputs Rscript script.R
```

Full run, the one behind the reported results:

```bash
SIM_SETTINGS=FULL SIM_OUTDIR=./outputs Rscript script.R
```

The full run is heavy. Start with `QUICK`.

## Output

The output directory is set by `SIM_OUTDIR`. The script writes `.tex` files
holding TikZ figures and LaTeX tables, which the report includes directly with
`\input`. Two options at the top of `script.R` add other formats:

```r
WRITE_PDF <- FALSE   # PDF copies of the figures
WRITE_CSV <- FALSE   # raw numerical results
```

## Reproducibility

The simulation size is read from the environment,

```r
SETTINGS <- Sys.getenv("SIM_SETTINGS", unset = "FULL")
```

and every replication derives its seed from a single master seed,

```r
BASESEED <- 2026L
```

so that `SIM_SETTINGS=FULL` reproduces the reported numbers exactly.

## What the study looks at

The effect on the integrated squared error of the truncation level $K$, of the
weighting scheme, of the visit-time density and of the choice of orthonormal
system; the decay of the design term for each scheme; the cost of estimating
$g$ rather than knowing it; and the selection of the pair (system, truncation)
by five-fold cross-validation at the subject level.
