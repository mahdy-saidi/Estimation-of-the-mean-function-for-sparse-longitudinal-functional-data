##########################################################################
#                           MONTE CARLO STUDY
#   Estimating the mean function in sparse longitudinal functional data
#
#   Layout:  CONFIGURATION -> NUMERICAL CORE -> CANDIDATE BASES CORE
#            -> REGULARITY CORE -> OUTPUT HELPERS -> E1..E15
#            -> VISUALISATION -> DRIVER
#
#   The E-numbers follow the order of the simulation section of the thesis:
#   E1 effective regularity (tab_decay); E2-E4 anatomy of the risk; E5 the
#   feasible truncation rule; E6-E7 the rates; E8-E11 the comparison of the
#   schemes; E12-E13 the augmented basis; E14 the choice of basis; E15 the
#   FULLY DATA-DRIVEN spacing estimator (basis and truncation both by
#   cross-validation) against the local-linear and the penalised-spline
#   benchmarks, one figure per mean.
#
#   Rates are reported as fitted EXPONENTS: the order of decay of a risk (or
#   of a tuning parameter) as a power of the pooled size M.
###########################################################################

# =========================================================================
#  CONFIGURATION
# =========================================================================

suppressPackageStartupMessages(library(splines))   # spline benchmark (E15)

## ---- Mean function  mu(u,t) --------------------------------------------
P   <- 8L                       # number of Fourier terms
a   <- 1;  b <- 1;  cc <- pi    # cc avoids clash with base::c()
MU_FUN <- function(u, t) {
  nu <- length(u)
  S  <- matrix(0, nu, length(t))                 # start from the zero matrix
  for (p in 1:P) {
    # outer(2*pi*p*u, cc*t, "+")[i,j] = 2*pi*p*u[i] + cc*t[j]   (an nu x nt grid)
    S <- S + ((-1)^p / p^2) * sin(outer(2 * pi * p * u, cc * t, "+"))
  }
  # scale column j by (a*t[j] + b*t[j]^2); 'byrow' repeats that vector down rows
  S * matrix(a * t + b * t^2, nu, length(t), byrow = TRUE)
}

## ---- Second mean: boundary derivatives ZERO  (paper: mu_C) --------------
## mu_C = w * mu with w(t) = 16 t^2 (1-t)^2, so w(0)=w(1)=w'(0)=w'(1)=0.
WIN_C    <- function(t) 16 * t^2 * (1 - t)^2
MU_FUN_C <- function(u, t)
  MU_FUN(u, t) * matrix(WIN_C(t), length(u), length(t), byrow = TRUE)

## ---- Third mean: boundary derivatives NON-ZERO but EQUAL (paper: mu_B) --
## mu_B = mu_C + lambda sin(2 pi u) t, so d_t mu_B(u,0) = d_t mu_B(u,1).
MU_B_AMP <- 5
MU_FUN_B <- function(u, t) MU_FUN_C(u, t) + MU_B_AMP * outer(sin(2 * pi * u), t)

## ---- Zero mean ----------------------------------------------------------
## Passed as 'muf' to gen_data() by E9, E14 and E15, which draw the design,
## the subject processes and the noise ONCE per replication and add each of
## the means afterwards; the comparisons are then paired.
MU_ZERO <- function(u, t) matrix(0, length(u), length(t))

## Registry used by the regularity core and by the basis-selection experiments.
## r_pred is the augmentation the theory predicts for each mean.
MEANS <- list(
  mu  = list(f = MU_FUN,   tex = "$\\mu$",   r_pred = 2L),
  muB = list(f = MU_FUN_B, tex = "$\\mu_B$", r_pred = 1L),
  muC = list(f = MU_FUN_C, tex = "$\\mu_C$", r_pred = 0L))

## ---- Observation design: density g(t) of the visit times ---------------
## Visit times are i.i.d. draws with density g on [0,1]. A design is a list
## with four entries:
##   g(t)  : density (must be > 0 on [0,1]);
##   G(t)  : its integral G(t) = \int_0^t g  (the c.d.f.);
##   gmax  : any number >= max_t g(t)  (used by the rejection sampler);
##   tex   : how the density is written inside LaTeX captions.
## Default: the non-uniform cosine density g(t) = 1 + 3 cos(2 pi t)/5.
G_DENSITY <- function(t) 1 + 0.6 * cos(2 * pi * t)
G_CDF     <- function(t) t + (0.6 / (2 * pi)) * sin(2 * pi * t)
G_MAX     <- 1.6
G_MIN     <- 0.4                                # exact min of g on [0,1]
G_TEX     <- "g(t)=1 + 3 \\cos(2\\pi t) / 5"    # how g is written in captions
DES       <- list(g = G_DENSITY, G = G_CDF, gmax = G_MAX, gmin = G_MIN, tex = G_TEX)
## Skewed alternative used by the robustness study (same g_min = 2/5):
DES_SKEW  <- list(g    = function(t) 0.4 + 1.2 * t,
                  G    = function(t) 0.4 * t + 0.6 * t^2,
                  gmax = 1.6, gmin = 0.4,
                  tex  = "g_2(t)=(2+6t)/5")
lambda_m <- 5L                  # Poisson mean of the visit counts
m_cap    <- 15L                 # visit counts truncated to [2, m_cap]
##       (P[Poisson(5) > 15] < 1e-4), so the
##       bounded-visits assumption holds exactly with m_max = m_cap

## ---- Subject-specific process  X_i -------------------------------------
## X_i(u,t) = sum_{l,l'} C[l,l'] cos(l*pi*u) cos(l'*pi*t), with
## C[l,l'] ~ N(0,1) / (l^k_decay * l'^k_prime_decay) (decay = smoothness).
L             <- 100L           # number of modes
k_decay       <- 2L             # decay exponent: l^k
k_prime_decay <- 2L             # decay exponent: l'^{k'}

## ---- Noise  eps_ij(u) = tau(T_ij) * eta_ij(u) ---------------------------
tau     <- 0.5                  # baseline noise standard deviation
TAU_FUN <- function(t) rep(tau, length(t))   ## baseline scale function
##  tau(.) == tau; the robustness study (E11)
##  passes genuinely t-dependent versions
L_prime <- 50L                  # number of terms in eta_ij

## ---- Truncation, kernel bandwidths and spline penalty grids -------------
Kmax     <- 20L                 # largest truncation level K ever considered
LL_HGRID <- exp(seq(log(0.03), log(0.35), length.out = 30))
# logarithmic h_t grid (sparse t-direction)
# over which the bivariate local-linear
# benchmark is given its oracle
LL_HGRID_U <- c(0, exp(seq(log(0.03), log(0.30), length.out = 6)))
# h_u grid, across the functional argument
# u; at h_u = 0 the joint fit reduces to
# the slice-wise (no u-smoothing) smoother
# of Yao-Mueller-Wang, nested as a section
SPL_NIK   <- 20L                # interior knots of the cubic B-spline basis
SPL_LGRID <- exp(seq(log(1e-9), log(1e1), length.out = 24))
# logarithmic grid of penalties lambda over
# which the spline benchmark (E15) is given
# its own empirical oracle

## ---- Regularity core settings ------------------------------------------
## The t-grid must resolve phi_k for every k that enters the decay regression:
## at k = ORC_KSUM = 300 the integrand oscillates 300 times, and 16001 Simpson
## nodes keep the relative error on ||beta_k||^2 below 1e-6 (checked against a
## 64001-node grid). The u-direction carries no oscillation beyond the P
## Fourier modes of mu, so 501 nodes suffice there.
ORC_NG    <- 16001L  # Simpson nodes in t for every coefficient integral (odd)
ORC_NU    <- 501L    # Simpson nodes in u                                (odd)
ORC_KSUM  <- 300L    # coefficients computed exactly (decay regression)
ORC_KFIT  <- 60L     # first index of the decay regression window

## ---- U/T grids for the VISUALISATION panel only ------------------------
## (the per-experiment evaluation grids are in 'cfg' below, keyed Nu_eval/Nt_eval)
Nu      <- 1001L                # fine grid for observed curves
Nu_viz  <- 201L                 # coarser grid for surface plots (u-axis)
Nt_viz  <- 101L                 # grid for surface plots (t-axis)
VIZ_MI  <- 3L                   # observed profiles shown in the right panel
u_fine <- seq(0, 1, length.out = Nu)
u_viz  <- seq(0, 1, length.out = Nu_viz)
t_viz  <- seq(0, 1, length.out = Nt_viz)

## ---- Output toggles ----------------------------------------------------
WRITE_PDF <- FALSE   # TRUE: also save each figure as fig_*.pdf (besides TikZ .tex)
WRITE_CSV <- FALSE   # TRUE: also write the res_*.csv files with the raw numbers
## (If tikzDevice is NOT installed, figures fall back to fig_*.pdf)

## ---- Run size, output folder, master seed ------------------------------
SETTINGS <- Sys.getenv("SIM_SETTINGS", unset = "FULL")  # "QUICK" | "FULL"
OUTDIR <- Sys.getenv("SIM_OUTDIR", unset = "./outputs/")
BASESEED <- 2026L

# create the output folder if needed
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---- TikZ availability -------------------------------------------------
## TikzDevice needs a working LaTeX install WITH the 'pgf' package (which
## provides TikZ). If 'pgf' is missing we fall back to PDF figures so the
## run still completes.
USE_TIKZ <- requireNamespace("tikzDevice", quietly = TRUE)
if (USE_TIKZ) {
  suppressMessages(requireNamespace("tikzDevice"))
  options(tikzDefaultEngine = "pdftex")
  options(tikzMetricsDictionary = file.path(OUTDIR, "tikzMetrics"))   # cache the metrics
  # append the maths packages to whatever tikzDevice loads by default
  options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
                                "\\usepackage{amsmath}\n", "\\usepackage{amssymb}\n"))
  tikz_works <- function() {                       # try one tiny metric computation
    f <- tempfile(fileext = ".tex")
    ok <- tryCatch({
      tikzDevice::tikz(f, width = 2, height = 2, standAlone = FALSE)
      plot.new(); text(0.5, 0.5, "metric test $x_1$"); TRUE
    }, error = function(e) FALSE)
    while (length(dev.list())) try(dev.off(), silent = TRUE)   # close any device left open
    isTRUE(ok)
  }
  if (!tikz_works()) {
    USE_TIKZ <- FALSE
    message("\n[!] tikzDevice cannot compute LaTeX metrics here, so figures will be\n",
            "    written as PDF instead of TikZ. To enable TikZ output, install the\n",
            "    LaTeX 'pgf' package, then re-run:\n",
            "      TeX Live : tlmgr install pgf preview\n",
            "      tinytex  : tinytex::tlmgr_install(c('pgf','preview'))\n",
            "    (Your final LaTeX document also needs \\usepackage{tikz}, which pgf provides.)\n")
  }
}

# Effort levels. Nu_eval/Nt_eval are the EVALUATION grids used to compute the
# integrated risk (distinct from the visualisation grids Nu/Nu_viz/Nt_viz above).
if (SETTINGS == "FULL") {
  cfg <- list(n_rate=c(100,200,500,1000,2000), R_rate=200L, Nu_eval=50L, Nt_eval=50L,
              n_design=c(125,250,500,1000), R_design=150L, K_design=8L,
              n_riskK=500L, R_riskK=100L, R_case12=80L, n_case12=150L,
              n_plug=c(50,100,200,500), R_plug=100L, R_unbias=2000L,
              n_cv=c(100,200,500,1000,2000), R_cv=80L,
              n_rob=200L, R_rob=150L,
              n_box=c(100,500,2000), R_box=150L,
              n_bcv=500L, R_bcv=150L,
              n_cvbox=c(100,500,2000), R_cvbox=100L)
} else {  # QUICK
  cfg <- list(n_rate=c(50,100,200), R_rate=15L, Nu_eval=30L, Nt_eval=30L,
              n_design=c(100,200,400), R_design=40L, K_design=6L,
              n_riskK=200L, R_riskK=10L, R_case12=20L, n_case12=120L,
              n_plug=c(80,160), R_plug=15L, R_unbias=300L,
              n_cv=c(60,120), R_cv=8L,
              n_rob=80L, R_rob=12L,
              n_box=c(50,200), R_box=10L,
              n_bcv=120L, R_bcv=10L,
              n_cvbox=c(50,200), R_cvbox=8L)
}
cat(sprintf("==== SETTINGS=%s | tikzDevice=%s | WRITE_PDF=%s | WRITE_CSV=%s ====\n",
            SETTINGS, USE_TIKZ, WRITE_PDF, WRITE_CSV))


# =========================================================================
#  NUMERICAL CORE  --  the data-generating process and the estimators
# =========================================================================

## simpson_w(): composite Simpson weights on a regular grid of n (odd) nodes
## spanning [lo, hi]. Every coefficient integral uses these.
simpson_w <- function(n, lo = 0, hi = 1) {
  stopifnot(n %% 2L == 1L)
  w <- c(1, rep(c(4, 2), length.out = n - 2L), 1)
  w * (hi - lo) / (3 * (n - 1))
}

## make_Xi(): returns A NEW FUNCTION representing one random draw of the
## subject process X. Calling make_Xi() again gives an independent subject.
## The returned function maps (u, t) to the matrix X(u[i], t[j]).
make_Xi <- function(L_modes = L, k = k_decay, kp = k_prime_decay) {
  # matrix of random coefficients, divided entrywise by the decay weights
  C <- matrix(rnorm(L_modes * L_modes), L_modes, L_modes) /
    outer(seq_len(L_modes), seq_len(L_modes), function(l, lp) l^k * lp^kp)
  function(u, t)
    outer(u, seq_len(L_modes), function(u, l)  cos(l  * pi * u)) %*% C %*%
    t(outer(t, seq_len(L_modes), function(t, lp) cos(lp * pi * t)))
}

## gen_noise(): one draw of the noise field eta on the points u, for 'nc'
## curves (columns). Returns a (length(u))-by-nc matrix.
## The multiplication by the scale tau(T_ij) happens in gen_data().
gen_noise <- function(u, nc, Lp = L_prime) {
  B  <- matrix(rnorm(Lp * nc), Lp, nc) / (seq_len(Lp)^4)  # decaying random coefficients
  Cu <- outer(u, seq_len(Lp), function(u, l) cos(l * pi * u))
  v  <- as.vector((Cu^2) %*% (seq_len(Lp)^(-8)))          # pointwise variance v(u) > 0
  (Cu %*% B) / sqrt(v)                                    # unit variance at every u
}

## sample_T(): draw m visit times with density des$g on [0,1] by rejection
## sampling: propose uniform candidates, keep each with probability
## g(candidate)/gmax; the kept points then follow g exactly.
sample_T <- function(m, des = DES) {
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m)                       # uniform candidate times
    acc <- runif(2 * m)                     # uniform acceptance variables
    out <- c(out, x[acc <= des$g(x) / des$gmax])
  }
  sort(out[1:m])                            # return the first m, in increasing order
}

## phi_mat(): evaluate the first K half-cosine basis functions at the points t
## (the orthonormal system B_0 of L^2([0,1])): phi_1 = 1 and
## phi_k(t) = sqrt(2) cos((k-1) pi t) for k >= 2. Returns length(t)-by-K.
phi_mat <- function(t, K) {
  Pm <- matrix(1, length(t), K)             # first column is phi_1 == 1
  if (K >= 2) for (k in 2:K) Pm[, k] <- sqrt(2) * cos((k - 1) * pi * t)
  Pm
}

## h_of_M(): the spacing window half-width h = ceil(1/2 + 1/2 * ln ln(M+20)).
## It grows so slowly that h = 2 for every realistic sample size.
h_of_M <- function(M) ceiling(0.5 + 0.5 * log(log(M + 20)))

## cn_weights(): the control-neighbors leave-one-out Voronoi weights
## omega_ij^MC = (1 + c_hat_ij - d_hat_ij)/M, computed exactly in 1-D from the
## sorted design ('s' = visit times, 'G' = c.d.f.). d_hat counts how many other
## points have T_ij as nearest neighbour; c_hat is the g-volume they vacate.
cn_weights <- function(s, G) {
  M <- length(s); ord <- order(s); ss <- s[ord]    # ss = sorted times
  mid <- (ss[-M] + ss[-1]) / 2                     # midpoints between neighbours
  Lb <- c(0, mid); Rb <- c(mid, 1)                 # left/right cell boundaries
  Vfull <- G(Rb) - G(Lb)                           # full Voronoi volume of each point
  nn <- integer(M); if (M >= 2) { nn[1] <- 2L; nn[M] <- (M - 1L) }
  if (M >= 3) {
    lg <- ss[2:(M-1)] - ss[1:(M-2)]; rg <- ss[3:M] - ss[2:(M-1)]
    nn[2:(M-1)] <- ifelse(lg <= rg, 1:(M-2), 3:M)  # nearest neighbour of each point
  }
  dh <- tabulate(nn, nbins = M)                    # d_hat: in-degree as nearest neighbour
  newL <- rep(NA_real_, M); if (M >= 2) { newL[2] <- 0; if (M >= 3) newL[3:M] <- (ss[1:(M-2)] + ss[3:M]) / 2 }
  Vl <- ifelse(seq_len(M) > 1, G(Rb) - G(newL), 0)
  newR <- rep(NA_real_, M); if (M >= 2) { newR[M-1] <- 1; if (M >= 3) newR[1:(M-2)] <- (ss[1:(M-2)] + ss[3:M]) / 2 }
  Vr <- ifelse(seq_len(M) < M, G(newR) - G(Lb), 0)
  nna <- (M - 1) - (seq_len(M) > 1) - (seq_len(M) < M)   # # points NOT adjacent to this one
  o <- numeric(M)
  o[ord] <- (1 + (nna * Vfull + Vl + Vr) - dh) / M       # put weights back in original order
  o
}

## beta_weighted(): generic coefficient estimator for a weight vector 'omega':
## beta_hat_k(u) = sum_ij omega_ij Y_ij(u) phi_k(T_ij)/g(T_ij). 'gfun' is the
## density used (true g, or a plug-in estimate). Returns nu-by-K.
beta_weighted <- function(Y, Tp, omega, gfun, K)
  Y %*% (omega * (phi_mat(Tp, K) / gfun(Tp)))

## beta_weighted_basis(): beta_weighted on an arbitrary candidate basis; for
## r = 0 it reproduces beta_weighted exactly, basis_phi(basis_r(0,K),.) being
## the half-cosine system phi_mat(.,K).
beta_weighted_basis <- function(Y, Tp, omega, gfun, co)
  Y %*% (omega * (basis_phi(co, Tp) / gfun(Tp)))

## beta_spacing(): the spacing estimator beta_hat_k(u) = sum_l Y_(l)(u) D_kl,
## with D_kl = (1/2h) \int_{T_(l-h)}^{T_(l+h)} phi_k. It uses NO density g.
## Window endpoints past the ends are obtained by symmetric reflection.
beta_spacing <- function(Y, Tp, K, h) {
  M <- length(Tp); ord <- order(Tp); Ts <- Tp[ord]; Yo <- Y[, ord, drop = FALSE]
  Text <- function(r) {                          # reflected order statistic at rank r
    out <- numeric(length(r)); lo <- r <= 0; hi <- r > M; in_ <- !lo & !hi
    out[in_] <- Ts[r[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - r[lo]]        # reflect at the left endpoint
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - r[hi]]    # reflect at the right endpoint
    out
  }
  a <- Text((1:M) - h); bnd <- Text((1:M) + h)   # window endpoints for every rank
  D <- matrix(0, M, K)
  D[, 1] <- (bnd - a) / (2 * h)                  # phi_1 == 1: the integral is the width
  if (K >= 2) for (k in 2:K)                     # closed form of the cosine integral
    D[, k] <- (1 / (2 * h)) * sqrt(2) / ((k - 1) * pi) *
    (sin((k - 1) * pi * bnd) - sin((k - 1) * pi * a))
  Yo %*% D
}

## kde_g_loo(): LEAVE-ONE-OUT feasible plug-in estimate of g from the observed
## times (to show the cost of not knowing g). Kernel density with reflection at
## 0 and 1, renormalised to integrate to 1, floored away from 0.
## The density entering observation (i,j)'s weight is built WITHOUT T_ij (and
## without its two reflected images), so no observation contributes to its own
## importance weight.
kde_g_loo <- function(Tp) {
  M  <- length(Tp)
  bw <- bw.nrd0(c(-Tp, Tp, 2 - Tp))                 # bandwidth on the reflected sample
  xg <- seq(0, 1, length.out = 512)
  P1 <- outer(xg, Tp, "+")
  cm <- colMeans(dnorm(P1 / bw) + dnorm(outer(xg, Tp, "-") / bw) + dnorm((P1 - 2) / bw))
  meanK <- sum(cm)                                  # grid-mean of the full reflected KDE shape
  S <- outer(Tp, Tp, "+"); Kat <- rowSums(dnorm(S / bw)) + rowSums(dnorm((S - 2) / bw)); rm(S)
  A <- outer(Tp, Tp, "-"); Kat <- Kat + rowSums(dnorm(A / bw)); rm(A)
  self <- dnorm(2 * Tp / bw) + dnorm(0) + dnorm((2 * Tp - 2) / bw)  # point i's 3 images at T_i
  pmax((Kat - self) / (meanK - cm), 0.05)           # leave-one-out, renormalised, floored
}

## beta_plugin_loo(): kernel plug-in coefficients with the leave-one-out density.
## Drop-in replacement for beta_weighted(., ., ., G_DENSITY, .).
beta_plugin_loo <- function(Y, Tp, omega, K)
  Y %*% (omega * (phi_mat(Tp, K) / kde_g_loo(Tp)))

## mu_hat_grid(): assemble the fitted surface mu_hat(u,t) on a t-grid from the
## coefficient matrix 'beta', using the first K terms. Returns nu-by-(len tgrid).
mu_hat_grid <- function(beta, tgrid, K) beta[, 1:K, drop = FALSE] %*% t(phi_mat(tgrid, K))

## gen_data(): simulate ONE data set of n subjects. Returns a list with the
## pooled visit times T, the (nu) x M matrix Y of observed curves, the subject
## label of each observation, and n. Options:
##   signal_only = TRUE -> drop X and noise (used to isolate the design term);
##   unbalanced  = TRUE -> a fraction p_big of subjects get m_big visits, the
##                         rest get 2 (used to contrast Case 1 vs Case 2);
##   m_fixed = m        -> every subject gets exactly m visits
##                         (balanced arm of E8: Case 1 == Case 2 exactly);
##   des, tau_fun       -> design density and noise-scale function
##                         (the robustness study E11 swaps these in);
##   muf = MU_ZERO      -> a pool with no mean at all, to which E9, E14 and E15
##                         add each mean in turn (common random numbers).
gen_data <- function(n, ug, signal_only = FALSE, unbalanced = FALSE,
                     p_big = 0.15, m_big = 20L, m_fixed = NULL,
                     des = DES, tau_fun = TAU_FUN, muf = MU_FUN) {
  Tp <- vector("list", n); Yc <- vector("list", n); subj <- integer(0)
  for (i in seq_len(n)) {
    if (!is.null(m_fixed)) {
      mi <- as.integer(m_fixed)
    } else if (unbalanced) {
      mi <- if (runif(1) < p_big) m_big else 2L
    } else {
      mi <- 0L                                            ## Poisson visits
      while (mi < 2L || mi > m_cap) mi <- rpois(1L, lambda_m)  # truncated to [2, m_cap]
    }
    Ti <- sample_T(mi, des)
    Yi <- muf(ug, Ti)
    if (!signal_only)                                     ## tau(T_ij) scaling
      Yi <- Yi + make_Xi()(ug, Ti) +
      gen_noise(ug, mi) * matrix(tau_fun(Ti), length(ug), mi, byrow = TRUE)
    Tp[[i]] <- Ti; Yc[[i]] <- Yi; subj <- c(subj, rep(i, mi))
  }
  list(T = unlist(Tp), Y = do.call(cbind, Yc), subj = subj, n = n)
}

## make_weights(): the three deterministic/random weight vectors for a data set.
make_weights <- function(d, des = DES) {
  M  <- length(d$T)
  nm <- table(d$subj)[as.character(d$subj)]         # m_i repeated for each observation
  list(unif = rep(1 / M, M),                        # Case 1
       bal  = 1 / (d$n * as.numeric(nm)),           # Case 2
       MC   = cn_weights(d$T, des$G))               # control-neighbors
}

## ise_curve(): integrated squared error of a fit as a function of K = 1..Kup.
## 'MUtrue' is mu on the (ug, tg) grid. Returns a vector of length Kup.
ise_curve <- function(beta, tg, MUtrue, Kup = Kmax)
  sapply(1:Kup, function(K) mean((mu_hat_grid(beta, tg, K) - MUtrue)^2))

## ---- the BIVARIATE LOCAL-LINEAR external benchmark -----------------
## ll2_ise(): ISE of the benchmark over the FULL bandwidth grid. Returns a
## length(hgrid_t) x length(hgrid_u) matrix; the oracle pair (h_t*, h_u*)
## is chosen downstream from the Monte Carlo average, exactly as the
## empirical oracle K* of the series estimators.
ll2_ise <- function(Y, Tp, ug, tg, MUtrue,
                    hgrid_t = LL_HGRID, hgrid_u = LL_HGRID_U) {
  M <- length(Tp); Nt <- length(tg); Nu <- length(ug)
  Dt <- outer(Tp, tg, "-")                          # M x Nt:  T_ij - t_s
  Du <- outer(ug, ug, "-")                          # Nu x Nu: u_r  - u_q
  ISE <- matrix(NA_real_, length(hgrid_t), length(hgrid_u))
  for (a in seq_along(hgrid_t)) {
    ht <- hgrid_t[a]
    Kt <- pmax(1 - (Dt / ht)^2, 0)                  # t-kernel weights (M x Nt)
    T0 <- colSums(Kt); T1 <- colSums(Kt * Dt); T2 <- colSums(Kt * Dt * Dt)
    A0 <- Y %*% Kt                                  # response t-moments,
    A1 <- Y %*% (Kt * Dt)                           # per slice (Nu x Nt each)
    for (b in seq_along(hgrid_u)) {
      hu <- hgrid_u[b]
      if (hu <= 0) {
        ## h_u = 0 section: univariate local-linear in t at each u_r,
        ## i.e. the pooled smoother of Yao-Mueller-Wang per slice.
        den1 <- T0 * T2 - T1 * T1
        F <- (A0 * rep(T2, each = Nu) - A1 * rep(T1, each = Nu)) /
          rep(pmax(den1, 1e-300), each = Nu)
        bad <- which(den1 <= 1e-10)                 # degenerate windows: NW
        if (length(bad))
          F[, bad] <- (A0 / rep(pmax(T0, 1e-300), each = Nu))[, bad]
      } else {
        ## exact joint plane fit of eq. (8), via factorised normal equations
        Ku <- pmax(1 - (Du / hu)^2, 0)              # u-kernel weights (Nu x Nu)
        C1 <- Ku * Du; C2 <- Ku * Du * Du
        U0 <- colSums(Ku); U1 <- colSums(C1); U2 <- colSums(C2)
        R00 <- crossprod(Ku, A0)                    # response (u,t)-moments,
        R10 <- crossprod(C1, A0)                    # one per target point
        R01 <- crossprod(Ku, A1)                    # (Nu x Nt each)
        S11 <- outer(U0, T0); S12 <- outer(U1, T0); S13 <- outer(U0, T1)
        S22 <- outer(U2, T0); S23 <- outer(U1, T1); S33 <- outer(U0, T2)
        dS <- S11 * (S22 * S33 - S23^2) -           # 3x3 Gram determinant
          S12 * (S12 * S33 - S23 * S13) +
          S13 * (S12 * S23 - S22 * S13)
        dM <- R00 * (S22 * S33 - S23^2) -           # Cramer: first column
          S12 * (R10 * S33 - S23 * R01) +       # replaced by responses
          S13 * (R10 * S23 - S22 * R01)
        F <- dM / pmax(dS, 1e-300)                  # a0-hat at every target
        bad <- which(dS <= 1e-10 * pmax(S11, 1)^3)  # degenerate windows: NW
        if (length(bad)) F[bad] <- (R00 / pmax(S11, 1e-300))[bad]
      }
      ISE[a, b] <- mean((F - MUtrue)^2)             # ISE on the eval grid
    }
  }
  ISE
}

## ---- the PENALISED-SPLINE external benchmark -----------------------
## spl_design(): cubic B-spline basis on the pooled visit times, with SPL_NIK
## equidistant interior knots, its evaluation on the t-grid, and the
## second-order difference penalty.
spl_design <- function(Tp, tg, nik = SPL_NIK) {
  kn <- seq(0, 1, length.out = nik + 2L)[-c(1L, nik + 2L)]
  B  <- splines::bs(Tp, knots = kn, degree = 3L, intercept = TRUE,
                    Boundary.knots = c(0, 1))
  Bt <- predict(B, tg)
  q  <- ncol(B)
  Dm <- diff(diag(q), differences = 2L)
  list(B = as.matrix(B), Bt = as.matrix(Bt), Pen = crossprod(Dm), q = q)
}

## spl_ise(): ISE of the spline benchmark over the whole penalty grid; the
## oracle lambda is chosen downstream from the Monte Carlo average, exactly as
## the empirical oracle K* of the series estimators. The fit is slice-wise in
## u (as the local-linear benchmark, which selects h_u = 0), and the normal
## system is factorised ONCE per lambda and applied to all u-slices at once.
spl_ise <- function(Y, Tp, tg, MUtrue, lgrid = SPL_LGRID) {
  ds  <- spl_design(Tp, tg); M <- length(Tp)
  BtB <- crossprod(ds$B)
  YB  <- Y %*% ds$B                                  # Nu x q
  sapply(lgrid, function(lam) {
    A <- tryCatch(chol2inv(chol(BtB + lam * M * ds$Pen + 1e-10 * diag(ds$q))),
                  error = function(e) NULL)
    if (is.null(A)) return(NA_real_)
    mean((YB %*% A %*% t(ds$Bt) - MUtrue)^2)
  })
}


# =========================================================================
#  CANDIDATE BASES CORE  --  the half-cosine system B_0 and its polynomial
#  augmentations B_1, B_2; cf. Section 5.4 and Remark "Augmenting with t alone"
# =========================================================================

## aug_coeffs_basis(K, mode): explicit Gram-Schmidt coefficients of the
## orthonormal basis phi_1..phi_K built from an ordered dictionary, for the
## three candidates B_r compared in the paper:
##
##     mode = "cos"  ->  (c_1, c_2, c_3, ...)          r = dim V = 0   (B_0)
##     mode = "t"    ->  (c_1, t, c_2, c_3, ...)       r = 1           (B_1)
##     mode = "t2"   ->  (c_1, t, t^2, c_2, c_3, ...)  r = 2           (B_2)
##
## At mode = "t2" this is the J = 1 augmented basis of Section 5.4 (r = 2J):
## phi_1 = 1, phi_2 = sqrt3 (2t-1), phi_3 = sqrt5 (6t^2-6t+1) (shifted
## Legendre); for k = 3+m (m >= 1),
## phi_k(t) = P_k(t) + sqrt2 sum_j theta_{k,j} cos(omega_j t), with P_k a
## scalar multiple of phi_2 (m odd) or phi_3 (m even). Returns, per k, the
## polynomial coefficients over {1,t,t^2} and the cosine part (frequencies
## omega_j and coefficients theta_{k,j}); the amplitude sqrt2 is applied at
## evaluation/integration time. The rank-one recursion uses the partial sums
## A_i = 1 - sum_{j<=i} a_j^2 and B_i = 1 - sum_{j<=i} b_j^2 with
## a_j = <c_{2j},phi_2> = -4 sqrt6/((2j-1)^2 pi^2) and
## b_j = <c_{2j+1},phi_3> = 3 sqrt10/(j^2 pi^2).
##
## The "t" branch is the one to read. Write k = 2 + m, m >= 1, for the position
## of the cosine c_{m+1} in the augmented order. Two things happen, according to
## the parity of m, and they are exactly the two blocks of Lemma cos_ibp:
##
##  * m EVEN. c_{m+1} is symmetric about t = 1/2, phi_2 (propto t - 1/2) is
##    antisymmetric, and every earlier phi is either symmetric with a different
##    cosine frequency or antisymmetric. So c_{m+1} is ALREADY orthogonal to
##    everything before it: phi_k = c_{m+1}, with no polynomial part.
##
##  * m ODD. c_{m+1} is antisymmetric and is NOT orthogonal to phi_2. The
##    Gram-Schmidt step is then the a_j / A_i one above, and the resulting
##    phi_k is the odd-index element of B_2, verbatim.
##
## Consequences used in the paper: every element of B_1 is either a cosine or an
## element of B_2, so Lemma uniform_bound applies to it without change; and the
## closed forms below are exact, so the window integrals D_kl stay closed-form
## for all three candidates.
aug_coeffs_basis <- function(K, mode = c("t2", "t", "cos")) {
  mode <- match.arg(mode)
  poly   <- matrix(0, K, 3)          # rows: coefficients of phi_k over 1, t, t^2
  freqs  <- vector("list", K)        # freqs[[k]] : cosine frequencies of phi_k
  thetas <- vector("list", K)        # thetas[[k]]: matching amplitudes
  poly[1, ] <- c(1, 0, 0)            # phi_1 = 1 in every candidate
  
  ## the two orthonormalised monomials, as polynomials in (1, t, t^2)
  P2 <- c(-sqrt(3),  2 * sqrt(3), 0)             # phi_2 = sqrt3 (2t - 1)
  P3 <- c( sqrt(5), -6 * sqrt(5), 6 * sqrt(5))   # phi_3 = sqrt5 (6t^2 - 6t + 1)
  
  ## ---- r = 0: the plain half-cosine basis, phi_k = c_k --------------------
  if (mode == "cos") {
    if (K >= 2) for (k in 2:K) { freqs[[k]] <- (k - 1) * pi; thetas[[k]] <- 1 }
    return(list(poly = poly, freqs = freqs, thetas = thetas,
                K = K, mode = mode, r = 0L))
  }
  
  ## ---- the two Gram-Schmidt recursions of Section "Explicit formulae" ----
  jm <- max(K, 4L)
  j  <- seq_len(jm)
  a  <- -4 * sqrt(6)  / ((2 * j - 1)^2 * pi^2)   # a_j = <c_{2j},   phi_2>
  b  <-  3 * sqrt(10) / (j^2 * pi^2)             # b_j = <c_{2j+1}, phi_3>
  A  <- c(1, 1 - cumsum(a^2))                    # A[i+1] = A_i (A_0 = 1)
  B  <- c(1, 1 - cumsum(b^2))                    # B[i+1] = B_i
  
  ## ---- r = 1: dictionary (c_1, t, c_2, c_3, ...) -------------------------
  if (mode == "t") {
    if (K >= 2) poly[2, ] <- P2
    if (K >= 3) for (m in seq_len(K - 2)) {
      k <- 2 + m
      if (m %% 2 == 0) {                          # symmetric: nothing to remove
        freqs[[k]]  <- m * pi
        thetas[[k]] <- 1
      } else {                                    # antisymmetric: a_j / A_i
        nu <- (m + 1) %/% 2
        th <- numeric(nu)
        th[nu] <- sqrt(A[nu] / A[nu + 1])
        if (nu >= 2) th[seq_len(nu - 1)] <-
          a[nu] * a[seq_len(nu - 1)] / sqrt(A[nu] * A[nu + 1])
        poly[k, ]   <- (-a[nu] / sqrt(A[nu] * A[nu + 1])) * P2
        freqs[[k]]  <- (2 * seq_len(nu) - 1) * pi
        thetas[[k]] <- th
      }
    }
    return(list(poly = poly, freqs = freqs, thetas = thetas,
                K = K, mode = mode, r = 1L))
  }
  
  ## ---- r = 2: dictionary (c_1, t, t^2, c_2, c_3, ...) --------------------
  if (K >= 2) poly[2, ] <- P2
  if (K >= 3) poly[3, ] <- P3
  if (K >= 4) for (m in seq_len(K - 3)) {
    k <- 3 + m
    if (m %% 2 == 1) {                            # odd m -> phi_2, odd freqs
      nu <- (m + 1) %/% 2
      th <- numeric(nu)
      th[nu] <- sqrt(A[nu] / A[nu + 1])           # theta_{k,nu}
      if (nu >= 2) th[seq_len(nu - 1)] <-
        a[nu] * a[seq_len(nu - 1)] / sqrt(A[nu] * A[nu + 1])
      poly[k, ]  <- (-a[nu] / sqrt(A[nu] * A[nu + 1])) * P2
      freqs[[k]] <- (2 * seq_len(nu) - 1) * pi
    } else {                                      # even m -> phi_3, even freqs
      nu <- m %/% 2
      th <- numeric(nu)
      th[nu] <- sqrt(B[nu] / B[nu + 1])
      if (nu >= 2) th[seq_len(nu - 1)] <-
        b[nu] * b[seq_len(nu - 1)] / sqrt(B[nu] * B[nu + 1])
      poly[k, ]  <- (-b[nu] / sqrt(B[nu] * B[nu + 1])) * P3
      freqs[[k]] <- (2 * seq_len(nu)) * pi
    }
    thetas[[k]] <- th
  }
  list(poly = poly, freqs = freqs, thetas = thetas, K = K, mode = mode, r = 2L)
}

## BMODE: r -> mode string, so that the paper's notation B_r drives the code.
BMODE <- c("cos", "t", "t2")                      # BMODE[r+1]
basis_r <- function(r, K) aug_coeffs_basis(K, BMODE[r + 1L])

## basis_phi(co, t): evaluate phi_1..phi_K of a candidate basis at the points t.
## Returns a length(t)-by-K matrix. It reads only 'poly', 'freqs', 'thetas' and
## 'K', so it accepts all three candidates unchanged.
basis_phi <- function(co, t) {
  Phi <- matrix(0, length(t), co$K)
  for (k in seq_len(co$K)) {
    v <- co$poly[k, 1] + co$poly[k, 2] * t + co$poly[k, 3] * t^2
    if (!is.null(co$freqs[[k]]))
      v <- v + sqrt(2) * as.vector(cos(outer(t, co$freqs[[k]])) %*% co$thetas[[k]])
    Phi[, k] <- v
  }
  Phi
}

## basis_int(co, lo, hi): the window integrals \int_lo^hi phi_k for equal-length
## vectors lo, hi. Polynomial part in closed form; cosine part via
## sqrt2 theta_{k,j} (sin(omega_j hi) - sin(omega_j lo))/omega_j. Returns
## length(lo)-by-K. This is the candidate-basis analogue of the cosine formula
## used in beta_spacing, so no numerical quadrature is required.
basis_int <- function(co, lo, hi) {
  D  <- matrix(0, length(lo), co$K)
  d1 <- hi - lo; d2 <- hi^2 - lo^2; d3 <- hi^3 - lo^3
  for (k in seq_len(co$K)) {
    v <- co$poly[k, 1] * d1 + co$poly[k, 2] * d2 / 2 + co$poly[k, 3] * d3 / 3
    if (!is.null(co$freqs[[k]])) {
      om <- co$freqs[[k]]
      S  <- (sin(outer(hi, om)) - sin(outer(lo, om))) %*% (co$thetas[[k]] / om)
      v  <- v + sqrt(2) * as.vector(S)
    }
    D[, k] <- v
  }
  D
}

## beta_spacing_basis(): the spacing estimator on an arbitrary candidate basis.
## Identical to beta_spacing except that the order-statistic coefficients D_kl
## are the window integrals of the candidate's phi_k rather than of the raw
## cosines. It still uses NO density g: the dependence on g is implicit in the
## geometry of the order statistics. Symmetric boundary reflection as above.
beta_spacing_basis <- function(Y, Tp, co, h) {
  M <- length(Tp); ord <- order(Tp); Ts <- Tp[ord]; Yo <- Y[, ord, drop = FALSE]
  Text <- function(r) {                          # reflected order statistic at rank r
    out <- numeric(length(r)); lo <- r <= 0; hi <- r > M; in_ <- !lo & !hi
    out[in_] <- Ts[r[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - r[lo]]
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - r[hi]]
    out
  }
  l <- seq_len(M)
  D <- basis_int(co, Text(l - h), Text(l + h)) / (2 * h)   # M x K
  Yo %*% D
}

## ise_curve_basis(): ISE of a candidate-basis fit as a function of K = 1..Kup,
## given the precomputed matrix Phi_t = basis_phi(co, tg). Analogue of
## ise_curve; it walks K by rank-one updates, so the whole curve costs one pass.
ise_curve_basis <- function(beta, Phi_t, MUtrue, Kup = Kmax) {
  out <- numeric(Kup)
  Fh  <- matrix(0, nrow(MUtrue), ncol(MUtrue))
  for (K in seq_len(Kup)) {
    Fh     <- Fh + outer(beta[, K], Phi_t[, K])
    out[K] <- mean((Fh - MUtrue)^2)
  }
  out
}

## sse_curve_basis(): held-out sum of squares of the CV criterion, i.e. the
## discretised \int_U (Y_ij(u) - mu^{(-v)}(u,T_ij))^2 du summed over the
## held-out visits (see Remark rem:cv_integral: the factor 1/N_u is common to
## all candidates and is dropped). Same rank-one walk over K.
sse_curve_basis <- function(beta, Phi_te, Y_te, Kup = Kmax) {
  out <- numeric(Kup)
  Rk  <- Y_te
  for (K in seq_len(Kup)) {
    Rk     <- Rk - outer(beta[, K], Phi_te[, K])
    out[K] <- sum(Rk * Rk)
  }
  out
}

## ---- names used by E12-E13 for the augmented basis B_2 -------------------
## B_2 IS the J = 1 augmented basis, so there is one implementation and these
## are aliases kept for readability of the augmentation experiments.
aug_coeffs       <- function(K) aug_coeffs_basis(K, "t2")
aug_phi          <- basis_phi
aug_integral     <- basis_int
beta_spacing_aug <- beta_spacing_basis
ise_curve_aug    <- ise_curve_basis

## ---- SELF-TESTS ---------------------------------------------------------
#  Run once, at the top of E14. They check the properties the paper claims for
#  the candidates rather than assuming them:
#    1. each candidate is orthonormal on [0,1]  (Simpson quadrature);
#    2. the closed forms agree with a brute-force Gram-Schmidt orthonormalisation
#       of the dictionary -- i.e. against the DEFINITION, not against another
#       implementation of the same algebra;
#    3. every element of B_1 is either a cosine or an odd-index element of B_2
#       -- the fact quoted in Section sec:simu_cvbasis to carry Lemma
#       uniform_bound over to the middle candidate;
#    4. the window integrals agree with numerical quadrature.
# -------------------------------------------------------------------------

check_cvbasis <- function(K = 24L, ng = 20001L, tol = 1e-8, verbose = TRUE) {
  stopifnot(ng %% 2L == 1L)                              # Simpson needs odd ng
  tq <- seq(0, 1, length.out = ng)                       # Simpson nodes
  wq <- simpson_w(ng)
  say <- function(...) if (verbose) cat(sprintf(...))
  
  ## 1. orthonormality
  for (r in 0:2) {
    co  <- basis_r(r, K)
    Phi <- basis_phi(co, tq)
    Gm  <- crossprod(Phi * wq, Phi)
    e   <- max(abs(Gm - diag(K)))
    say("  [B_%d] max|Gram - I| = %.2e   sup_k ||phi_k||_inf = %.4f  (2sqrt2 = %.4f)\n",
        r, e, max(abs(Phi)), 2 * sqrt(2))
    stopifnot(e < 1e-6, max(abs(Phi)) < 4 * sqrt(2))
  }
  
  ## 2. closed forms against a brute-force Gram-Schmidt of the dictionary
  gs_ref <- function(Kg, r) {
    ip   <- function(x, y) sum(wq * x * y)
    dict <- list(rep(1, ng))
    if (r >= 1) dict <- c(dict, list(tq))
    if (r >= 2) dict <- c(dict, list(tq^2))
    m <- 1L
    while (length(dict) < Kg) { dict <- c(dict, list(sqrt(2) * cos(m * pi * tq))); m <- m + 1L }
    Q <- matrix(0, ng, Kg)
    for (k in seq_len(Kg)) {
      v <- dict[[k]]
      for (jj in seq_len(k - 1)) v <- v - ip(v, Q[, jj]) * Q[, jj]
      Q[, k] <- v / sqrt(ip(v, v))
    }
    Q
  }
  for (r in 0:2) {
    Pc <- basis_phi(basis_r(r, 16L), tq)
    Qr <- gs_ref(16L, r)
    Qr <- Qr * rep(sign(colSums(Pc * Qr)), each = ng)    # fix the Gram-Schmidt sign
    e  <- max(abs(Pc - Qr))
    say("  [B_%d] closed form vs brute-force Gram-Schmidt = %.2e\n", r, e)
    stopifnot(e < 1e-9)
  }
  
  ## 3. the B_1 elements are cosines (m even) or B_2 elements (m odd)
  P1 <- basis_phi(basis_r(1L, K), tq)
  P2 <- basis_phi(basis_r(2L, K), tq)
  dev <- 0
  for (m in seq_len(K - 2)) {
    if (m %% 2 == 0) {
      ref <- sqrt(2) * cos(m * pi * tq)          # symmetric block: a cosine
    } else if (3 + m <= K) {
      ref <- P2[, 3 + m]                         # antisymmetric: the B_2 element
    } else next
    dev <- max(dev, max(abs(P1[, 2 + m] - ref)))
  }
  say("  [B_1] max deviation from {cosine, B_2 odd-index} = %.2e\n", dev)
  stopifnot(dev < tol)
  
  ## 4. window integrals against quadrature
  co <- basis_r(1L, K)
  lo <- c(0.03, 0.21, 0.44, 0.66); hi <- lo + c(0.11, 0.19, 0.28, 0.33)
  Dq <- t(sapply(seq_along(lo), function(i) {
    tt <- seq(lo[i], hi[i], length.out = 4001)
    ww <- simpson_w(4001L, lo[i], hi[i])
    as.vector(ww %*% basis_phi(co, tt))
  }))
  e <- max(abs(basis_int(co, lo, hi) - Dq))
  say("  [B_1] max|window integral - quadrature| = %.2e\n", e)
  stopifnot(e < 1e-9)
  invisible(TRUE)
}


# =========================================================================
#  REGULARITY CORE
#
#  Nothing here is a theoretical constant and nothing here is a bound: the
#  only quantities computed are (i) the true coefficients beta_k(u) of a mean
#  in a candidate basis, needed by the unbiasedness diagnostic and by the
#  ORACLE PROJECTION mu*_K (an analytical device, not an estimator), and
#  (ii) the decay exponent of ||beta_k||^2, which fixes the EFFECTIVE
#  REGULARITY cap s* of Assumption decay -- a joint property of the mean and
#  of the basis, and the exponent the rate experiments are read against.
# =========================================================================

ORC_TQ <- seq(0, 1, length.out = ORC_NG); ORC_WQ <- simpson_w(ORC_NG)
ORC_UQ <- seq(0, 1, length.out = ORC_NU); ORC_WU <- simpson_w(ORC_NU)

## orc_beta_norm2(): ||beta_k||^2_{L2(U)} for k = 1..K, in a candidate basis.
## 'Phi' may be supplied when several means share one basis, the evaluation of
## phi_1..phi_K on the grid being the expensive part.
orc_beta_norm2 <- function(co, muf, Phi = basis_phi(co, ORC_TQ)) {
  Bc <- muf(ORC_UQ, ORC_TQ) %*% (Phi * ORC_WQ)                     # nu x K
  as.vector(ORC_WU %*% Bc^2)
}

## orc_decay(): fitted decay exponent p of ||beta_k||^2 ~ c k^{-p}, per parity
## (a polynomial augmentation can accelerate one parity subsequence only), and
## the implied regularity cap s* = (p-1)/2 driven by the SLOWEST subsequence.
## The regression uses the hi - lo + 1 coefficients k = lo,...,hi, split in two
## parity subsequences of about half that size each; lo is taken large enough
## for the asymptotic regime to have set in, hi as large as the quadrature
## resolves accurately.
orc_decay <- function(bn, lo = ORC_KFIT, hi = length(bn)) {
  k <- lo:hi; y <- log(pmax(bn[k], 1e-300)); x <- log(k)
  f <- function(id) -unname(coef(lm(y[id] ~ x[id]))[2])
  pe <- f(seq(1, length(k), by = 2)); po <- f(seq(2, length(k), by = 2))
  p  <- min(pe, po)
  list(p_even = pe, p_odd = po, p = p, s_cap = (p - 1) / 2,
       lo = lo, hi = hi, n_fit = length(k))
}

## beta_true_basis(): the TRUE coefficients beta_k(u) in a candidate basis, by
## Simpson quadrature. Used by the unbiasedness diagnostic and by the oracle
## PROJECTION mu*_K(u,t) = sum_{k<=K} beta_k(u) phi_k(t).
beta_true_basis <- function(u, co, K, muf = MU_FUN) {
  muf(u, ORC_TQ) %*% (basis_phi(co, ORC_TQ)[, 1:K, drop = FALSE] * ORC_WQ)
}
beta_true <- function(u, K) beta_true_basis(u, basis_r(0L, K), K, MU_FUN)


# =========================================================================
#  OUTPUT HELPERS  --  figures (TikZ/PDF), tables (booktabs), formatting
# =========================================================================

## Plotting dictionaries: point symbol, line type and LaTeX label of each
## estimator (keys: orac, unif, bal, MC, spac, plug, ll, spl).
EST_PCH <- c(orac = 8,  unif = 19, bal = 17, MC = 15, spac = 18, plug = 4, ll = 1, spl = 6)
EST_LTY <- c(orac = 3,  unif = 1,  bal = 2,  MC = 3,  spac = 5,  plug = 4, ll = 6, spl = 2)
EST_TEX <- c(orac = "$\\widehat\\mu^{(*)}$",
             unif = "$\\widehat\\mu^{(1)}$",
             bal  = "$\\widehat\\mu^{(2)}$",
             MC   = "$\\widehat\\mu^{(\\mathrm{MC})}$",
             spac = "$\\widehat\\mu^{(\\mathrm{sp})}$",
             plug = "$\\widehat\\mu^{(1)}_{\\widehat g}$",
             ll   = "$\\widehat{\\mu}^{(\\mathrm{LL})}$",
             spl  = "$\\widehat{\\mu}^{(\\mathrm{spl})}$")
PT_CEX <- 1.05; LN_LWD <- 1.5; FIT_LWD <- 0.9

## setup_panel(): open one axis system with light margins and a faint grid.
setup_panel <- function(xlab, ylab, xlim, ylim) {
  par(mar = c(4.0, 4.4, 1.0, 0.9), mgp = c(2.5, 0.8, 0), tcl = -0.3)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       las = 1, bty = "l", cex.axis = 0.9, cex.lab = 1.0)
  faint_grid()
}
faint_grid <- function() grid(col = "grey85", lty = 3, lwd = 0.6)

## pts_mark(): a larger marker used to flag the oracle minimiser on a curve.
pts_mark <- function(x, y, pch) points(x, y, pch = pch, cex = PT_CEX * 1.45, lwd = 1.4)

## legend_box(): a compact, frameless legend.
legend_box <- function(pos, labels, pch, lty, ...)
  legend(pos, legend = labels, pch = pch, lty = lty, lwd = LN_LWD,
         bty = "n", cex = 0.82, seg.len = 2.8, pt.cex = PT_CEX, ...)

## emit_fig(): draw the figure 'plotfun' into OUTDIR/<name>.tex as TikZ code,
## or -- if TikZ is unavailable -- into <name>.pdf plus a one-line stub
## <name>.tex that \includegraphics's the PDF, so the LaTeX document compiles
## either way. With WRITE_PDF = TRUE a PDF copy is always written too.
emit_fig <- function(name, width, height, plotfun) {
  if (USE_TIKZ) {
    tikzDevice::tikz(file.path(OUTDIR, paste0(name, ".tex")),
                     width = width, height = height, standAlone = FALSE)
    plotfun(); dev.off()
    if (WRITE_PDF) {
      pdf(file.path(OUTDIR, paste0(name, ".pdf")), width = width, height = height)
      plotfun(); dev.off()
    }
  } else {
    pdf(file.path(OUTDIR, paste0(name, ".pdf")), width = width, height = height)
    plotfun(); dev.off()
    writeLines(sprintf("\\includegraphics[width=.92\\linewidth]{%s.pdf}", name),
               file.path(OUTDIR, paste0(name, ".tex")))
  }
  cat(sprintf("  wrote %s.tex (%s)\n", name, if (USE_TIKZ) "TikZ" else "PDF stub"))
}

## save_csv(): raw numbers behind a table/figure (only if WRITE_CSV = TRUE).
save_csv <- function(name, df) if (WRITE_CSV) {
  write.csv(df, file.path(OUTDIR, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("  wrote %s.csv\n", name))
}

## latex_table(): one booktabs table float per file. 'headers' is the header
## row (character vector), 'body' a CHARACTER matrix of formatted cells.
latex_table <- function(name, headers, body, caption, label, align, extra = NULL,
                        stretch = 1.15) {
  body <- as.matrix(body)
  lines <- c("\\begin{table}[H]", "\\centering", extra,
             if (!is.null(stretch)) sprintf("\\renewcommand{\\arraystretch}{%s}", stretch) else NULL,
             sprintf("\\begin{tabular}{%s}", align), "\\toprule",
             paste0(paste(headers, collapse = " & "), " \\\\"), "\\midrule",
             paste0(apply(body, 1, paste, collapse = " & "), " \\\\"),
             "\\bottomrule", "\\end{tabular}",
             sprintf("\\caption{%s}", caption),
             sprintf("\\label{%s}", label), "\\end{table}")
  writeLines(lines, file.path(OUTDIR, paste0(name, ".tex")))
  cat(sprintf("  wrote %s.tex\n", name))
}

## add_group_rules(): insert a \midrule before every \multirow row but the
## first, so that a table grouped by n is visually separated. Acts in place on
## OUTDIR/<name>.tex, which need not sit under the working directory.
add_group_rules <- function(name) {
  f <- file.path(OUTDIR, paste0(name, ".tex"))
  if (!file.exists(f)) { warning(sprintf("Could not find %s.tex", name)); return(invisible(FALSE)) }
  z <- readLines(f, warn = FALSE)
  ml <- grep("\\\\multirow\\{", z)
  if (length(ml) > 1L) {
    off <- 0L
    for (pos in ml[-1]) {
      q <- pos + off
      if (q > 1L && !grepl("^\\s*\\\\midrule\\s*$", z[q - 1L])) {
        z <- append(z, "\\midrule", after = q - 1L); off <- off + 1L
      }
    }
    writeLines(z, f)
  }
  invisible(TRUE)
}

## multirow_col(): the leading column of a table grouped by 'by', with a
## \multirow cell on the first row of each group and blanks below it.
multirow_col <- function(by, lab = by) {
  out <- rep("", length(by))
  for (idx in split(seq_along(by), by))
    out[idx[1]] <- sprintf("\\multirow{%d}{*}{%s}", length(idx), lab[idx[1]])
  out
}

## ---- number formatting ---------------------------------------------------
fmt_int <- function(v) sprintf("%d", as.integer(round(v)))
fmt_num <- function(v, d = 3) sprintf(paste0("%.", d, "f"), v)

## common_exp(): the single power of ten a whole table can be expressed in,
## chosen so that every cell is at least one unit: e = floor(log10(min|v|)).
## The exponent is reported in the caption, which keeps the cells narrow.
common_exp <- function(v) {
  v <- v[is.finite(v) & v != 0]
  as.integer(floor(log10(min(abs(v)))))
}

## fmt_scaled(): "value (s.e.)" in units of 10^e, with no per-cell exponent.
fmt_scaled <- function(v, se = NULL, e = common_exp(v), d = 2) {
  f <- 10^e
  if (is.null(se)) sprintf(paste0("$%.", d, "f$"), v / f)
  else sprintf(paste0("$%.", d, "f\\,(%.", d, "f)$"), v / f, se / f)
}

## fmt_sci(): plain scientific notation $m.mm \times 10^{e}$.
fmt_sci <- function(v) {
  out <- character(length(v))
  for (i in seq_along(v)) {
    if (!is.finite(v[i])) { out[i] <- "--"; next }
    if (v[i] == 0)        { out[i] <- "$0$"; next }
    e <- floor(log10(abs(v[i])))
    out[i] <- sprintf("$%.2f\\times 10^{%d}$", v[i] / 10^e, e)
  }
  out
}

## fmt_val_se(): the TABLE CONVENTION of the paper -- "mean (s.e.) in units
## of the common power of ten": $m.mm\,(s.ss)\times 10^{e}$, or plain
## $m.mm\,(s.ss)$ when no exponent is needed. Decimals are widened (up to 4)
## so the s.e. always shows at least one significant digit.
fmt_val_se <- function(v, se, d = 2) {
  out <- character(length(v))
  for (i in seq_along(v)) {
    if (!is.finite(v[i])) { out[i] <- "--"; next }
    av <- abs(v[i])
    if (av == 0) { out[i] <- "$0$"; next }
    if (av >= 0.1 && av < 1000) {                       # no exponent needed
      dd <- d
      if (is.finite(se[i]) && se[i] > 0)
        dd <- min(max(d, ceiling(-log10(se[i]))), 4L)
      out[i] <- sprintf(paste0("$%.", dd, "f\\,(%.", dd, "f)$"), v[i], se[i])
    } else {                                            # common power of ten
      e <- floor(log10(av)); f <- 10^e
      dd <- d
      if (is.finite(se[i]) && se[i] > 0)
        dd <- min(max(d, ceiling(-log10(se[i] / f))), 4L)
      out[i] <- sprintf(paste0("$%.", dd, "f\\,(%.", dd, "f)\\times 10^{%d}$"),
                        v[i] / f, se[i] / f, e)
    }
  }
  out
}

## rate_exp(): the fitted EXPONENT of a power law y ~ c x^gamma, obtained by
## least squares of log y on log x, together with its regression standard
## error. Every rate reported in the paper is an exponent of this kind: what
## is estimated is the ORDER OF DECAY of a risk (or of a tuning parameter) in
## the sample size, so "exponent" is used throughout and never "slope".
rate_exp <- function(x, y) {
  s <- suppressWarnings(summary(lm(y ~ x)))$coefficients
  c(exponent = unname(s[2, 1]), se = unname(s[2, 2]))
}
fmt_exp_se <- function(b, s) sprintf("$%.2f\\,(%.2f)$", b, s)

## ratio_se(): mean(x)/mean(y) for PAIRED samples, with the delta-method
## standard error (the covariance term captures the pairing).
ratio_se <- function(x, y) {
  R <- length(x); mx <- mean(x); my <- mean(y)
  va <- (var(x) / my^2 + mx^2 * var(y) / my^4 - 2 * mx * cov(x, y) / my^3) / R
  c(ratio = mx / my, se = sqrt(max(va, 0)))
}


# =========================================================================
#  E1 -- EFFECTIVE REGULARITY: decay exponent of ||beta_k||^2 and the induced
#        cap s* = (p-1)/2, for each mean and each candidate basis
#        (Table tab_decay)
# =========================================================================
exp_decay <- function() {
  cat("\n[E1] effective regularity of each mean in each basis ...\n")
  PHI <- lapply(0:2, function(r) basis_phi(basis_r(r, ORC_KSUM), ORC_TQ))
  dec0 <- orc_decay(orc_beta_norm2(basis_r(0L, ORC_KSUM), MU_FUN, PHI[[1]]))
  
  dtab <- NULL
  for (r in 0:2) {
    cor <- basis_r(r, ORC_KSUM)
    for (nm in names(MEANS)) {
      d <- orc_decay(orc_beta_norm2(cor, MEANS[[nm]]$f, PHI[[r + 1L]]))
      dtab <- rbind(dtab, data.frame(mean = nm, r = r, p_even = d$p_even,
                                     p_odd = d$p_odd, p = d$p, s_cap = d$s_cap))
    }
  }
  dtab <- dtab[order(match(dtab$mean, names(MEANS)), dtab$r), ]
  save_csv("res_decay", dtab)
  
  bd <- character(0)
  for (nm in names(MEANS)) {
    ro <- dtab[dtab$mean == nm, ]
    bd <- c(bd, paste(c(MEANS[[nm]]$tex, fmt_int(MEANS[[nm]]$r_pred),
                        as.vector(rbind(fmt_num(ro$p, 2), fmt_num(ro$s_cap, 2)))),
                      collapse = " & "))
  }
  writeLines(c("\\begin{table}[H]", "\\centering", "\\begin{tabular}{lccccccc}", "\\toprule",
               paste0(" & & \\multicolumn{2}{c}{$\\mathcal{B}_0$} & ",
                      "\\multicolumn{2}{c}{$\\mathcal{B}_1$} & ",
                      "\\multicolumn{2}{c}{$\\mathcal{B}_2$} \\\\"),
               "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
               "mean & $r$ & $p$ & $s^*$ & $p$ & $s^*$ & $p$ & $s^*$ \\\\",
               "\\midrule", paste0(bd, " \\\\"), "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Decay exponent $p$ of $\\|\\beta_k\\|^2_{L^2(\\mathcal U)}$, ",
                              "of order $k^{-p}$ along its slowest parity subsequence, and induced ",
                              "cap $s^*=(p-1)/2$ of Assumption~\\ref{ass:decay}. Each $p$ is a ",
                              "least-squares fit of $\\log\\|\\beta_k\\|^2$ on $\\log k$ over the %d ",
                              "coefficients $k=%d,\\ldots,%d$, about %d per parity subsequence. ",
                              "Column $r$ is the augmentation predicted by ",
                              "Proposition~\\ref{prop:cosine_cap} and Remark~\\ref{rem:aug_t_only}.}"),
                       dec0$n_fit, dec0$lo, dec0$hi, dec0$n_fit %/% 2L),
               "\\label{tab:decay}", "\\end{table}"),
             file.path(OUTDIR, "tab_decay.tex"))
  cat("  wrote tab_decay.tex\n")
  print(dtab, digits = 4)
  invisible(list(decay = dtab))
}


# =========================================================================
#  E2 -- MONTE CARLO DIAGNOSTIC for unbiasedness of the coefficient
#        estimators (Table tab_unbias; p-values reported in the text)
# =========================================================================
exp_unbias <- function() {
  cat("\n[E2] unbiasedness of beta_hat_k(u) on signal-only data ...\n")
  set.seed(BASESEED + 2L)
  ug <- seq(0, 1, length.out = 50); K <- 8L; n <- 30L; R <- cfg$R_unbias
  uidx <- c(10L, 25L, 40L)                    # u ~ 0.18, 0.49, 0.80
  keys <- c("unif", "bal", "MC", "spac")
  A  <- setNames(lapply(keys, function(k) matrix(0, length(uidx), K)), keys)
  A2 <- A                                     # running sums of B and B^2
  for (r in seq_len(R)) {
    d <- gen_data(n, ug, signal_only = TRUE)  # mu only: removes X and noise,
    # leaving just the design randomness
    w <- make_weights(d)
    B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, K),
              bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, K),
              MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, K),
              spac = beta_spacing(d$Y, d$T, K, h_of_M(length(d$T))))
    for (k in keys) {
      A[[k]]  <- A[[k]]  + B[[k]][uidx, ]
      A2[[k]] <- A2[[k]] + B[[k]][uidx, ]^2
    }
  }
  bt <- beta_true(ug, K)[uidx, ]
  d_cells <- length(uidx) * K                  # number of (u,k) grid cells
  df <- R - 1L
  res <- do.call(rbind, lapply(keys, function(k) {
    m   <- A[[k]] / R                                      # MC mean of beta_hat
    se  <- sqrt(pmax(A2[[k]] / R - m^2, 0) * (R / (R - 1)) / R)  # its MC s.e.
    tval <- (m - bt) / se                                  # one-sample t per cell:
    ## for each cell (u,k), t_{u,k} is the standard one-sample Student t-statistic
    ## (df = R-1) of H0: E[beta_hat_k(u)] = beta_k(u). Under exact unbiasedness
    ## each t_{u,k} ~ t_{R-1}. We test each estimator with the two-sided t-test at
    ## its most extreme cell, Bonferroni-corrected over the d cells (log scale so
    ## the tiny spacing p-value does not underflow). The p-values go in the text.
    maxt <- max(abs(tval))
    logp <- min(0, log(d_cells) + log(2) + pt(-maxt, df = df, log.p = TRUE))
    data.frame(est = k,
               mean_abs_bias = mean(abs(m - bt)),          # average over cells
               max_abs_bias  = max(abs(m - bt)),           # worst cell (scale)
               pval          = exp(logp),                  # Bonferroni t-test p
               log10p        = logp / log(10))
  }))
  save_csv("res_unbias", res)
  latex_table("tab_unbias",
              headers = c("estimator",
                          "$\\operatorname{avg}_{u,k}\\bigl|\\widehat{\\mathbb{E}}\\,\\widehat{\\beta}_k(u)-\\beta_k(u)\\bigr|$",
                          "$\\max_{u,k}\\bigl|\\widehat{\\mathbb{E}}\\,\\widehat{\\beta}_k(u)-\\beta_k(u)\\bigr|$"),
              body = cbind(EST_TEX[res$est], fmt_sci(res$mean_abs_bias),
                           fmt_sci(res$max_abs_bias)),
              caption = sprintf(paste0("Diagnostic for the unbiasedness of ",
                                       "the coefficient estimators on signal-only data ",
                                       "($Y_{ij}=\\mu(\\cdot,T_{ij})$), over ",
                                       "$u\\in\\{0.18,\\,0.49,\\,0.80\\}$ and $k\\le%d$ ($d=%d$ cells; $n=%d$, ",
                                       "$R=%d$): mean and maximum absolute bias ",
                                       "$|\\widehat{\\mathbb{E}}\\,\\widehat{\\beta}_k(u)-\\beta_k(u)|$."),
                                K, d_cells, n, R),
              label = "tab:unbias", align = "lcc")
  ## p-values for the text (two-sided one-sample t-test, Bonferroni over d cells)
  fmtp <- function(l10, p) if (l10 <= -16) "<1e-16"
  else if (p >= 1e-3) sprintf("%.3f", p)
  else sprintf("%.1fe%d", 10^(l10 - floor(l10)), floor(l10))
  cat(sprintf("  Bonferroni t-test p-values (d=%d cells, df=%d) -- for the text:\n", d_cells, df))
  for (i in seq_len(nrow(res)))
    cat(sprintf("    %-4s  p = %s\n", res$est[i], fmtp(res$log10p[i], res$pval[i])))
  print(res[, c("est", "mean_abs_bias", "max_abs_bias", "pval")], digits = 3)
  invisible(res)
}

# =========================================================================
#  E3 -- RISK vs TRUNCATION K, incl. the kernel plug-in
#        (Figure fig_riskK, Table tab_riskK)
# =========================================================================
exp_riskK <- function() {
  cat("\n[E3] risk vs truncation K (incl. kernel plug-in) ...\n")
  set.seed(BASESEED + 3L)
  n <- cfg$n_riskK; R <- cfg$R_riskK
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("unif", "bal", "MC", "spac", "plug")
  IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
  for (r in seq_len(R)) {
    d <- gen_data(n, ug); w <- make_weights(d)
    B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
              bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
              MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
              spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))),
              plug = beta_plugin_loo(d$Y, d$T, w$unif, Kmax))
    for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
  }
  mc   <- sapply(keys, function(k) colMeans(IS[[k]]))      # Kmax x 5 matrix
  Kst  <- apply(mc, 2, which.min)
  Imin <- mc[cbind(Kst, seq_along(keys))]
  Ise  <- sapply(seq_along(keys), function(j) sd(IS[[keys[j]]][, Kst[j]]) / sqrt(R))
  res  <- data.frame(est = keys, Kstar = Kst, MISE = Imin, se = Ise)
  save_csv("res_riskK", res)
  latex_table("tab_riskK",
              headers = c("estimator", "$\\widehat K^*$",
                          "$\\overline{\\mathrm{MISE}}$ at $\\widehat K^*$ (s.e.)"),
              body = cbind(EST_TEX[keys], fmt_int(Kst), fmt_val_se(Imin, Ise)),
              caption = sprintf(paste0("Empirical oracle truncation and the risk attained there ",
                                       "($n=%d$, $R=%d$; MC s.e.\\ in parentheses)."),
                                n, R),
              label = "tab:riskK", align = "lcc")
  emit_fig("fig_riskK", 6.4, 5.8, function() {
    Y <- log10(mc)
    setup_panel("truncation $K$", "$\\log_{10}\\overline{\\mathrm{MISE}}$",
                c(1, Kmax), range(Y) + c(-0.05, 0.40))
    for (j in seq_along(keys)) {
      k <- keys[j]
      lines(1:Kmax, Y[, j], lty = EST_LTY[k], lwd = LN_LWD)
      points(1:Kmax, Y[, j], pch = EST_PCH[k], cex = 0.62)
      pts_mark(Kst[j], Y[Kst[j], j], EST_PCH[k])   # flag the empirical oracle
    }
    legend_box("topright", EST_TEX[keys], EST_PCH[keys], EST_LTY[keys])
  })
  print(res, digits = 3)
  invisible(res)
}


# =========================================================================
#  E4 -- DECAY OF THE DESIGN TERM (signal-only data)
#        (Figure fig_design, Table tab_design)
# =========================================================================
exp_design <- function() {
  cat("\n[E4] decay of the design term (signal-only data) ...\n")
  set.seed(BASESEED + 6L)
  ug <- seq(0, 1, length.out = 40); K <- cfg$K_design; R <- cfg$R_design
  keys <- c("unif", "bal", "spac", "MC")
  res <- NULL
  for (n in cfg$n_design) {
    ## On signal-only data beta_hat_k(u) = M_k(u): its variance ACROSS the
    ## replications, summed over k and integrated in u, IS the design term.
    s1 <- setNames(lapply(keys, function(k) matrix(0, length(ug), K)), keys)
    s2 <- s1; Mv <- numeric(R)                # running sums of B and B^2
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, signal_only = TRUE); Mv[r] <- length(d$T)
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, K),
                bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, K),
                spac = beta_spacing(d$Y, d$T, K, h_of_M(length(d$T))),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, K))
      for (k in keys) { s1[[k]] <- s1[[k]] + B[[k]]; s2[[k]] <- s2[[k]] + B[[k]]^2 }
    }
    row <- list(n = n, M = mean(Mv))
    for (k in keys) {
      V <- (s2[[k]] / R - (s1[[k]] / R)^2) * R / (R - 1)   # unbiased variance
      row[[paste0("D_", k)]] <- sum(colMeans(V))           # sum_k int Var(M_k(u)) du
    }
    res <- rbind(res, as.data.frame(row))
    cat(sprintf("  n=%5d (M~%5.0f) D: %s\n", n, row$M,
                paste(sprintf("%.3g", unlist(row[paste0("D_", keys)])), collapse = " ")))
  }
  save_csv("res_design", res)
  lM <- log10(res$M)
  sl <- sapply(keys, function(k) rate_exp(lM, log10(res[[paste0("D_", k)]])))
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                            fmt_sci(unlist(res[i, paste0("D_", keys)]))),
                          collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rrcccc}", "\\toprule",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys], collapse = " & "), " \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{fitted exponent (s.e.)} & ",
                      paste(fmt_exp_se(sl["exponent", ], sl["se", ]), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Design term $\\mathcal{R}^{\\mathrm{design}}=",
                              "\\sum_{k\\le%d}\\int_{\\mathcal{U}}\\mathrm{Var}(\\Lambda_k(u))\\,\\mathrm{d}u$, ",
                              "isolated from signal-only data, versus $M$ ($R=%d$ per design ",
                              "point). Bottom row: fitted exponent in $M$, regression s.e.}"),
                       K, R),
               "\\label{tab:design}", "\\end{table}"),
             file.path(OUTDIR, "tab_design.tex"))
  cat("  wrote tab_design.tex\n")
  emit_fig("fig_design", 6.4, 4.7, function() {
    x  <- log10(res$M)
    ys <- lapply(keys, function(k) log10(res[[paste0("D_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\mathcal{R}^{\\mathrm{design}}$",
                range(x), range(unlist(ys)) + c(-0.10, 0.55))
    for (j in seq_along(keys)) {
      k <- keys[j]
      abline(lm(ys[[j]] ~ x), lty = EST_LTY[k], lwd = FIT_LWD, col = "grey60")
      lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD)
      points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX)
    }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys], sl["exponent", keys]),
               EST_PCH[keys], EST_LTY[keys])
  })
  invisible(list(res = res, exponents = sl))
}


# =========================================================================
#  E5 -- THE FEASIBLE TRUNCATION RULE           (Table tab_trunc)
# =========================================================================
#  Can a practitioner, who sees one sample and no true mu, reach the
#  empirical oracle Khat*? The subject-level V-fold criterion answers it:
#  tab_trunc puts Khat* and K_CV side by side with the risk inflation
#  IR_CV = ISE(K_CV)/ISE(Khat*) that the feasible rule costs.
# =========================================================================
exp_cv <- function() {
  cat("\n[E5] the feasible truncation rule ...\n")
  V <- 5L; des <- DES
  set.seed(BASESEED + 8L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("unif", "bal", "MC", "spac")
  res <- NULL
  for (n in cfg$n_cv) {
    R   <- cfg$R_cv
    IS  <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    KCV <- setNames(lapply(keys, function(k) integer(R)), keys)
    Mv  <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, des = des); Mv[r] <- length(d$T)
      w <- make_weights(d, des)
      ## Full-sample ISE curves, to translate K_CV into a realised risk
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, des$g, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  des$g, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   des$g, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      ## CV criterion, accumulated fold by fold
      fold <- sample(rep_len(1:V, n))         # random subject-level folds
      sse  <- setNames(lapply(keys, function(k) numeric(Kmax)), keys)
      for (v in 1:V) {
        tr  <- fold[d$subj] != v              # TRUE = training observation
        Ttr <- d$T[tr];  Ytr <- d$Y[, tr,  drop = FALSE]
        Tte <- d$T[!tr]; Yte <- d$Y[, !tr, drop = FALSE]
        ntr  <- sum(fold != v)                # number of training subjects
        nmtr <- table(d$subj[tr])[as.character(d$subj[tr])]
        wtr <- list(unif = rep(1 / length(Ttr), length(Ttr)),
                    bal  = 1 / (ntr * as.numeric(nmtr)),
                    MC   = cn_weights(Ttr, des$G))
        Btr <- list(unif = beta_weighted(Ytr, Ttr, wtr$unif, des$g, Kmax),
                    bal  = beta_weighted(Ytr, Ttr, wtr$bal,  des$g, Kmax),
                    MC   = beta_weighted(Ytr, Ttr, wtr$MC,   des$g, Kmax),
                    spac = beta_spacing(Ytr, Ttr, Kmax, h_of_M(length(Ttr))))
        Pte <- phi_mat(Tte, Kmax)
        for (k in keys) sse[[k]] <- sse[[k]] + sse_curve_basis(Btr[[k]], Pte, Yte, Kmax)
      }
      for (k in keys) KCV[[k]][r] <- which.min(sse[[k]])
    }
    for (k in keys) {
      mc   <- colMeans(IS[[k]])
      Ks   <- which.min(mc)                   # empirical oracle Khat*
      rcv  <- IS[[k]][cbind(seq_len(R), KCV[[k]])]   # ISE realised at K_CV
      infl <- rcv / IS[[k]][, Ks]
      res <- rbind(res, data.frame(
        n = n, M = mean(Mv), est = k, Kstar = Ks,
        K_med = unname(quantile(KCV[[k]], 0.50)),
        K_q1  = unname(quantile(KCV[[k]], 0.25)),
        K_q3  = unname(quantile(KCV[[k]], 0.75)),
        R_cv = mean(rcv), se_cv = sd(rcv) / sqrt(R),
        infl_med = median(infl),
        infl_q90 = unname(quantile(infl, 0.90))))
    }
    cat(sprintf("  n=%5d (M~%6.0f) | med IR_CV: %s\n", n, mean(Mv),
                paste(sprintf("%s %.3f", keys,
                              sapply(keys, function(k)
                                res$infl_med[res$n == n & res$est == k])),
                      collapse = "  ")))
  }
  save_csv("res_cv", res)
  
  ## ---- tab_trunc: the two truncations, and the cost of the feasible one ---
  sel <- res[res$n %in% range(cfg$n_cv), ]
  latex_table("tab_trunc",
              headers = c("$n$", "estimator", "$\\widehat K^*$",
                          "$\\mathrm{med}\\,K_{\\mathrm{CV}}\\;[\\mathrm{IQR}]$",
                          "med.\\ $\\mathrm{IR}_{\\mathrm{CV}}$",
                          "$q_{0.9}(\\mathrm{IR}_{\\mathrm{CV}})$"),
              body = cbind(multirow_col(sel$n, fmt_int(sel$n)),
                           EST_TEX[sel$est],
                           fmt_int(sel$Kstar),
                           sprintf("$%g\\;[%g,\\,%g]$", sel$K_med, sel$K_q1, sel$K_q3),
                           fmt_num(sel$infl_med, 3), fmt_num(sel$infl_q90, 3)),
              caption = sprintf(paste0(
                "The feasible truncation against the empirical oracle, at the two ends of the ",
                "sample-size grid. $\\widehat K^*$ minimises the Monte Carlo risk, ",
                "$K_{\\mathrm{CV}}$ the subject-level five-fold criterion on a single sample; ",
                "$\\mathrm{IR}_{\\mathrm{CV}}=\\mathrm{ISE}(K_{\\mathrm{CV}})/",
                "\\mathrm{ISE}(\\widehat K^*)$ ($R=%d$ replications)."), cfg$R_cv),
              label = "tab:trunc", align = "clcccc")
  add_group_rules("tab_trunc")
  print(res[, c("n", "est", "Kstar", "K_med", "infl_med", "R_cv")], digits = 3)
  invisible(res)
}


# =========================================================================
#  E6 -- TOTAL RISK vs M: convergence rate, empirical oracle Khat* and the
#        local-linear benchmark
#        (Tables tab_rate, tab_Kstar; Figure fig_rate_total)
# =========================================================================
exp_rate <- function() {
  cat("\n[E6] total risk vs M: rate, empirical oracle, LL benchmark ...\n")
  set.seed(BASESEED)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys  <- c("unif", "bal", "MC", "spac")      # the four series estimators
  keys5 <- c(keys, "ll")                       # ... plus the LL benchmark
  res <- NULL; Kres <- NULL
  ISbn <- list(); LLbn <- list()               # per-n ISE blocks (exponent bootstrap)
  for (n in cfg$n_rate) {
    R  <- cfg$R_rate
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    LL <- array(NA_real_, c(R, length(LL_HGRID), length(LL_HGRID_U)))
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug); Mv[r] <- length(d$T)
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      LL[r, , ] <- ll2_ise(d$Y, d$T, ug, tg, MU)           ## benchmark
    }
    Mbar <- mean(Mv)
    row <- list(n = n, M = Mbar); krow <- list(n = n, M = Mbar)
    for (k in keys) {   # EMPIRICAL oracle: minimiser of the MC estimate of the RISK
      mc <- colMeans(IS[[k]]); Ks <- which.min(mc)
      row[[paste0("MISE_", k)]] <- mc[Ks]
      row[[paste0("se_",   k)]] <- sd(IS[[k]][, Ks]) / sqrt(R)
      krow[[paste0("Kstar_", k)]] <- Ks
    }
    mcLL <- apply(LL, c(2, 3), mean)                       # mean ISE per (h_t, h_u)
    iLL  <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]  # joint oracle pair
    row$MISE_ll <- mcLL[iLL[1], iLL[2]]
    row$se_ll   <- sd(LL[, iLL[1], iLL[2]]) / sqrt(R)
    krow$ht_ll  <- LL_HGRID[iLL[1]]                        # oracle h_t (sparse dir.)
    krow$hu_ll  <- LL_HGRID_U[iLL[2]]                      # oracle h_u (0 = slice-wise)
    res  <- rbind(res,  as.data.frame(row))
    Kres <- rbind(Kres, as.data.frame(krow))
    ISbn[[length(ISbn) + 1L]] <- IS        # keep blocks for the exponent bootstrap
    LLbn[[length(LLbn) + 1L]] <- LL
    cat(sprintf("  n=%5d (M~%5.0f) MISE %s | LL %.3g  K* %s  h_t*=%.3f h_u*=%.3f\n",
                n, Mbar,
                paste(sprintf("%.3g", unlist(row[paste0("MISE_", keys)])), collapse = " "),
                row$MISE_ll,
                paste(unlist(krow[paste0("Kstar_", keys)]), collapse = "/"),
                krow$ht_ll, krow$hu_ll))
  }
  save_csv("res_rate_total", res); save_csv("res_Kstar", Kres)
  
  ## Fitted rate EXPONENTS. The POINT estimates come from least squares of the
  ## logged design-point means on log M; the STANDARD ERRORS are obtained by a
  ## nonparametric bootstrap over the R replications (see below), NOT from the
  ## regression. The regression s.e. would only measure scatter of the means
  ## about the line and badly understates the true uncertainty, which is
  ## dominated by Monte Carlo noise in each mean and by the integer/grid
  ## oscillation of the oracle tuning (Khat* jumps in unit steps, ht-hat* on a
  ## grid).
  lM  <- log10(res$M)
  sl  <- sapply(keys5, function(k) rate_exp(lM, log10(res[[paste0("MISE_", k)]])))
  ksl <- sapply(keys,  function(k) rate_exp(lM, log10(Kres[[paste0("Kstar_", k)]])))
  hsl <- rate_exp(lM, log10(Kres$ht_ll))      # prediction -1/5 concerns h_t;
  # h_u has no rate prediction (it is
  # expected to sit at 0, see the tex)
  ## ---- honest bootstrap standard errors of all fitted exponents ----------
  set.seed(BASESEED + 100L)                   # reproducible, after all data gen
  Bb <- 1000L; nN <- length(ISbn)
  xc <- lM - mean(lM); den <- sum(xc^2)       # closed-form LS exponent = (xc.y)/den
  ls_exp <- function(y) sum(xc * y) / den
  nht <- length(LL_HGRID)                     # flatten each LL block to R x (nht*nhu)
  LLflat <- lapply(LLbn, function(A) matrix(A, nrow = dim(A)[1]))
  htidx  <- rep(seq_len(nht), times = length(LL_HGRID_U))   # h_t index of each column
  boot_sl <- matrix(NA_real_, Bb, length(keys5), dimnames = list(NULL, keys5))
  boot_k  <- matrix(NA_real_, Bb, length(keys),  dimnames = list(NULL, keys))
  boot_ht <- numeric(Bb)
  for (b in seq_len(Bb)) {
    yR <- matrix(NA_real_, nN, length(keys5)); yK <- matrix(NA_real_, nN, length(keys))
    yH <- numeric(nN)
    for (i in seq_len(nN)) {
      Ri <- nrow(ISbn[[i]][[1]]); id <- sample.int(Ri, Ri, replace = TRUE)
      for (j in seq_along(keys)) {
        mc <- colMeans(ISbn[[i]][[keys[j]]][id, , drop = FALSE])
        yR[i, j] <- log10(min(mc)); yK[i, j] <- log10(which.min(mc))
      }
      mcLL <- colMeans(LLflat[[i]][id, , drop = FALSE])   # mean ISE per (h_t,h_u) cell
      jmin <- which.min(mcLL)
      yR[i, length(keys5)] <- log10(mcLL[jmin])
      yH[i] <- log10(LL_HGRID[htidx[jmin]])
    }
    for (j in seq_along(keys5)) boot_sl[b, j] <- ls_exp(yR[, j])
    for (j in seq_along(keys))  boot_k[b, j]  <- ls_exp(yK[, j])
    boot_ht[b] <- ls_exp(yH)
  }
  sl["se", ]  <- apply(boot_sl, 2, sd)        # overwrite with bootstrap s.e.
  ksl["se", ] <- apply(boot_k,  2, sd)
  hsl["se"]   <- sd(boot_ht)
  ## the caption should not assert h_u* = 0 unless it holds in this run
  hu_note <- if (all(Kres$hu_ll == 0))
    "The selected $\\widehat h_u^*$ is $0$ at every $n$ and is omitted. " else
      sprintf("The selected $\\widehat h_u^*$ is not always $0$ here (%s). ",
              paste(format(Kres$hu_ll, digits = 3), collapse = ", "))
  
  ## ---- tab_rate: risks with s.e., exponent row at the bottom --------------
  ## All risks are reported in units of 10^{-3} (a single factor for the whole
  ## table, stated in the caption), so each cell is the compact "value (s.e.)"
  ## without a per-cell power of ten -- this keeps the table narrow.
  sc <- 1e3
  fmt_r <- function(v, se) sprintf("$%.2f\\,(%.2f)$", v * sc, se * sc)
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                            fmt_r(unlist(res[i, paste0("MISE_", keys5)]),
                                  unlist(res[i, paste0("se_",   keys5)]))),
                          collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rrccccc}", "\\toprule",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys5], collapse = " & "), " \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{fitted exponent (s.e.)} & ",
                      paste(fmt_exp_se(sl["exponent", ], sl["se", ]), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Total integrated risk at the empirical oracle tuning ",
                              "versus the pooled sample size ($R=%d$ replications per design point). ",
                              "Risks in units of $10^{-3}$, MC s.e.\\ in parentheses. Bottom row: ",
                              "fitted exponent of $\\overline{\\mathrm{MISE}}$ as a power of $M$, with ",
                              "bootstrap s.e.}"),
                       cfg$R_rate),
               "\\label{tab:rate}", "\\end{table}"),
             file.path(OUTDIR, "tab_rate.tex"))
  cat("  wrote tab_rate.tex\n")
  
  ## ---- tab_Kstar: empirical oracle truncation and bandwidth ---------------
  kbody <- character(0)
  for (i in seq_len(nrow(Kres)))
    kbody <- c(kbody, paste(c(fmt_int(Kres$n[i]), fmt_int(Kres$M[i]),
                              fmt_int(unlist(Kres[i, paste0("Kstar_", keys)])),
                              fmt_num(Kres$ht_ll[i], 3)), collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rrccccc}", "\\toprule",
               paste0(" & & \\multicolumn{4}{c}{empirical oracle $\\widehat K^*$} & ",
                      "LL \\\\"),
               "\\cmidrule(lr){3-6}\\cmidrule(lr){7-7}",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys], collapse = " & "),
                      " & $\\widehat h_t^*$ \\\\"),
               "\\midrule", paste0(kbody, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{fitted exponent (s.e.)} & ",
                      paste(fmt_exp_se(ksl["exponent", ], ksl["se", ]), collapse = " & "),
                      " & ", fmt_exp_se(hsl["exponent"], hsl["se"]), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Empirical oracle truncation and $t$-bandwidth versus $M$ ",
                              "(same runs as Table~\\ref{tab:rate}, $R=%d$). %sBottom row: ",
                              "fitted exponent in $M$, bootstrap s.e.}"),
                       cfg$R_rate, hu_note),
               "\\label{tab:Kstar}", "\\end{table}"),
             file.path(OUTDIR, "tab_Kstar.tex"))
  cat("  wrote tab_Kstar.tex\n")
  
  ## ---- fig_rate_total: 5 log-log curves, fitted exponents in the legend ---
  emit_fig("fig_rate_total", 6.4, 5.8, function() {
    x  <- log10(res$M)
    ys <- lapply(keys5, function(k) log10(res[[paste0("MISE_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{MISE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.45))
    for (j in seq_along(keys5)) {
      k <- keys5[j]
      lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD)
      points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX)
    }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys5], sl["exponent", keys5]),
               EST_PCH[keys5], EST_LTY[keys5])
  })
  
  invisible(list(res = res, Kres = Kres, exponents = sl, kexponents = ksl,
                 hexponent = hsl))
}


# =========================================================================
#  E7 -- CONVERGENCE RATE UNDER THE BOUNDARY-SMOOTH MEAN mu_C
#         Same protocol as E6, the local-linear benchmark included, so that
#         the two rate experiments can be read against one another.
#         (Table tab_rate_reg, Figure fig_rate_reg)
# =========================================================================
exp_rate_reg <- function() {
  cat("\n[E7] convergence rate under the boundary-smooth mean mu_C ...\n")
  set.seed(BASESEED + 10L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN_C(ug, tg)
  keys  <- c("unif", "bal", "MC", "spac")      # the four series estimators
  keys5 <- c(keys, "ll")                       # ... plus the LL benchmark
  res <- NULL; ISbn <- list(); LLbn <- list()
  for (n in cfg$n_rate) {
    R  <- cfg$R_rate
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    LL <- array(NA_real_, c(R, length(LL_HGRID), length(LL_HGRID_U)))
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, muf = MU_FUN_C); Mv[r] <- length(d$T)   # second mean
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      LL[r, , ] <- ll2_ise(d$Y, d$T, ug, tg, MU)           ## benchmark
    }
    Mbar <- mean(Mv); row <- list(n = n, M = Mbar)
    for (k in keys) {
      mc <- colMeans(IS[[k]]); Ks <- which.min(mc)
      row[[paste0("MISE_", k)]]  <- mc[Ks]
      row[[paste0("se_",   k)]]  <- sd(IS[[k]][, Ks]) / sqrt(R)
      row[[paste0("Kstar_", k)]] <- Ks
    }
    mcLL <- apply(LL, c(2, 3), mean)                       # mean ISE per (h_t,h_u)
    iLL  <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]  # joint empirical oracle
    row$MISE_ll <- mcLL[iLL[1], iLL[2]]
    row$se_ll   <- sd(LL[, iLL[1], iLL[2]]) / sqrt(R)
    row$ht_ll   <- LL_HGRID[iLL[1]]
    row$hu_ll   <- LL_HGRID_U[iLL[2]]
    res <- rbind(res, as.data.frame(row))
    ISbn[[length(ISbn) + 1L]] <- IS; LLbn[[length(LLbn) + 1L]] <- LL
    cat(sprintf("  n=%5d (M~%5.0f) MISE %s | LL %.3g  Khat* %s  h_t*=%.3f h_u*=%.3f\n",
                n, Mbar,
                paste(sprintf("%.3g", unlist(row[paste0("MISE_", keys)])), collapse = " "),
                row$MISE_ll,
                paste(unlist(row[paste0("Kstar_", keys)]), collapse = "/"),
                row$ht_ll, row$hu_ll))
  }
  save_csv("res_rate_reg", res)
  
  ## fitted exponents with bootstrap s.e. (same convention as tab_rate)
  lM <- log10(res$M)
  sl <- sapply(keys5, function(k) rate_exp(lM, log10(res[[paste0("MISE_", k)]])))
  set.seed(BASESEED + 110L)
  Bb <- 1000L; nN <- length(ISbn)
  xc <- lM - mean(lM); den <- sum(xc^2); ls_exp <- function(y) sum(xc * y) / den
  LLflat  <- lapply(LLbn, function(A) matrix(A, nrow = dim(A)[1]))
  boot_sl <- matrix(NA_real_, Bb, length(keys5), dimnames = list(NULL, keys5))
  for (b in seq_len(Bb)) {
    yR <- matrix(NA_real_, nN, length(keys5))
    for (i in seq_len(nN)) {
      Ri <- nrow(ISbn[[i]][[1]]); id <- sample.int(Ri, Ri, replace = TRUE)
      for (j in seq_along(keys))
        yR[i, j] <- log10(min(colMeans(ISbn[[i]][[keys[j]]][id, , drop = FALSE])))
      yR[i, length(keys5)] <- log10(min(colMeans(LLflat[[i]][id, , drop = FALSE])))
    }
    for (j in seq_along(keys5)) boot_sl[b, j] <- ls_exp(yR[, j])
  }
  sl["se", ] <- apply(boot_sl, 2, sd)
  
  ## table (risks in 10^{-3}, exponent row at the bottom)
  sc <- 1e3; fmt_r <- function(v, se) sprintf("$%.2f\\,(%.2f)$", v * sc, se * sc)
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                            fmt_r(unlist(res[i, paste0("MISE_", keys5)]),
                                  unlist(res[i, paste0("se_",   keys5)]))), collapse = " & "))
  hu_note <- if (all(res$hu_ll == 0))
    "The benchmark again selects $\\widehat h_u^*=0$ at every sample size. " else ""
  writeLines(c("\\begin{table}[H]", "\\centering", "\\begin{tabular}{rrccccc}", "\\toprule",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys5], collapse = " & "), " \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{fitted exponent (s.e.)} & ",
                      paste(fmt_exp_se(sl["exponent", ], sl["se", ]), collapse = " & "), " \\\\"),
               paste0("\\multicolumn{2}{r}{predicted exponent} & ",
                      paste(c(rep("$-7/8$", length(keys)), "$-4/5$"), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{As Table~\\ref{tab:rate}, for the boundary-smooth mean ",
                              "$\\mu_C$ ($R=%d$). %sBottom rows: fitted and predicted exponents, the ",
                              "latter from $s^*=7/2$ (Table~\\ref{tab:decay}) for the series ",
                              "estimators and from second-order local-linear theory for the ",
                              "benchmark.}"),
                       cfg$R_rate, hu_note),
               "\\label{tab:rate_reg}", "\\end{table}"),
             file.path(OUTDIR, "tab_rate_reg.tex"))
  cat("  wrote tab_rate_reg.tex\n")
  
  emit_fig("fig_rate_reg", 6.4, 5.0, function() {
    x <- log10(res$M); ys <- lapply(keys5, function(k) log10(res[[paste0("MISE_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{MISE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.45))
    for (j in seq_along(keys5)) { k <- keys5[j]
    lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD); points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX) }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys5], sl["exponent", keys5]),
               EST_PCH[keys5], EST_LTY[keys5])
  })
  invisible(list(res = res, exponents = sl))
}


# =========================================================================
#  E8 -- CASE 1 vs CASE 2 in balanced / unbalanced designs
#        (Figure fig_case12, Table tab_case12) -- PAIRED comparison
# =========================================================================
exp_case12 <- function() {
  cat("\n[E8] Case 1 vs Case 2: balanced vs unbalanced designs ...\n")
  n <- cfg$n_case12; R <- cfg$R_case12
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  one <- function(type, seed) {               # run one design arm
    set.seed(seed)
    I1 <- matrix(NA_real_, R, Kmax); I2 <- I1
    mA <- numeric(R); mH <- numeric(R)
    for (r in seq_len(R)) {
      d <- if (type == "balanced") gen_data(n, ug, m_fixed = 5L)
      else                    gen_data(n, ug, unbalanced = TRUE)
      w <- make_weights(d)
      if (type == "balanced" && r == 1L)      ## internal consistency:
        stopifnot(max(abs(w$unif - w$bal)) < 1e-15)  # 1/M == 1/(n*5) exactly
      B1 <- beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax)
      B2 <- beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax)
      I1[r, ] <- ise_curve(B1, tg, MU); I2[r, ] <- ise_curve(B2, tg, MU)
      mi <- as.numeric(table(d$subj))
      mA[r] <- mean(mi); mH[r] <- 1 / mean(1 / mi)   # arithmetic/harmonic means
    }
    K1 <- which.min(colMeans(I1)); K2 <- which.min(colMeans(I2))
    x1 <- I1[, K1]; x2 <- I2[, K2]            # paired: same replications
    list(v1 = mean(x1), se1 = sd(x1) / sqrt(R),
         v2 = mean(x2), se2 = sd(x2) / sqrt(R),
         ratio = mean(x1) / mean(x2), se_diff = sd(x1 - x2) / sqrt(R),
         mA = mean(mA), mH = mean(mH))
  }
  bal <- one("balanced",   BASESEED + 4L)
  unb <- one("unbalanced", BASESEED + 5L)
  res <- data.frame(design = c("balanced", "unbalanced"),
                    mA = c(bal$mA, unb$mA), mH = c(bal$mH, unb$mH),
                    v1 = c(bal$v1, unb$v1), se1 = c(bal$se1, unb$se1),
                    v2 = c(bal$v2, unb$v2), se2 = c(bal$se2, unb$se2),
                    ratio = c(bal$ratio, unb$ratio),
                    se_diff = c(bal$se_diff, unb$se_diff))
  save_csv("res_case12", res)
  latex_table("tab_case12",
              headers = c("design", "$\\bar m_A$", "$\\bar m_H$",
                          EST_TEX["unif"], EST_TEX["bal"],
                          "$\\widehat\\mu^{(1)}/\\widehat\\mu^{(2)}$",
                          "s.e.\\ of difference"),
              body = cbind(res$design, fmt_num(res$mA, 2), fmt_num(res$mH, 2),
                           fmt_val_se(res$v1, res$se1), fmt_val_se(res$v2, res$se2),
                           fmt_num(res$ratio, 3), fmt_sci(res$se_diff)),
              caption = sprintf(paste0("Integrated risk at the empirical oracle truncation for ",
                                       "$\\widehat{\\mu}^{(1)}$ and $\\widehat{\\mu}^{(2)}$ under a balanced design ",
                                       "($m_i\\equiv5$) and an unbalanced two-point mixture ($15\\%%$ of subjects ",
                                       "with $m_i=20$, the rest $m_i=2$); $n=%d$, $R=%d$ paired replications. ",
                                       "$\\bar m_A$, $\\bar m_H$: arithmetic and harmonic mean visit counts; ",
                                       "MC s.e.\\ in parentheses, last column that of the paired difference."), n, R),
              label = "tab:case12", align = "lcccccc")
  emit_fig("fig_case12", 6.0, 4.7, function() {
    H  <- rbind(c(bal$v1, unb$v1), c(bal$v2, unb$v2))
    SE <- rbind(c(bal$se1, unb$se1), c(bal$se2, unb$se2))
    colnames(H) <- c("balanced", "unbalanced")
    par(mar = c(2.6, 4.4, 1.0, 0.9), mgp = c(3.3, 0.8, 0))
    bp <- barplot(H, beside = TRUE, col = c("grey35", "grey80"),
                  ylab = "$\\overline{\\mathrm{MISE}}$ at $\\widehat K^*$",
                  ylim = c(0, 1.30 * max(H + 2 * SE)), las = 1)
    suppressWarnings(arrows(bp, H - 2 * SE, bp, H + 2 * SE,
                            angle = 90, code = 3, length = 0.04))
    text(colMeans(bp), apply(H + 2 * SE, 2, max) + 0.10 * max(H),
         sprintf("ratio $=%.3f$", c(bal$ratio, unb$ratio)), cex = 0.88)
    legend("topleft", c(EST_TEX["unif"], EST_TEX["bal"]),
           fill = c("grey35", "grey80"), bty = "n", cex = 0.9)
  })
  print(res[, c("design", "mA", "mH", "v1", "v2", "ratio")], digits = 4)
  invisible(list(res = res))
}


# =========================================================================
#  E9 -- PAIRED LOG-RATIO COMPARISON against the spacing estimator, on the
#        two means mu (boundary derivatives non-zero) and mu_C (zero), in
#        the half-cosine basis B_0 AND in the augmented basis B_2
#        Delta_e = log( ISE(mu^e) / ISE(mu^sp) ), replication by replication
#        (Figures fig_boxratio, fig_boxratio_aug; Table tab_boxratio)
# =========================================================================
## Every competing estimator is taken at its OWN empirical oracle tuning
## (Khat* for the series estimators, the bandwidth pair for the local-linear
## benchmark), so the comparison is between estimation principles and not
## between selectors. The oracle PROJECTION is the exception, and must be:
## its error decreases in K, so it has no oracle truncation of its own, and it
## is read at the truncation of the estimator it is compared with -- at which
## point Delta_star is the log of the approximation share of that estimator's
## risk. Within a replication, ONE draw of the design, of the subject
## processes and of the noise serves the two means (gen_data with MU_ZERO,
## each mean added afterwards) and the two bases, so every log-ratio reported
## here is paired -- across estimators, across means and across bases. The
## local-linear benchmark is basis-free and is fitted once.
exp_boxratio <- function(rbas = c(0L, 2L)) {
  cat("\n[E9] paired log-ratio boxplots against the spacing estimator ...\n")
  set.seed(BASESEED + 13L)
  ug  <- seq(0, 1, length.out = cfg$Nu_eval)
  tg  <- seq(0, 1, length.out = cfg$Nt_eval)
  CO   <- lapply(rbas, function(r) basis_r(r, Kmax))
  PHIt <- lapply(CO, function(co) basis_phi(co, tg))
  ser  <- c("unif", "bal", "MC", "spac")        # the four series estimators
  keys <- c("unif", "bal", "MC", "ll", "orac")  # everything compared to spac
  MS <- MEANS[c("mu", "muC")]                   # the two means compared here
  ns <- cfg$n_box
  nn <- length(ns); nm <- length(MS); ne <- length(keys); nb <- length(rbas)
  MUev <- lapply(MS, function(z) z$f(ug, tg))
  ## error curve of the ORACLE PROJECTION in each basis: true coefficients,
  ## hence truncation bias only.
  ISor <- lapply(seq_len(nb), function(bq) lapply(seq_len(nm), function(q)
    ise_curve_basis(beta_true_basis(ug, CO[[bq]], Kmax, MS[[q]]$f), PHIt[[bq]], MUev[[q]])))
  
  D  <- array(list(), c(nb, nm, ne, nn))        # paired log-ratios
  Mb <- numeric(nn)
  for (j in seq_len(nn)) {
    n <- ns[j]; R <- cfg$R_box
    IS <- array(NA_real_, c(nb, nm, length(ser), R, Kmax))
    LL <- array(NA_real_, c(nm, R, length(LL_HGRID), length(LL_HGRID_U)))
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      po <- gen_data(n, ug, muf = MU_ZERO)      # one pool for the two means
      Tv <- po$T; Mv[r] <- length(Tv); h <- h_of_M(length(Tv))
      w  <- make_weights(po)
      for (q in seq_len(nm)) {
        Y <- po$Y + MS[[q]]$f(ug, Tv)           # same noise, another mean
        LL[q, r, , ] <- ll2_ise(Y, Tv, ug, tg, MUev[[q]])   ## basis-free
        for (bq in seq_len(nb)) {
          co <- CO[[bq]]
          B <- list(unif = beta_weighted_basis(Y, Tv, w$unif, G_DENSITY, co),
                    bal  = beta_weighted_basis(Y, Tv, w$bal,  G_DENSITY, co),
                    MC   = beta_weighted_basis(Y, Tv, w$MC,   G_DENSITY, co),
                    spac = beta_spacing_basis(Y, Tv, co, h))
          for (e in seq_along(ser))
            IS[bq, q, e, r, ] <- ise_curve_basis(B[[ser[e]]], PHIt[[bq]], MUev[[q]])
        }
      }
      if (r %% 25 == 0 || r == 1)
        cat(sprintf("    n=%5d  rep %3d / %3d\n", n, r, R))
    }
    Mb[j] <- mean(Mv)
    isp <- which(ser == "spac")
    for (bq in seq_len(nb)) for (q in seq_len(nm)) {
      xs <- list()
      for (e in seq_along(ser)) {               # empirical oracle Khat* per scheme
        mc <- colMeans(IS[bq, q, e, , ])
        xs[[ser[e]]] <- IS[bq, q, e, , which.min(mc)]
      }
      mcLL    <- apply(LL[q, , , ], c(2, 3), mean)         # mean ISE per (h_t,h_u)
      iL      <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]
      xs$ll   <- LL[q, , iL[1], iL[2]]
      Ksp     <- which.min(colMeans(IS[bq, q, isp, , ]))
      xs$orac <- rep(ISor[[bq]][[q]][Ksp], R)   # SAME truncation, no estimation
      for (e in seq_len(ne)) D[[bq, q, e, j]] <- log(xs[[keys[e]]] / xs$spac)
    }
    cat(sprintf("  n=%5d (M~%5.0f) median Delta by basis x mean x competitor:\n", n, Mb[j]))
    for (bq in seq_len(nb)) for (q in seq_len(nm))
      cat(sprintf("     B_%d %-4s %s\n", rbas[bq], names(MS)[q],
                  paste(sprintf("%s %+.2f", keys,
                                sapply(seq_len(ne), function(e) median(D[[bq, q, e, j]]))),
                        collapse = "  ")))
  }
  
  ## ---- long-format summary ------------------------------------------------
  res <- NULL
  for (bq in seq_len(nb)) for (q in seq_len(nm)) for (e in seq_len(ne)) for (j in seq_len(nn)) {
    v <- D[[bq, q, e, j]]
    res <- rbind(res, data.frame(
      r = rbas[bq], mean = names(MS)[q], est = keys[e], n = ns[j], M = Mb[j],
      med = median(v), q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)),
      win = 100 * mean(v > 0), se = sd(v) / sqrt(length(v))))
  }
  rownames(res) <- NULL
  save_csv("res_boxratio", res)
  
  ## ---- tab_boxratio: median log-ratio and win rate, basis B_0 -------------
  r0 <- res[res$r == 0L, ]
  lines <- character(0)
  for (q in seq_len(nm)) {
    if (q > 1L) lines <- c(lines, "\\midrule")
    for (j in seq_len(nn)) {
      cells <- sapply(keys, function(k) {
        z <- r0[r0$mean == names(MS)[q] & r0$est == k & r0$n == ns[j], ]
        if (k == "orac") sprintf("$%.2f$", z$med)
        else sprintf("$%+.2f$ (%.0f\\%%)", z$med, z$win)
      })
      lines <- c(lines, paste(c(if (j == 1L)
        sprintf("\\multirow{%d}{*}{%s}", nn, MS[[q]]$tex) else "",
        fmt_int(Mb[j]), cells), collapse = " & "))
    }
  }
  writeLines(c("\\begin{table}[H]", "\\centering", "\\small",
               sprintf("\\begin{tabular}{ll%s}", strrep("c", ne)), "\\toprule",
               paste0("mean & $\\overline M$ & ", paste(EST_TEX[keys], collapse = " & "),
                      " \\\\"),
               "\\midrule", paste0(lines, " \\\\"), "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Median paired log-ratio $\\Delta_e$ against the spacing ",
                              "estimator over $R=%d$ replications and, in parentheses, the percentage ",
                              "of replications on which the spacing estimator wins ($\\Delta_e>0$); ",
                              "$\\Delta_*<0$ by construction, hence no win rate. All fits use ",
                              "$\\mathcal{B}_0$; Figure~\\ref{fig:boxratio} shows the distributions ",
                              "and Figure~\\ref{fig:boxratio_aug} their counterpart in ",
                              "$\\mathcal{B}_2$.}"),
                       cfg$R_box),
               "\\label{tab:boxratio}", "\\end{table}"),
             file.path(OUTDIR, "tab_boxratio.tex"))
  cat("  wrote tab_boxratio.tex\n")
  
  ## ---- fig_boxratio (B_0) and fig_boxratio_aug (B_2) ----------------------
  shade <- grey.colors(nn, start = 0.88, end = 0.55)
  draw <- function(bq) function() {
    ## a legend strip spanning the width, then one panel per mean
    layout(matrix(c(rep(1, nm), 1 + seq_len(nm)), nrow = 2, byrow = TRUE),
           heights = c(0.14, 1))
    par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", sprintf("$n=%d$", ns), fill = shade, horiz = TRUE,
           bty = "n", cex = 0.9)
    yl <- range(sapply(seq_len(nm), function(q) sapply(seq_len(ne), function(e)
      sapply(seq_len(nn), function(j) quantile(D[[bq, q, e, j]], c(.02, .98))))))
    at <- as.vector(sapply(seq_len(ne), function(e) (e - 1) * (nn + 1) + seq_len(nn)))
    ctr <- sapply(seq_len(ne), function(e) mean((e - 1) * (nn + 1) + seq_len(nn)))
    for (q in seq_len(nm)) {
      par(mar = c(3.4, if (q == 1L) 4.4 else 2.6, 1.8, 0.5),
          mgp = c(3.0, 0.7, 0), tcl = -0.3)
      ## unlist(., recursive = FALSE), NOT as.vector(sapply(.)): the latter
      ## leaves a dim attribute on the list, which boxplot() cannot digest.
      boxplot(unlist(lapply(seq_len(ne), function(e)
        lapply(seq_len(nn), function(j) D[[bq, q, e, j]])), recursive = FALSE),
        at = at, col = rep(shade, ne), boxwex = 0.7, outline = FALSE,
        xaxt = "n", las = 1, bty = "l", cex.axis = 0.85, ylim = yl,
        ylab = if (q == 1L) "$\\Delta_e$" else "")
      abline(h = 0, lty = 2, col = "grey35")
      mtext(EST_TEX[keys], side = 1, at = ctr, line = 0.7, cex = 0.72)
      mtext(MS[[q]]$tex, side = 3, line = 0.4, cex = 0.85)
    }
  }
  for (bq in seq_len(nb))
    emit_fig(if (rbas[bq] == 0L) "fig_boxratio" else "fig_boxratio_aug",
             9.2, 4.1, draw(bq))
  
  print(res[res$est != "orac" & res$r == 0L, c("mean", "est", "n", "med", "win")], digits = 3)
  invisible(list(res = res, D = D, ns = ns, M = Mb, rbas = rbas))
}


# =========================================================================
#  E10 -- ORACLE g vs KERNEL PLUG-IN vs SPACING
#        (Figure fig_plugin, Table tab_plugin)
# =========================================================================
exp_plugin <- function() {
  cat("\n[E10] true g vs kernel plug-in vs spacing ...\n")
  set.seed(BASESEED + 7L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("trueg", "plug", "spac")
  res <- NULL
  for (n in cfg$n_plug) {
    R  <- cfg$R_plug
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug); Mv[r] <- length(d$T)
      wu <- rep(1 / length(d$T), length(d$T))
      B <- list(trueg = beta_weighted(d$Y, d$T, wu, G_DENSITY, Kmax),
                plug  = beta_plugin_loo(d$Y, d$T, wu, Kmax),
                spac  = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
    }
    Ks <- sapply(keys, function(k) which.min(colMeans(IS[[k]])))
    xo <- IS$trueg[, Ks["trueg"]]; xp <- IS$plug[, Ks["plug"]]
    xs <- IS$spac[,  Ks["spac"]]
    rp <- ratio_se(xp, xo); rs <- ratio_se(xs, xo)   # paired delta-method s.e.
    res <- rbind(res, data.frame(n = n, M = mean(Mv),
                                 MISE_trueg = mean(xo), se_trueg = sd(xo) / sqrt(R),
                                 MISE_plug  = mean(xp), se_plug  = sd(xp) / sqrt(R),
                                 MISE_spac  = mean(xs), se_spac  = sd(xs) / sqrt(R),
                                 r_plug = unname(rp["ratio"]), se_rplug = unname(rp["se"]),
                                 r_spac = unname(rs["ratio"]), se_rspac = unname(rs["se"])))
    cat(sprintf("  n=%5d  ratios: plug-in %.3f  spacing %.3f\n",
                n, rp["ratio"], rs["ratio"]))
  }
  save_csv("res_plugin", res)
  latex_table("tab_plugin",
              headers = c("$n$", "$\\overline M$", "true $g$",
                          "plug-in $\\widehat g$ $/$ true $g$", "spacing $/$ true $g$"),
              body = cbind(fmt_int(res$n), fmt_int(res$M),
                           fmt_val_se(res$MISE_trueg, res$se_trueg),
                           fmt_val_se(res$r_plug, res$se_rplug),
                           fmt_val_se(res$r_spac, res$se_rspac)),
              caption = sprintf(paste0("Cost of (not) knowing $%s$: risk of the infeasible ",
                                       "estimator that uses the true $g$, and risk ratios to it of the ",
                                       "kernel plug-in and of the spacing estimator ($R=%d$ per design ",
                                       "point; s.e.\\ of each ratio by the delta method on the paired ",
                                       "replications)."), G_TEX, cfg$R_plug),
              label = "tab:plugin", align = "rrccc")
  emit_fig("fig_plugin", 6.4, 4.7, function() {
    x    <- log10(res$M)
    ys   <- list(trueg = log10(res$MISE_trueg),
                 plug  = log10(res$MISE_plug),
                 spac  = log10(res$MISE_spac))
    pchv <- c(trueg = 1, plug = 4, spac = 18)
    ltyv <- c(trueg = 1, plug = 4, spac = 5)
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{MISE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.35))
    for (k in names(ys)) {
      lines(x, ys[[k]], lty = ltyv[k], lwd = LN_LWD)
      points(x, ys[[k]], pch = pchv[k], cex = PT_CEX)
    }
    legend_box("topright",
               c(paste0(EST_TEX["unif"], " (true $g$)"),
                 EST_TEX["plug"],
                 paste0(EST_TEX["spac"], " (no $g$)")),
               pchv, ltyv, y.intersp = 1.4)
  })
  print(res[, c("n", "M", "MISE_trueg", "r_plug", "r_spac")], digits = 3)
  invisible(res)
}


# =========================================================================
#  E11 -- ROBUSTNESS: heteroscedastic noise / skewed design / high noise
#         (Table tab_robust)
# =========================================================================
exp_robust <- function() {
  cat("\n[E11] robustness: heteroscedastic / skewed design / high noise ...\n")
  set.seed(BASESEED + 9L)
  n <- cfg$n_rob; R <- cfg$R_rob
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("unif", "bal", "MC", "spac")
  scens <- list(
    base  = list(des = DES,      tau_fun = TAU_FUN,                        lab = "baseline"),
    het   = list(des = DES,      tau_fun = function(t) 0.2 + 0.6 * t,      lab = "heterosc."),
    skew  = list(des = DES_SKEW, tau_fun = TAU_FUN,                        lab = "skewed design"),
    hinoi = list(des = DES,      tau_fun = function(t) rep(1, length(t)),  lab = "high noise"))
  res <- NULL
  for (s in names(scens)) {
    sc <- scens[[s]]
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, des = sc$des, tau_fun = sc$tau_fun)
      w <- make_weights(d, sc$des)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, sc$des$g, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  sc$des$g, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   sc$des$g, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
    }
    row <- list(scen = s, lab = sc$lab)
    for (k in keys) {
      mc <- colMeans(IS[[k]]); Ks <- which.min(mc)
      row[[paste0("MISE_", k)]] <- mc[Ks]
      row[[paste0("se_",   k)]] <- sd(IS[[k]][, Ks]) / sqrt(R)
    }
    res <- rbind(res, as.data.frame(row, stringsAsFactors = FALSE))
    cat(sprintf("  %-6s MISE: %s\n", s,
                paste(sprintf("%.3g", unlist(row[paste0("MISE_", keys)])), collapse = " ")))
  }
  rownames(res) <- res$scen
  save_csv("res_robust", res[, setdiff(names(res), "lab")])
  e10 <- common_exp(unlist(res[, paste0("MISE_", keys)]))   # one scale, whole table
  latex_table("tab_robust",
              headers = c("scenario", EST_TEX[keys]),
              body = cbind(res$lab,
                           fmt_scaled(res$MISE_unif, res$se_unif, e10),
                           fmt_scaled(res$MISE_bal,  res$se_bal,  e10),
                           fmt_scaled(res$MISE_MC,   res$se_MC,   e10),
                           fmt_scaled(res$MISE_spac, res$se_spac, e10)),
              caption = sprintf(paste0(
                "Risk at the empirical oracle truncation under four scenarios ($n=%d$, $R=%d$), in ",
                "units of $10^{%d}$ with MC s.e.\\ in parentheses: baseline ($\\tau=1/2$, ",
                "density $g$); heteroscedastic noise ($\\tau(t)=(1+3t)/5$); skewed design ",
                "($g_2(t)=(2+6t)/5$); high noise ($\\tau=1$)."), n, R, e10),
              label = "tab:robust", align = "lcccc")
  print(res[, c("scen", paste0("MISE_", keys))], digits = 3)
  invisible(res)
}


# =========================================================================
#  E12 -- THE AUGMENTED BASIS B_2: functions and sup-norm
#         (Figure fig_augbasis)
# =========================================================================
exp_augbasis <- function() {
  cat("\n[E12] augmented basis B_2: functions and sup-norm (Remark 2 sqrt2) ...\n")
  Kbig <- 50L                                    # sup-norm panel range
  co   <- basis_r(2L, Kbig)
  co1  <- basis_r(1L, Kbig)
  tt   <- seq(0, 1, length.out = 4001)           # fine grid incl. endpoints
  Phi  <- basis_phi(co, tt); Phi1 <- basis_phi(co1, tt)
  supn <- apply(abs(Phi), 2, max)                # ||phi_k||_inf on the grid
  supn1 <- apply(abs(Phi1), 2, max)
  argt <- tt[apply(abs(Phi), 2, which.max)]      # location of the sup
  
  ## ---- Figure fig_augbasis: (left) first phi_k ; (right) sup-norm vs k ----
  emit_fig("fig_augbasis", 7.6, 3.7, function() {
    par(mfrow = c(1, 2))
    ## left panel: the first augmented basis functions
    tp <- seq(0, 1, length.out = 400); Pp <- basis_phi(co, tp)
    setup_panel("$t$", "$\\phi_k(t)$", c(0, 1), range(Pp[, 1:5]) + c(-0.1, 0.4))
    lab5 <- c("$\\phi_1$", "$\\phi_2$", "$\\phi_3$", "$\\phi_4$", "$\\phi_5$")
    lt5  <- c(1, 2, 3, 4, 5)
    for (k in 1:5) lines(tp, Pp[, k], lty = lt5[k], lwd = LN_LWD)
    legend_box("top", lab5, pch = rep(NA, 5), lty = lt5, horiz = TRUE)
    ## right panel: sup-norm vs k, with the 2 sqrt2 reference level
    setup_panel("$k$", "$\\|\\phi_k\\|_\\infty$", c(1, Kbig), c(0.95, 2 * sqrt(2) + 0.2))
    abline(h = 2 * sqrt(2), col = "grey55", lty = 2, lwd = 1.0)
    lines(1:Kbig, supn, lwd = LN_LWD); points(1:Kbig, supn, pch = 20, cex = 0.55)
    text(Kbig, 2 * sqrt(2) + 0.06, "$2\\sqrt2$", pos = 2, cex = 0.8, col = "grey35")
  })
  cat(sprintf("  max_k ||phi_k||_inf (k<=%d): B_2 %.4f, B_1 %.4f (2 sqrt2 = %.4f)\n",
              Kbig, max(supn), max(supn1), 2 * sqrt(2)))
  cat(sprintf("  argmax at an endpoint for all k: %s\n",
              all(argt < 1e-6 | argt > 1 - 1e-6)))
  invisible(list(supn = supn, supn1 = supn1, argt = argt))
}


# =========================================================================
#  E13 -- RATE COMPARISON  spacing(B_0) / spacing(B_2) / local-linear,
#         for the baseline mean mu and for the boundary-smooth mean mu_C
#         (Figure fig_augrate, Tables tab_augrate_mu, tab_augrate_muC)
# =========================================================================
exp_augrate <- function() {
  cat("\n[E13] rate: spacing B_0 vs B_2 vs local-linear (mu and mu_C) ...\n")
  set.seed(BASESEED + 11L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  co0  <- basis_r(0L, Kmax); co2 <- basis_r(2L, Kmax)
  Phi0 <- basis_phi(co0, tg); Phi2 <- basis_phi(co2, tg)
  keys <- c("sp0", "sp2", "ll")
  KEY_TEX <- c(sp0 = "$\\widehat\\mu^{(\\mathrm{sp})}_{\\mathcal{B}_0}$",
               sp2 = "$\\widehat\\mu^{(\\mathrm{sp})}_{\\mathcal{B}_2}$",
               ll  = "$\\widehat{\\mu}^{(\\mathrm{LL})}$")
  KEY_LTY <- c(sp0 = 5, sp2 = 1, ll = 6)
  KEY_PCH <- c(sp0 = 18, sp2 = 15, ll = 1)
  
  means <- list(
    list(muf = MU_FUN,   tag = "mu",
         lab = "$\\mu$ (non-vanishing $\\partial_t\\mu$ at the boundary)",
         pred0 = "$-3/4$", pred2 = "$-7/8$", predll = "$-4/5$"),
    list(muf = MU_FUN_C, tag = "muC",
         lab = "$\\mu_C$ (vanishing $\\partial_t\\mu$ at the boundary)",
         pred0 = "$-7/8$", pred2 = "$-7/8$", predll = "$-4/5$"))
  
  store <- list()
  for (mm in means) {
    MU <- mm$muf(ug, tg)
    res <- NULL; ISb0 <- list(); ISb2 <- list(); LLbn <- list()
    for (n in cfg$n_rate) {
      R  <- cfg$R_rate
      I0 <- matrix(NA_real_, R, Kmax); I2 <- matrix(NA_real_, R, Kmax)
      LL <- array(NA_real_, c(R, length(LL_HGRID), length(LL_HGRID_U)))
      Mv <- numeric(R)
      for (r in seq_len(R)) {
        d <- gen_data(n, ug, muf = mm$muf); Mv[r] <- length(d$T)
        h <- h_of_M(length(d$T))
        I0[r, ] <- ise_curve_basis(beta_spacing_basis(d$Y, d$T, co0, h), Phi0, MU)
        I2[r, ] <- ise_curve_basis(beta_spacing_basis(d$Y, d$T, co2, h), Phi2, MU)
        LL[r, , ] <- ll2_ise(d$Y, d$T, ug, tg, MU)                     # LL benchmark
      }
      Mbar <- mean(Mv)
      mc0 <- colMeans(I0); K0 <- which.min(mc0)
      mc2 <- colMeans(I2); K2 <- which.min(mc2)
      mcLL <- apply(LL, c(2, 3), mean)
      iLL  <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]
      res <- rbind(res, data.frame(n = n, M = Mbar,
                                   MISE_sp0 = mc0[K0], se_sp0 = sd(I0[, K0]) / sqrt(R), Kstar_sp0 = K0,
                                   MISE_sp2 = mc2[K2], se_sp2 = sd(I2[, K2]) / sqrt(R), Kstar_sp2 = K2,
                                   MISE_ll  = mcLL[iLL[1], iLL[2]], se_ll = sd(LL[, iLL[1], iLL[2]]) / sqrt(R),
                                   ht_ll = LL_HGRID[iLL[1]], hu_ll = LL_HGRID_U[iLL[2]]))
      ISb0[[length(ISb0) + 1L]] <- I0
      ISb2[[length(ISb2) + 1L]] <- I2
      LLbn[[length(LLbn) + 1L]] <- LL
      cat(sprintf("  [%s] n=%5d M~%6.0f | B0 K*=%2d %.3e | B2 K*=%2d %.3e | LL %.3e (h_u*=%.2f)\n",
                  mm$tag, n, Mbar, K0, mc0[K0], K2, mc2[K2],
                  mcLL[iLL[1], iLL[2]], LL_HGRID_U[iLL[2]]))
    }
    
    ## ---- fitted rate exponents with nonparametric bootstrap s.e. -----------
    lM  <- log10(res$M)
    xc <- lM - mean(lM); den <- sum(xc^2); ls_exp <- function(y) sum(xc * y) / den
    sl  <- c(sp0 = ls_exp(log10(res$MISE_sp0)),      # centred-LS exponent,
             sp2 = ls_exp(log10(res$MISE_sp2)),      # matches the bootstrap below
             ll  = ls_exp(log10(res$MISE_ll)))
    Bb <- 1000L; nN <- length(ISb0)
    LLflat <- lapply(LLbn, function(A) matrix(A, nrow = dim(A)[1]))
    boot <- matrix(NA_real_, Bb, 3, dimnames = list(NULL, keys))
    for (bb in seq_len(Bb)) {
      y0 <- numeric(nN); y2 <- numeric(nN); yL <- numeric(nN)
      for (i in seq_len(nN)) {
        Ri <- nrow(ISb0[[i]]); id <- sample.int(Ri, Ri, replace = TRUE)
        y0[i] <- log10(min(colMeans(ISb0[[i]][id, , drop = FALSE])))
        y2[i] <- log10(min(colMeans(ISb2[[i]][id, , drop = FALSE])))
        yL[i] <- log10(min(colMeans(LLflat[[i]][id, , drop = FALSE])))
      }
      boot[bb, ] <- c(ls_exp(y0), ls_exp(y2), ls_exp(yL))
    }
    slse <- apply(boot, 2, sd)
    
    ## ---- table tab_augrate_<tag> (risks in 10^-3, exponent + K* rows) ------
    hu_note <- if (all(res$hu_ll == 0))
      "The empirical oracle $u$-bandwidth of the benchmark is $\\widehat h_u^*=0$ throughout, so it reduces to the slice-wise local-linear smoother. " else ""
    sc <- 1e3; fr <- function(v, se) sprintf("$%.2f\\,(%.2f)$", v * sc, se * sc)
    tb <- character(0)
    for (i in seq_len(nrow(res)))
      tb <- c(tb, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                          fr(res$MISE_sp0[i], res$se_sp0[i]), fmt_int(res$Kstar_sp0[i]),
                          fr(res$MISE_sp2[i], res$se_sp2[i]), fmt_int(res$Kstar_sp2[i]),
                          fr(res$MISE_ll[i],  res$se_ll[i])), collapse = " & "))
    writeLines(c("\\begin{table}[H]", "\\centering",
                 "\\begin{tabular}{rrcccccc}", "\\toprule",
                 paste0(" & & \\multicolumn{2}{c}{", KEY_TEX["sp0"], " (half-cosine)} & ",
                        "\\multicolumn{2}{c}{", KEY_TEX["sp2"], " (augmented)} & ",
                        KEY_TEX["ll"], " \\\\"),
                 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-7}",
                 paste0("$n$ & $\\overline M$ & MISE & $\\widehat K^*$ & MISE & $\\widehat K^*$ & MISE \\\\"),
                 "\\midrule", paste0(tb, " \\\\"), "\\midrule",
                 paste0("\\multicolumn{2}{r}{fitted exponent (s.e.)} & ",
                        "\\multicolumn{2}{c}{", fmt_exp_se(sl["sp0"], slse["sp0"]), "} & ",
                        "\\multicolumn{2}{c}{", fmt_exp_se(sl["sp2"], slse["sp2"]), "} & ",
                        fmt_exp_se(sl["ll"], slse["ll"]), " \\\\"),
                 paste0("\\multicolumn{2}{r}{predicted rate} & ",
                        "\\multicolumn{2}{c}{", mm$pred0, "} & ",
                        "\\multicolumn{2}{c}{", mm$pred2, "} & ", mm$predll, " \\\\"),
                 "\\bottomrule", "\\end{tabular}",
                 sprintf(paste0("\\caption{Rate comparison for the mean %s, at empirical oracle ",
                                "tuning ($R=%d$; risks in units of $10^{-3}$, MC s.e.\\ in ",
                                "parentheses). %sBottom rows: fitted exponent in $M$ with bootstrap ",
                                "s.e., and the rate predicted by Proposition~\\ref{prop:cosine_cap} ",
                                "and Corollary~\\ref{cor:aug_rate} for the series estimators, ",
                                "$M^{-4/5}$ for the benchmark.}"),
                         mm$lab, cfg$R_rate, hu_note),
                 sprintf("\\label{tab:augrate_%s}", mm$tag), "\\end{table}"),
               file.path(OUTDIR, sprintf("tab_augrate_%s.tex", mm$tag)))
    cat(sprintf("  wrote tab_augrate_%s.tex\n", mm$tag))
    store[[mm$tag]] <- list(res = res, sl = sl, slse = slse, lab = mm$lab)
  }
  
  ## ---- Figure fig_augrate: two log-log panels (mu | mu_C) ----------------
  emit_fig("fig_augrate", 7.8, 4.0, function() {
    par(mfrow = c(1, 2))
    for (tg_ in c("mu", "muC")) {
      st <- store[[tg_]]; res <- st$res; x <- log10(res$M)
      ys <- list(sp0 = log10(res$MISE_sp0), sp2 = log10(res$MISE_sp2), ll = log10(res$MISE_ll))
      setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{MISE}}$",
                  range(x), range(unlist(ys)) + c(-0.05, 0.55))
      for (k in keys) { lines(x, ys[[k]], lty = KEY_LTY[k], lwd = LN_LWD)
        points(x, ys[[k]], pch = KEY_PCH[k], cex = PT_CEX) }
      title(main = if (tg_ == "mu") "$\\mu$: boundary derivative $\\neq 0$"
            else "$\\mu_C$: boundary derivative $=0$", cex.main = 0.95)
      legend_box("topright", sprintf("%s  ($%.2f$)", KEY_TEX[keys], st$sl[keys]), y.intersp = 1.8,
                 KEY_PCH[keys], KEY_LTY[keys])
    }
  })
  invisible(store)
}


# =========================================================================
#  E14 -- CHOOSING THE BASIS BY CROSS-VALIDATION (Table tab_cvbasis)
# =========================================================================
#  The basis and the truncation enter the fitted surface in the same way, so
#  the selector of E5 needs no modification: it suffices to minimise the
#  subject-level V-fold criterion over the 3 x Kmax PAIRS (basis, truncation)
#  instead of over the truncation alone. Three means whose correct basis is
#  known in advance are used as witnesses:
#
#     mu    boundary derivatives non-zero, DIFFERENT ->  B_2   ({t, t^2})
#     mu_B  boundary derivatives non-zero, EQUAL     ->  B_1   ({t})
#     mu_C  boundary derivatives ZERO                ->  B_0   (half-cosine)
#
#  Per replication and per mean we record the whole ISE curve in each basis
#  (whose Monte Carlo average gives the empirical oracle K* and the reference
#  risk of that basis) and, independently, the basis and truncation selected by
#  CV together with the ISE actually realised there. The reported risk inflation
#  IR_CV = ISE(bhat, K_CV) / ISE(b*, K*_{b*}) compares the fully data-driven fit
#  with the empirical oracle in BOTH arguments. The table reports the SELECTION
#  FREQUENCIES and IR_CV only; the reference risks are printed to the console.
# -------------------------------------------------------------------------
exp_cvbasis <- function(n = cfg$n_bcv, R = cfg$R_bcv, V = 5L, des = DES) {
  cat("\n[E14] choosing the basis by cross-validation ...\n")
  cat("  self-tests:\n"); check_cvbasis()
  set.seed(BASESEED + 12L)
  u_ev <- seq(0, 1, length.out = cfg$Nu_eval)
  t_ev <- seq(0, 1, length.out = cfg$Nt_eval)
  
  blab   <- c("$\\mathcal{B}_0$", "$\\mathcal{B}_1$", "$\\mathcal{B}_2$")
  co     <- lapply(0:2, function(r) basis_r(r, Kmax))
  Phi_t  <- lapply(co, function(x) basis_phi(x, t_ev))
  means  <- list(mu = MU_FUN, muB = MU_FUN_B, muC = MU_FUN_C)
  mlab   <- c("$\\mu$", "$\\mu_B$", "$\\mu_C$")
  rpred  <- c(2L, 1L, 0L)                             # the basis theory predicts
  MU_ev  <- lapply(means, function(f) f(u_ev, t_ev))
  
  nm <- length(means); nb <- 3L
  ISE    <- array(NA_real_, c(nm, nb, R, Kmax))       # ISE curves
  SELb   <- matrix(NA_integer_, nm, R)                # basis chosen by CV
  SELK   <- matrix(NA_integer_, nm, R)                # truncation chosen by CV
  ISEsel <- matrix(NA_real_,    nm, R)                # ISE realised at that pair
  
  for (r in seq_len(R)) {
    ## one pool (design, subject processes, noise) shared by the three means
    po   <- gen_data(n, u_ev, des = des, muf = MU_ZERO)
    Tv   <- po$T; M <- length(Tv); h <- h_of_M(M)
    fold <- sample(rep_len(seq_len(V), n))            # folds are SUBJECT-level
    for (q in seq_len(nm)) {
      Y <- po$Y + means[[q]](u_ev, Tv)                # same noise, another mean
      for (bq in seq_len(nb)) {
        Bq <- beta_spacing_basis(Y, Tv, co[[bq]], h)
        ISE[q, bq, r, ] <- ise_curve_basis(Bq, Phi_t[[bq]], MU_ev[[q]], Kmax)
      }
      cvm <- matrix(0, nb, Kmax)                      # CV(b, K)
      for (v in seq_len(V)) {
        tr  <- fold[po$subj] != v
        Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]
        Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
        htr <- h_of_M(length(Ttr))
        for (bq in seq_len(nb)) {
          Btr <- beta_spacing_basis(Ytr, Ttr, co[[bq]], htr)
          cvm[bq, ] <- cvm[bq, ] +
            sse_curve_basis(Btr, basis_phi(co[[bq]], Tte), Yte, Kmax)
        }
      }
      ij <- which(cvm == min(cvm), arr.ind = TRUE)[1, ]      # joint minimiser
      SELb[q, r] <- ij[1]; SELK[q, r] <- ij[2]
      ISEsel[q, r] <- ISE[q, ij[1], r, ij[2]]
    }
    if (r %% 10 == 0 || r == 1) cat(sprintf("    rep %3d / %3d  (M = %d, h = %d)\n",
                                            r, R, M, h))
  }
  
  ## ---- summaries ---------------------------------------------------------
  mISE <- apply(ISE, c(1, 2, 4), mean)                # nm x nb x Kmax
  Kst  <- matrix(0L, nm, nb); Rst <- Sst <- matrix(0, nm, nb)
  for (q in seq_len(nm)) for (bq in seq_len(nb)) {
    k <- which.min(mISE[q, bq, ]); Kst[q, bq] <- k
    Rst[q, bq] <- mISE[q, bq, k]
    Sst[q, bq] <- sd(ISE[q, bq, , k]) / sqrt(R)
  }
  bstar <- apply(Rst, 1, which.min)                   # oracle-best basis per row
  freq  <- t(sapply(seq_len(nm), function(q) tabulate(SELb[q, ], nb) / R * 100))
  IRq   <- t(sapply(seq_len(nm), function(q)
    quantile(ISEsel[q, ] / ISE[q, bstar[q], , Kst[q, bstar[q]]], c(0.25, 0.5, 0.75))))
  
  cat("\n  reference risk at the oracle K, selection frequency, IR_CV\n")
  for (q in seq_len(nm)) {
    cat(sprintf("   %-4s (theory r = %d):", names(means)[q], rpred[q]))
    for (bq in seq_len(nb))
      cat(sprintf("  B_%d %9.3e (K*=%2d, %3.0f%%)%s",
                  bq - 1L, Rst[q, bq], Kst[q, bq], freq[q, bq],
                  if (bq == bstar[q]) "*" else " "))
    cat(sprintf("   IR_CV median %.3f\n", IRq[q, 2]))
  }
  cat(sprintf("  diagonal recovered: reference risk %s, modal CV choice %s\n",
              if (all(bstar == rpred + 1L)) "YES" else "NO",
              if (all(apply(freq, 1, which.max) == rpred + 1L)) "YES" else "NO"))
  
  ## ---- LaTeX table: selection frequencies and IR_CV only ------------------
  fh  <- file.path(OUTDIR, "tab_cvbasis.tex")
  con <- file(fh, "w")
  wl  <- function(...) writeLines(sprintf(...), con)
  wl("\\begin{table}[H]")
  wl("\\centering\\small\\setlength{\\tabcolsep}{8pt}")
  wl("\\begin{tabular}{lccccc}")
  wl("\\toprule")
  wl("& & \\multicolumn{3}{c}{selection frequency (\\%%)} & \\\\")
  wl("\\cmidrule(lr){3-5}")
  wl("mean & $r$ & %s & %s & %s & $\\mathrm{IR}_{\\mathrm{CV}}$ \\\\",
     blab[1], blab[2], blab[3])
  wl("\\midrule")
  for (q in seq_len(nm)) {
    if (q > 1L) wl("\\midrule")
    cells <- sapply(seq_len(nb), function(bq) {
      s <- sprintf("%.0f", freq[q, bq])
      if (bq == bstar[q]) sprintf("$\\mathbf{%s}^{*}$", s) else sprintf("$%s$", s)
    })
    wl("%s & %d & %s & $%.3f$ \\\\", mlab[q], rpred[q],
       paste(cells, collapse = " & "), IRq[q, 2])
    wl(" & & & & & {\\scriptsize $[%.2f,\\,%.2f]$} \\\\", IRq[q, 1], IRq[q, 3])
  }
  wl("\\bottomrule")
  wl("\\end{tabular}")
  wl("\\caption{Choosing the basis by cross-validation, $n=%d$ ($M\\approx%d$), $R=%d$,", n, M, R)
  wl("\t$V=%d$ folds, $K\\le K_{\\max}=%d$. Column $r$ is the augmentation predicted", V, Kmax)
  wl("\tby the theory, the three middle columns the percentage of replications on")
  wl("\twhich each basis is selected by the subject-level criterion, and a star marks")
  wl("\tthe basis with the smallest Monte Carlo risk of the row. The last column is")
  wl("\tthe median risk inflation $\\mathrm{IR}_{\\mathrm{CV}}=")
  wl("\t\\mathrm{ISE}(\\widehat b,K_{\\mathrm{CV}})/")
  wl("\t\\mathrm{ISE}(b^{\\star},\\widehat K^{*}_{b^{\\star}})$, with its interquartile")
  wl("\trange below.}")
  wl("\\label{tab:cvbasis}")
  wl("\\end{table}")
  close(con)
  cat(sprintf("  wrote tab_cvbasis.tex\n"))
  
  res <- list(n = n, M = M, R = R, V = V, rpred = rpred,
              risk = Rst, se = Sst, Kstar = Kst, freq = freq, IR = IRq,
              bstar = bstar, SELb = SELb, SELK = SELK, ISE = ISE)
  if (WRITE_CSV)
    save_csv("res_cvbasis",
             data.frame(mean = rep(names(means), each = nb),
                        basis = rep(paste0("B", 0:2), nm),
                        r_pred = rep(rpred, each = nb),
                        mise = as.vector(t(Rst)), se = as.vector(t(Sst)),
                        Kstar = as.vector(t(Kst)), sel_pct = as.vector(t(freq))))
  invisible(res)
}


# =========================================================================
#  E15 -- THE FULLY DATA-DRIVEN SPACING ESTIMATOR AGAINST THE TWO EXTERNAL
#         BENCHMARKS (local-linear and penalised spline), one figure per mean
#         (Figures fig_cvbox_mu, fig_cvbox_muB, fig_cvbox_muC)
# =========================================================================
#  Same paired log-ratios as E9,
#        Delta_e = log( ISE(mu^e) / ISE(mu^sp) ),
#  but the spacing estimator no longer receives any oracle: on each
#  replication the subject-level V-fold criterion of E14 selects the PAIR
#  (basis, truncation) over the 3 x Kmax candidates, the estimator is refitted
#  on the whole sample at that pair, and its ISE is the denominator. Two
#  references, BOTH kept at their own empirical oracle so that the comparison
#  is the conservative one:
#     ll    the bivariate local-linear benchmark, at its oracle (h_t, h_u);
#     spl   the penalised cubic B-spline benchmark, at its oracle lambda.
#  The oracle PROJECTION is still computed at (selected basis, K_CV) and its
#  median is printed to the console -- it is the approximation share of the
#  risk, near -log(1+2s*) -- but it does not enter the figures.
#  Within a replication one draw of the design, of the subject processes and
#  of the noise serves the three means (gen_data with MU_ZERO, each mean added
#  afterwards), so the three figures are paired and differ only through mu.
# -------------------------------------------------------------------------
exp_cvbox <- function(ns = cfg$n_cvbox, R = cfg$R_cvbox, V = 5L, des = DES) {
  cat("\n[E15] data-driven spacing (basis + K by CV) vs LL and spline ...\n")
  set.seed(BASESEED + 14L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  co    <- lapply(0:2, function(r) basis_r(r, Kmax))
  Phi_t <- lapply(co, function(x) basis_phi(x, tg))
  MS    <- MEANS                       # mu, mu_B, mu_C: the three regimes
  nm <- length(MS); nb <- 3L; nn <- length(ns); ne <- 2L
  keys <- c("ll", "spl")
  MUev <- lapply(MS, function(z) z$f(ug, tg))
  ## error curve of the ORACLE PROJECTION in every (mean, basis): true
  ## coefficients, hence truncation bias alone. Console diagnostic only.
  ISor <- lapply(seq_len(nm), function(q)
    sapply(seq_len(nb), function(bq)
      ise_curve_basis(beta_true_basis(ug, co[[bq]], Kmax, MS[[q]]$f),
                      Phi_t[[bq]], MUev[[q]], Kmax)))
  D    <- array(list(), c(nm, ne, nn))          # paired log-ratios
  freq <- array(0, c(nm, nn, nb))               # CV basis-selection frequencies
  Mb   <- numeric(nn)
  for (j in seq_len(nn)) {
    n   <- ns[j]
    Isp <- matrix(NA_real_, nm, R)              # ISE of the CV-selected fit
    Ior <- matrix(NA_real_, nm, R)              # oracle projection, SAME pair
    LL  <- array(NA_real_, c(nm, R, length(LL_HGRID), length(LL_HGRID_U)))
    SP  <- array(NA_real_, c(nm, R, length(SPL_LGRID)))
    Mv  <- numeric(R)
    for (r in seq_len(R)) {
      po   <- gen_data(n, ug, des = des, muf = MU_ZERO)   # one pool, three means
      Tv   <- po$T; Mv[r] <- length(Tv); h <- h_of_M(length(Tv))
      fold <- sample(rep_len(seq_len(V), n))              # SUBJECT-level folds
      for (q in seq_len(nm)) {
        Y <- po$Y + MS[[q]]$f(ug, Tv)                     # same noise, another mean
        LL[q, r, , ] <- ll2_ise(Y, Tv, ug, tg, MUev[[q]]) # basis-free benchmarks
        SP[q, r, ]   <- spl_ise(Y, Tv, tg, MUev[[q]])
        cvm <- matrix(0, nb, Kmax)                        # CV(basis, K)
        for (v in seq_len(V)) {
          tr  <- fold[po$subj] != v
          Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]
          Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
          htr <- h_of_M(length(Ttr))
          for (bq in seq_len(nb)) {
            Btr <- beta_spacing_basis(Ytr, Ttr, co[[bq]], htr)
            cvm[bq, ] <- cvm[bq, ] +
              sse_curve_basis(Btr, basis_phi(co[[bq]], Tte), Yte, Kmax)
          }
        }
        ij <- which(cvm == min(cvm), arr.ind = TRUE)[1, ] # joint minimiser
        bs <- ij[1]; Ks <- ij[2]
        freq[q, j, bs] <- freq[q, j, bs] + 1
        Bfull <- beta_spacing_basis(Y, Tv, co[[bs]], h)   # refit, full sample
        Isp[q, r] <- ise_curve_basis(Bfull, Phi_t[[bs]], MUev[[q]], Kmax)[Ks]
        Ior[q, r] <- ISor[[q]][Ks, bs]
      }
      if (r %% 25 == 0 || r == 1)
        cat(sprintf("    n=%5d  rep %3d / %3d\n", n, r, R))
    }
    Mb[j] <- mean(Mv); freq[, j, ] <- freq[, j, ] / R * 100
    for (q in seq_len(nm)) {
      mcLL <- apply(LL[q, , , ], c(2, 3), mean)           # its own empirical oracle
      iL   <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]
      mcSP <- colMeans(SP[q, , ])                         # its own empirical oracle
      iS   <- which.min(mcSP)
      D[[q, 1, j]] <- log(LL[q, , iL[1], iL[2]] / Isp[q, ])
      D[[q, 2, j]] <- log(SP[q, , iS]           / Isp[q, ])
      if (iS %in% c(1L, length(SPL_LGRID)))
        warning(sprintf("spline oracle lambda at a grid edge (mean %s, n=%d)",
                        names(MS)[q], ns[j]))
      cat(sprintf("     %-4s  Delta_LL %+.2f (win %3.0f%%)  Delta_spl %+.2f (win %3.0f%%)",
                  names(MS)[q], median(D[[q, 1, j]]), 100 * mean(D[[q, 1, j]] > 0),
                  median(D[[q, 2, j]]), 100 * mean(D[[q, 2, j]] > 0)))
      cat(sprintf("  |  modal basis B_%d (%.0f%%)  [Delta_* %+.2f]\n",
                  which.max(freq[q, j, ]) - 1L, max(freq[q, j, ]),
                  median(log(Ior[q, ] / Isp[q, ]))))
    }
    cat(sprintf("  n=%5d (M~%5.0f) done\n", ns[j], Mb[j]))
  }
  
  ## ---- long-format summary ------------------------------------------------
  res <- NULL
  for (q in seq_len(nm)) for (e in seq_len(ne)) for (j in seq_len(nn)) {
    v <- D[[q, e, j]]
    res <- rbind(res, data.frame(
      mean = names(MS)[q], est = keys[e], n = ns[j], M = Mb[j],
      med = median(v), q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)),
      win = 100 * mean(v > 0), se = sd(v) / sqrt(length(v)),
      basis_modal = which.max(freq[q, j, ]) - 1L,
      basis_pct   = max(freq[q, j, ])))
  }
  rownames(res) <- NULL
  save_csv("res_cvbox", res)
  
  ## ---- one figure per mean ------------------------------------------------
  ## Two competitor groups (local-linear, penalised spline), one box per sample
  ## size, exactly the layout of fig_boxratio; the modal basis chosen by the
  ## criterion at each sample size is printed above the panel.
  shade <- grey.colors(nn, start = 0.88, end = 0.55)
  for (q in seq_len(nm)) {
    emit_fig(sprintf("fig_cvbox_%s", names(MS)[q]), 6.8, 4.2, function() {
      yl <- range(sapply(seq_len(ne), function(e)
        sapply(seq_len(nn), function(j) quantile(D[[q, e, j]], c(.02, .98)))))
      yl <- yl + c(-0.05, 0.35) * max(diff(yl), 0.1)     # room for the legend
      at  <- as.vector(sapply(seq_len(ne), function(e) (e - 1) * (nn + 1) + seq_len(nn)))
      ctr <- sapply(seq_len(ne), function(e) mean((e - 1) * (nn + 1) + seq_len(nn)))
      par(mar = c(3.4, 4.4, 2.2, 0.6), mgp = c(3.0, 0.7, 0), tcl = -0.3)
      ## unlist(., recursive = FALSE), NOT as.vector(sapply(.)): the latter
      ## leaves a dim attribute on the list, which boxplot() cannot digest.
      boxplot(unlist(lapply(seq_len(ne), function(e)
        lapply(seq_len(nn), function(j) D[[q, e, j]])), recursive = FALSE),
        at = at, col = rep(shade, ne), boxwex = 0.7, outline = FALSE,
        xaxt = "n", las = 1, bty = "l", cex.axis = 0.85, ylim = yl,
        ylab = "$\\Delta_e$")
      abline(h = 0, lty = 2, col = "grey35")
      mtext(EST_TEX[keys], side = 1, at = ctr, line = 0.7, cex = 0.8)
      fq <- matrix(freq[q, , ], nn, nb)
      mtext(paste0(MS[[q]]$tex, " -- CV basis: ",
                   paste(sprintf("$n{=}%d$: $\\mathcal{B}_%d$ (%.0f\\%%)", ns,
                                 apply(fq, 1, which.max) - 1L,
                                 apply(fq, 1, max)), collapse = ", ")),
            side = 3, line = 0.35, cex = 0.68)
      legend("top", sprintf("$n=%d$", ns), fill = shade, horiz = TRUE,
             bty = "n", cex = 0.82)
    })
  }
  
  print(res[, c("mean", "est", "n", "med", "win", "basis_modal", "basis_pct")],
        digits = 3)
  invisible(list(res = res, D = D, ns = ns, M = Mb, freq = freq))
}


# =========================================================================
#  VISUALISATION  --  one realisation of the data-generating mechanism,
#  three panels on a single row (Figure fig_data)
# =========================================================================
make_viz <- function() {
  cat("\n[viz] data-mechanism figure ...\n")
  set.seed(BASESEED + 1L)
  MUv <- MU_FUN(u_viz, t_viz)
  Xi  <- make_Xi()                              # one subject, VIZ_MI visits
  mi  <- VIZ_MI
  Ti  <- sample_T(mi)
  Yi  <- MU_FUN(u_fine, Ti) + Xi(u_fine, Ti) +
    gen_noise(u_fine, mi) * matrix(TAU_FUN(Ti), length(u_fine), mi, byrow = TRUE)
  Zi  <- Xi(u_viz, t_viz)
  wire <- function(Z, zlim, main = "")          # a grey wireframe surface
    persp(u_viz, t_viz, Z, theta = 35, phi = 25, expand = 0.72, zlim = zlim,
          xlab = "u", ylab = "t", zlab = "", ticktype = "detailed",
          cex.axis = 0.55, cex.lab = 0.8, border = NA, col = "grey85",
          shade = 0.55, main = main, cex.main = 0.95)
  emit_fig("fig_data", 10.2, 3.6, function() {
    par(mfrow = c(1, 3))
    par(mar = c(1.0, 1.0, 1.6, 0.6))
    wire(MUv, range(MUv), "$\\mu(u,t)$")
    wire(Zi,  range(Zi),  "$X_i(u,t)$")
    par(mar = c(4.0, 4.2, 1.6, 0.8), mgp = c(2.4, 0.8, 0))
    cols <- grey.colors(mi, start = 0.10, end = 0.60)
    matplot(u_fine, Yi, type = "l", lty = 1, lwd = 0.9, col = cols,
            xlab = "$u$", ylab = "$Y_{ij}(u)$", las = 1, bty = "l",
            main = sprintf("observed curves ($N_u=%d$)", length(u_fine)),
            cex.main = 0.95)
    legend("bottomright", sprintf("$T_{ij}=%.2f$", Ti), col = cols,
           lty = 1, lwd = 1.2, bty = "n", cex = 0.70)
  })
}

# =========================================================================
#  DRIVER  --  run everything in the order of the paper, then print a
#              one-screen summary
# =========================================================================
#  E1 is deterministic (no data at all): it reads the effective regularity of
#  each mean in each basis off the true coefficients, which fixes the rate
#  exponents E6, E7 and E13 are read against. Everything else follows the
#  order of Section~\ref{sec:simulation}: anatomy of the risk (E2-E4), the
#  feasible truncation (E5), the rates (E6-E7), the comparison of the schemes
#  (E8-E11), the augmented basis (E12-E13), the choice of basis (E14) and the
#  fully data-driven estimator against the two benchmarks (E15).
# =========================================================================
t0 <- proc.time()[3]
make_viz()
E1  <- exp_decay()       # effective regularity                (tab_decay)
E2  <- exp_unbias()      # unbiasedness of beta_hat_k          (tab_unbias)
E3  <- exp_riskK()       # risk vs K, kernel plug-in           (tab_riskK, fig_riskK)
E4  <- exp_design()      # design-term decay                   (tab_design, fig_design)
E5  <- exp_cv()          # feasible truncation rule            (tab_trunc)
E6  <- exp_rate()        # rate, Khat*, LL benchmark           (tab_rate, tab_Kstar, fig_rate_total)
E7  <- exp_rate_reg()    # rate under mu_C, LL included        (tab_rate_reg, fig_rate_reg)
E8  <- exp_case12()      # Case 1 vs Case 2                    (tab_case12, fig_case12)
E9  <- exp_boxratio()    # paired log-ratios, B_0 and B_2      (tab_boxratio, fig_boxratio[_aug])
E10 <- exp_plugin()      # true g / plug-in / spacing          (tab_plugin, fig_plugin)
E11 <- exp_robust()      # robustness scenarios                (tab_robust)
E12 <- exp_augbasis()    # augmented basis + sup-norm          (fig_augbasis)
E13 <- exp_augrate()     # rates B_0 / B_2 / LL (mu, mu_C)     (fig_augrate, tab_augrate_mu, tab_augrate_muC)
E14 <- exp_cvbasis()     # basis chosen by CV                  (tab_cvbasis)
E15 <- exp_cvbox()       # data-driven spacing vs LL / spline  (fig_cvbox_mu, _muB, _muC)

cat("\n==================== SUMMARY ====================\n")
cat("[E1] effective regularity s* by (mean, basis):\n")
print(reshape(E1$decay[, c("mean", "r", "s_cap")], idvar = "mean",
              timevar = "r", direction = "wide"), digits = 3)
cat("[E5] the feasible truncation against the empirical oracle:\n")
print(E5[, c("n", "est", "Kstar", "K_med", "infl_med", "R_cv")], digits = 3)
cat("[E6] total-risk exponents (series -> -0.75; LL: leading -0.80):\n")
print(round(E6$exponents, 3))
cat("[E6] Khat* exponents (prediction +0.25) and h_t* exponent (prediction -0.20):\n")
print(round(E6$kexponents, 3)); print(round(E6$hexponent, 3))
cat("[E6] selected h_u* per n (0 = no u-smoothing helps):",
    paste(format(E6$Kres$hu_ll, digits = 3), collapse = "  "), "\n")
cat("[E7] rate exponents under the boundary-smooth mean mu_C (prediction -7/8 = -0.875; LL -0.80):\n")
print(round(E7$exponents["exponent", ], 3))
cat(sprintf("[E8] risk ratio mu^(1)/mu^(2): balanced %.6f (exactly 1), unbalanced %.3f\n",
            E8$res$ratio[1], E8$res$ratio[2]))
cat("[E9] median Delta_e, basis B_0 (positive = spacing wins):\n")
print(reshape(E9$res[E9$res$est != "orac" & E9$res$r == 0L, c("mean", "est", "n", "med")],
              idvar = c("mean", "est"), timevar = "n", direction = "wide"), digits = 3)
cat("[E9] the same in the augmented basis B_2:\n")
print(reshape(E9$res[E9$res$est != "orac" & E9$res$r == 2L, c("mean", "est", "n", "med")],
              idvar = c("mean", "est"), timevar = "n", direction = "wide"), digits = 3)
cat("[E10] plug-in/true-g and spacing/true-g risk ratios:\n")
print(round(E10[, c("n", "r_plug", "r_spac")], 3))
cat("[E11] best scheme per scenario (by MISE at Khat*):\n")
imse_cols <- paste0("MISE_", c("unif", "bal", "MC", "spac"))
print(data.frame(scenario = rownames(E11),
                 best = c("unif", "bal", "MC", "spac")[apply(E11[, imse_cols], 1, which.min)]))
cat("[E4] design-term exponents (prediction -1 / -1 / -2 / -3):\n")
print(round(E4$exponents, 3))
cat(sprintf("[E12] augmented basis: max_k ||phi_k||_inf = %.4f (B_2), %.4f (B_1); 2 sqrt2 = %.4f; argmax at an endpoint for all k: %s\n",
            max(E12$supn), max(E12$supn1), 2 * sqrt(2),
            all(E12$argt < 1e-6 | E12$argt > 1 - 1e-6)))
cat("[E13] rate exponents  spacing B_0 / B_2 / LL  (mu, prediction -0.75 / -0.875 / -0.80):\n")
print(round(E13$mu$sl, 3))
cat("[E13] rate exponents  spacing B_0 / B_2 / LL  (mu_C, prediction -0.875 / -0.875 / -0.80):\n")
print(round(E13$muC$sl, 3))
cat(sprintf("[E13] mu  empirical oracle Khat*: B_0 %d -> %d (grows ~M^{1/4}); B_2 ~%d (frozen ~M^{1/8})\n",
            E13$mu$res$Kstar_sp0[1], tail(E13$mu$res$Kstar_sp0, 1),
            as.integer(round(median(E13$mu$res$Kstar_sp2)))))
cat("[E14] basis selected by CV (rows mu / mu_B / mu_C; * = oracle-best):\n")
print(data.frame(mean = c("mu", "muB", "muC"),
                 r_theory = E14$rpred,
                 r_oracle = E14$bstar - 1L,
                 r_modalCV = apply(E14$freq, 1, which.max) - 1L,
                 pct_modal = round(apply(E14$freq, 1, max), 1),
                 IR_median = round(E14$IR[, 2], 3)))
cat("[E15] fully data-driven spacing (basis + K by CV): median Delta_e\n")
cat("      (positive = the data-driven spacing estimator wins):\n")
print(reshape(E15$res[, c("mean", "est", "n", "med")],
              idvar = c("mean", "est"), timevar = "n", direction = "wide"), digits = 3)
cat("[E15] win rates against the two benchmarks and modal basis chosen:\n")
print(reshape(E15$res[, c("mean", "est", "n", "win", "basis_modal")],
              idvar = c("mean", "est"), timevar = "n", direction = "wide"), digits = 3)
cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
            proc.time()[3] - t0, normalizePath(OUTDIR)))
