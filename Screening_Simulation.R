## =============================================================================
##  Screening_Simulation.R
##
##  Supplementary material to
##  Md. Belaiet Hossain and A Z M Shafiullah "Identifying Active Effects in Unreplicated Two-Level Designs: 
##  A Comparative Simulation Study and a Label-Invariant Qn Scale Estimator", 
##   Submitted to the Journal of Statistical Computation and Simulation (2026).
##    
##
##  Authors       Md. Belaiet Hossain
##                Monetary Policy Department, Bangladesh Bank (the central bank of Bangladesh)
##                Department of Statistics, University of Dhaka, Bangladesh;

##                A Z M Shafiullah
##                Department of Statistics, University of Dhaka, Bangladesh
##                 
##  Correspondence Md. Belaiet Hossain <mdbelaiethossain@gmail.com>
##  Version        1.0, October 2026
##
## -----------------------------------------------------------------------------
##  PURPOSE
##    Reproduces all numerical results of the article: implementation and
##    verification of the screening methods, calibration of every decision rule
##    under the null model, the worked example (Section 5), the Monte Carlo
##    study (Sections 6-7), and all tables and figures.
##
##  REQUIREMENTS
##    R >= 4.0 with the base packages stats, graphics, grDevices and parallel.
##    Optional: BsMD (version 2023.920 used) for the factor-level analysis of
##    Box and Meyer (1993) in scenario S3. Without BsMD that comparator is
##    skipped and all other results are unchanged.
##
##  USAGE
##    1. Set OUT_DIR in Section 0 (default: "Screening_Output" in the
##       working directory).
##    2. Run from a shell:  Rscript Screening_Simulation.R
##       or from R:         source("Screening_Simulation.R")
##    Setting QUICK <- TRUE in Section 0 performs a reduced run for testing the
##    installation; its numbers are not those reported in the article.
##
##  COMPUTING TIME
##    The full run (N_SIM = 10,000; N_NULL = 200,000) took about 75 minutes on
##    two cores (R 4.3.3, Linux): about 45 minutes for calibration and
##    15 minutes for the simulation study.
##
##  REPRODUCIBILITY
##    All random numbers are generated from BASE_SEED. Each simulation
##    condition sets its own seed, so the results do not depend on the number
##    of cores used. Results were reproduced exactly on Linux (R 4.3.3) and on
##    Windows. The script performs 198 automatic checks against published
##    values (Appendix A); the count of passed checks is printed and saved in
##    tables/table_S1_unit_checks.csv.
##
##  OUTPUT (OUT_DIR/tables, OUT_DIR/figures) AND CORRESPONDENCE WITH THE ARTICLE
##    Table 2    table_S4_efficiency.csv, table_S3_breakdown.csv
##    Table 3    table_1_BHH_effects_pvalues.csv, table_2_BHH_decisions.csv
##    Table 4    table_4_null_offtheshelf.csv
##    Table 5    table_S10_k11_delta4.csv
##    Table 6    table_5_share_not_worse_than_best.csv
##    Table 7    table_6_scale_bias_k11.csv
##    Table 8    table_7_interaction.csv
##    Table 9    table_S8_BoxMeyer_sensitivity.csv
##    Figure 1   fig1_BHH_halfnormal_and_BoxMeyer.png
##    Figure 2   fig2_null_error_offtheshelf.png
##    Figure 3   fig3_k11_Power.png
##    Figure 4   fig5_k11_ExactID.png
##    Figure 5   fig6_operating_curves_k11.png
##    Figure 6   fig7_power_by_k_p2.png
##    Figure 7   fig8_interaction_Power_CAL_IER.png
##    Figure 8   fig9_interaction_IER_CAL_IER.png
##    Appendix A table_S1_unit_checks.csv
##    Full results: simulation_results_all.csv, operating_curves_all.csv,
##    paired_vs_best_all.csv. Additional tables (table_S*.csv) and figures
##    (fig4, fig10, figS*.png) report results summarised in the text.
##
##  METHODS (label used in the code -> source)
##    DAN59   Daniel (1959): standardised half-normal plot with guardrails
##    ZAHN75  Zahn (1975a), version S
##    LEN89   Lenth (1989): pseudo standard error
##    SDLEN01 Ye, Hamada and Wu (2001): step-down Lenth method
##    JP92    Juan and Pena (1992): iterated MAD
##    DONG93  Dong (1993): trimmed root mean square, non-iterated form as in
##            Hamada and Balakrishnan (1998)
##    BM86    Box and Meyer (1986): posterior probability (alpha = 0.2, K = 10)
##    MAD0    1.4826 * median|c| (baseline; Lenth's s0 up to a constant)
##    QNS     Proposed: Rousseeuw-Croux Qn on the reflected batch {+c, -c}
##    BM93    Box and Meyer (1993): factor-level analysis via BsMD::BsProb
##            (p = 0.25, gamma = 2.5, up to 3 factors with their two-factor
##            interactions); scenario S3 only
## =============================================================================

## ---------------------------------------------------------------------
## 0. SETTINGS
## ---------------------------------------------------------------------
OUT_DIR   <- "Screening_Output"   # output directory
N_SIM     <- 10000                # replicates per simulation condition
N_NULL    <- 200000               # null replicates used for calibration
N_NULL_BM93 <- 20000              # null replicates to calibrate BM93 (per design)
ALPHA     <- 0.05                 # target error rate for calibrated rules
ALPHA_GRID <- c(0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40)
N_CORES   <- max(1, parallel::detectCores() - 1)   # set to 1 to run serially
BASE_SEED <- 20261011
QUICK     <- FALSE                # TRUE = reduced run for testing only

if (QUICK) { N_SIM <- 400; N_NULL <- 20000; N_NULL_BM93 <- 2000 }
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(OUT_DIR, "figures"); tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, showWarnings = FALSE); dir.create(tab_dir, showWarnings = FALSE)
HAVE_BSMD <- requireNamespace("BsMD", quietly = TRUE)
if (!HAVE_BSMD) message("Package BsMD not installed: the BM93 comparator in scenario S3 is skipped.",
                        " install.packages('BsMD') to include it.")

## tuning constants, as recommended in the original papers
LEN_TRIM <- 2.5; JP_W <- 3.5; JP_AW <- 0.6578; JP_BETA <- 0.05
DONG_TRIM <- 2.5; DONG_SIM <- 0.98
BM_ALPHA <- 0.20; BM_K <- 10; BM_CUT <- 0.5
SEQ_STEP_ALPHA <- 0.05          # per-step level of the published step-down rules
BM93_P <- 0.25; BM93_G <- 2.5   # Box, Hunter & Hunter (2005, Sec. 7.1) values
BM93_MFAC <- 3; BM93_MINT <- 2

METHODS <- c("DAN59","ZAHN75","LEN89","SDLEN01","JP92","DONG93","BM86","MAD0","QNS")
METHOD_LABEL <- c(DAN59 = "Daniel (1959)", ZAHN75 = "Zahn (1975)", LEN89 = "Lenth (1989)",
                  SDLEN01 = "Step-down Lenth (2001)", JP92 = "Juan & Pena (1992)",
                  DONG93 = "Dong (1993)", BM86 = "Box & Meyer (1986)",
                  MAD0 = "MAD0 (baseline)", QNS = "Qn-sym (proposed)",
                  BM93 = "Box & Meyer (1993), factor level")
SEQ_METHODS <- c("DAN59", "ZAHN75", "SDLEN01")

## ---------------------------------------------------------------------
## 1. DESIGNS
## ---------------------------------------------------------------------
make_cyclic_PB <- function(gen){
  k <- length(gen); M <- matrix(NA_real_, k, k); M[, 1] <- gen
  for (j in 2:k) M[, j] <- c(M[k, j - 1], M[1:(k - 1), j - 1])
  rbind(M, rep(-1, k))
}
yates <- function(q){
  n <- 2^q
  sapply(0:(q - 1), function(b) ifelse(((0:(n - 1)) %/% 2^b) %% 2 == 1, 1, -1))
}
yates_cols <- function(F){
  q <- ncol(F); out <- NULL; nm <- NULL
  for (s in 1:(2^q - 1)){
    bits <- which(bitwAnd(s, 2^(0:(q - 1))) > 0)
    out <- cbind(out, apply(F[, bits, drop = FALSE], 1, prod))
    nm <- c(nm, paste(LETTERS[bits], collapse = ""))
  }
  colnames(out) <- nm; out
}
F3 <- yates(3)
X8 <- cbind(A = F3[,1], B = F3[,2], C = F3[,3], D = F3[,1]*F3[,2],
            E = F3[,1]*F3[,3], F = F3[,2]*F3[,3], G = F3[,1]*F3[,2]*F3[,3])
X12 <- make_cyclic_PB(c(1,1,-1,1,1,1,-1,-1,-1,1,-1))         # BHH (2005) Tables 7.1-7.2
colnames(X12) <- c("A","B","C","D","E","F","G","H","J","K","L")
X16 <- yates_cols(yates(4))                                   # Box & Meyer (1986) Table 3
X20 <- make_cyclic_PB(c(1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1))
colnames(X20) <- paste0("X", 1:19)
DESIGNS <- list("n8_2^(7-4)" = X8, "n12_PB" = X12, "n16_2^(15-11)" = X16, "n20_PB" = X20)
DESIGN_K <- sapply(DESIGNS, ncol)
KS <- sort(unique(DESIGN_K))
check_orthogonal <- function(X){
  XtX <- crossprod(X)
  all(XtX[upper.tri(XtX)] == 0) && all(diag(XtX) == nrow(X)) && all(colSums(X) == 0)
}

## ---------------------------------------------------------------------
## 2. EFFECTS
## ---------------------------------------------------------------------
##   contrast sum L_j = x_j'y ; effect estimate c_j = 2 L_j / n ;
##   regression coefficient b_j = L_j / n = c_j / 2 ; Var(c_j) = 4 sigma^2 / n
get_effects <- function(X, y){
  n <- nrow(X); L <- as.vector(crossprod(X, y))
  data.frame(factor = colnames(X), contrast_sum = L, effect = 2 * L / n,
             coefficient = L / n, row.names = NULL)
}
effects_mat <- function(X, Y) t(crossprod(X, Y)) * (2 / nrow(X))

## ---------------------------------------------------------------------
## 3. METHODS (vectorised: rows = data sets, columns = effects)
## ---------------------------------------------------------------------
sort_abs <- function(C){
  C <- if (is.null(dim(C))) matrix(C, nrow = 1) else C
  S <- t(apply(abs(C), 1, sort.int)); if (ncol(C) == 1) S <- t(S); S
}
prefix_median <- function(S, m){
  r <- seq_len(nrow(S)); i1 <- (m + 1) %/% 2; i2 <- m %/% 2 + 1
  (S[cbind(r, i1)] + S[cbind(r, i2)]) / 2
}
row_median <- function(S) prefix_median(S, rep(ncol(S), nrow(S)))
prefix_meansq <- function(S, m){
  cs <- t(apply(S^2, 1, cumsum)); if (ncol(S) == 1) cs <- t(cs)
  cs[cbind(seq_len(nrow(S)), m)] / m
}
## expected half-normal order statistics E|Z|_(i:k) (Zahn 1975a)
half_normal_os_mean <- function(k){
  sapply(1:k, function(i){
    f <- function(x){
      lG <- log(pmax(2 * pnorm(x) - 1, 1e-300))
      l1G <- log(2) + pnorm(x, lower.tail = FALSE, log.p = TRUE)
      t1 <- if (i > 1) (i - 1) * lG else 0
      t2 <- if (k > i) (k - i) * l1G else 0
      x * exp(lchoose(k, i) + log(i) + t1 + t2) * 2 * dnorm(x)
    }
    integrate(f, 0, Inf, rel.tol = 1e-10)$value
  })
}
## Daniel (1959, p. 325): order statistic "numbered nearest (0.683n + 0.5)"
## (11th of 15, 22nd of 31). Note: unrepx (Lenth 2022) uses floor(.683k + .5),
## which gives the 10th of 15; the two rules agree for k = 7, 11, 19.
daniel_index <- function(k) floor(0.683 * k + 1)

tau_lenth <- function(S){
  s0 <- 1.5 * row_median(S); m <- rowSums(S < LEN_TRIM * s0)
  1.5 * prefix_median(S, m)
}
tau_jp <- function(S, w = JP_W, a_w = JP_AW){
  mad <- row_median(S)
  for (it in 1:500){
    m <- rowSums(S <= w * mad); new <- prefix_median(S, m)
    if (isTRUE(all.equal(new, mad, tolerance = 0))) break
    mad <- new
  }
  mad / a_w
}
dong_parts <- function(S){
  s0 <- 1.5 * row_median(S); m <- rowSums(S < DONG_TRIM * s0)
  list(tau = sqrt(prefix_meansq(S, m)), m = m)
}
tau_mad0 <- function(S) row_median(S) / qnorm(0.75)
tau_rms  <- function(S) sqrt(rowMeans(S^2))          # efficient under the null (reference only)
## Qn on the reflected batch: distances ||c_i|-|c_j||, |c_i|+|c_j| (each twice), 2|c_i|
qns_raw <- function(S, chunk = 10000){
  k <- ncol(S); pr <- combn(k, 2); K <- choose(k + 1, 2); out <- numeric(nrow(S))
  for (st in seq(1, nrow(S), by = chunk)){
    idx <- st:min(nrow(S), st + chunk - 1); Si <- S[idx, , drop = FALSE]
    D <- abs(Si[, pr[1, ], drop = FALSE] - Si[, pr[2, ], drop = FALSE])
    P <- Si[, pr[1, ], drop = FALSE] + Si[, pr[2, ], drop = FALSE]
    M <- cbind(D, D, P, P, 2 * Si)
    out[idx] <- apply(M, 1, function(v) sort.int(v, partial = K)[K])
  }
  out
}
qn_signed_raw <- function(c){
  k <- length(c); d <- abs(outer(c, c, "-"))[upper.tri(diag(k))]
  sort(d)[choose(k %/% 2 + 1, 2)]
}
CONST <- new.env()
tau_qns <- function(S) CONST[[paste0("qns_", ncol(S))]] * qns_raw(S)

## Box & Meyer (1986): (i) exact enumeration, (ii) grid integration (eqs 3-5)
bm_post_enum <- function(C, alpha = BM_ALPHA, K = BM_K){
  C <- if (is.null(dim(C))) matrix(C, nrow = 1) else C
  k <- ncol(C); G <- as.matrix(expand.grid(rep(list(0:1), k))); storage.mode(G) <- "double"
  m <- rowSums(G); W <- 1 / (1 + (K^2 - 1) * G); Ssum <- W %*% t(C^2)
  lp <- m * log(alpha / K) + (k - m) * log(1 - alpha) - (k / 2) * log(Ssum)
  lp <- sweep(lp, 2, apply(lp, 2, max)); P <- exp(lp); P <- sweep(P, 2, colSums(P), "/")
  t(crossprod(G, P))
}
bm_post_grid <- function(C, alpha = BM_ALPHA, K = BM_K, G = 500, lo = -9, hi = 2.5, chunk = 20000){
  C <- if (is.null(dim(C))) matrix(C, nrow = 1) else C
  R <- nrow(C); k <- ncol(C); out <- matrix(NA_real_, R, k); ls <- seq(lo, hi, length.out = G)
  for (st in seq(1, R, by = chunk)){
    idx <- st:min(R, st + chunk - 1); Z2 <- C[idx, , drop = FALSE]^2; Z2 <- Z2 / rowMeans(Z2)
    LP <- matrix(NA_real_, length(idx), G)
    for (g in 1:G){
      s2 <- exp(2 * ls[g])
      la <- log(1 - alpha) - Z2 / (2 * s2); lb <- log(alpha / K) - Z2 / (2 * K^2 * s2)
      mx <- pmax(la, lb); LP[, g] <- -k * ls[g] + rowSums(mx + log(exp(la - mx) + exp(lb - mx)))
    }
    W <- exp(LP - apply(LP, 1, max)); W <- W / rowSums(W); acc <- matrix(0, length(idx), k)
    for (g in 1:G){
      s2 <- exp(2 * ls[g])
      la <- log(1 - alpha) - Z2 / (2 * s2); lb <- log(alpha / K) - Z2 / (2 * K^2 * s2)
      acc <- acc + W[, g] / (1 + exp(la - lb))
    }
    out[idx, ] <- acc
  }
  out
}

## ---- sequential (step-down) methods ----
## Each returns an R x q matrix whose column j is the statistic for the j-th
## largest |c|; a method declares the j-th largest active only if steps
## 1..j all exceed their critical values.
seq_q <- function(k, method) if (method == "SDLEN01") floor(k / 2) else k - daniel_index(k)
zahn_den <- function(Sa, mu) as.vector(Sa %*% mu) / sum(mu^2)
seq_stats <- function(S, method){
  k <- ncol(S); q <- seq_q(k, method); a <- daniel_index(k)
  if (method == "DAN59") return(S[, k:(a + 1), drop = FALSE] / S[, a])
  if (method == "ZAHN75") return(S[, k:(a + 1), drop = FALSE] /
                                 zahn_den(S[, 1:a, drop = FALSE], CONST[[paste0("mu_", k)]][1:a]))
  ## SDLEN01, Ye, Hamada & Wu (2001): at step s the PSE is recomputed from the
  ## k - s + 1 contrasts not yet declared active (largest ones removed first)
  matrix(sapply(0:(q - 1), function(s){ j <- k - s; S[, j] / tau_lenth(S[, 1:j, drop = FALSE]) }),
         nrow = nrow(S))
}
## null distribution of each step statistic
seq_null <- function(k, method, N){
  q <- seq_q(k, method); a <- daniel_index(k); out <- matrix(NA_real_, N, q)
  if (method == "DAN59"){                    # all k inert (H&B 1998, Sec. 4.1)
    S0 <- sort_abs(matrix(rnorm(N * k), N, k)); return(S0[, k:(a + 1), drop = FALSE] / S0[, a])
  }
  for (j in 1:q){
    i <- k - j + 1; Si <- sort_abs(matrix(rnorm(N * i), N, i))
    out[, j] <- if (method == "ZAHN75")      # max of i inert (Zahn 1975a)
      Si[, i] / zahn_den(Si[, 1:a, drop = FALSE], CONST[[paste0("mu_", k)]][1:a]) else
      Si[, i] / tau_lenth(Si)                # max of i inert / PSE of those i (YHW 2001)
  }
  out
}
step_down <- function(Tm, crit){
  pass <- Tm[, 1] > crit[1]
  out <- matrix(FALSE, nrow(Tm), ncol(Tm)); out[, 1] <- pass
  if (ncol(Tm) > 1) for (j in 2:ncol(Tm)) { pass <- pass & (Tm[, j] > crit[j]); out[, j] <- pass }
  out
}

## ---- all statistics for a matrix of effects ----
all_stats <- function(C, bm_engine = c("grid", "enum")){
  bm_engine <- match.arg(bm_engine)
  C <- if (is.null(dim(C))) matrix(C, nrow = 1) else C
  S <- sort_abs(C); A <- abs(C); k <- ncol(C)
  ord <- t(apply(A, 1, order, decreasing = TRUE)); if (k == 1) ord <- t(ord)
  tl <- tau_lenth(S); tj <- tau_jp(S); dg <- dong_parts(S); tm <- tau_mad0(S); tq <- tau_qns(S)
  a <- daniel_index(k)
  tz <- zahn_den(S[, 1:a, drop = FALSE], CONST[[paste0("mu_", k)]][1:a])
  list(k = k, ord = ord, dong_m = dg$m,
       tau = cbind(DAN59 = S[, a], ZAHN75 = tz, LEN89 = tl, JP92 = tj, DONG93 = dg$tau,
                   MAD0 = tm, QNS = tq, RMS = tau_rms(S)),
       stat = list(LEN89 = A / tl, JP92 = A / tj, DONG93 = A / dg$tau, MAD0 = A / tm, QNS = A / tq,
                   BM86 = if (bm_engine == "grid") bm_post_grid(C) else bm_post_enum(C)),
       seq = lapply(setNames(SEQ_METHODS, SEQ_METHODS), function(m) seq_stats(S, m)))
}
## decisions (R x k logical, original column order)
decide <- function(st, method, rule){
  k <- st$k; R <- nrow(st$ord); cv <- rule$crit[[method]]
  if (is.null(cv)) return(matrix(NA, R, k))
  if (method %in% SEQ_METHODS){
    q <- ncol(st$seq[[method]]); pass <- step_down(st$seq[[method]], cv)
    D <- matrix(FALSE, R, k)
    D[cbind(rep(seq_len(R), q), as.vector(st$ord[, 1:q, drop = FALSE]))] <- as.vector(pass)
    return(D)
  }
  if (method == "DONG93" && length(cv) > 1) return(st$stat$DONG93 > cv[pmax(st$dong_m, 1)])
  st$stat[[method]] > cv
}

## ---------------------------------------------------------------------
## 3b. Rules and calibration
## ---------------------------------------------------------------------
CAL <- list()
ots_rule <- function(k, null_seq){
  list(name = "OTS", crit = c(list(
    LEN89 = qt(0.975, k / 3),                               # Lenth ME
    JP92 = qnorm((1 + (1 - JP_BETA)^(1 / k)) / 2),          # J&P eq. (8)
    DONG93 = qt((1 + DONG_SIM^(1 / k)) / 2, df = 1:k),      # Dong: t(gamma, m)
    BM86 = BM_CUT),
    lapply(null_seq, function(N0) apply(N0, 2, quantile, 1 - SEQ_STEP_ALPHA, names = FALSE))))
}
ots_extra <- function(k) list(LEN89_SME = qt((1 + 0.95^(1 / k)) / 2, k / 3),
                              DONG93_IND = qt(0.975, df = 1:k))
## rule calibrated so that the null IER equals a  (or null EER, type = "EER")
calib_rule <- function(st0, null_seq, a, type = c("IER", "EER")){
  type <- match.arg(type); crit <- list()
  for (m in setdiff(METHODS, SEQ_METHODS)){
    T0 <- st0$stat[[m]]
    crit[[m]] <- if (type == "IER") quantile(as.vector(T0), 1 - a, names = FALSE) else
                   quantile(apply(T0, 1, max), 1 - a, names = FALSE)
  }
  for (m in SEQ_METHODS){
    crit_at <- function(a1) apply(null_seq[[m]], 2, quantile, 1 - a1, names = FALSE)
    if (type == "EER") { crit[[m]] <- crit_at(a); attr(crit[[m]], "step_alpha") <- a; next }
    ier_at <- function(a1) mean(decide(st0, m, list(crit = setNames(list(crit_at(a1)), m))))
    lo <- 1e-5; hi <- 0.95
    if (ier_at(hi) < a) a1 <- hi else {
      for (it in 1:35){ mid <- (lo + hi) / 2; if (ier_at(mid) > a) hi <- mid else lo <- mid }
      a1 <- (lo + hi) / 2
    }
    crit[[m]] <- crit_at(a1); attr(crit[[m]], "step_alpha") <- a1
  }
  list(name = paste0("CAL_", type), level = a, crit = crit)
}
calibrate_k <- function(k, N = N_NULL, seed = BASE_SEED + k){
  set.seed(seed)
  CONST[[paste0("mu_", k)]] <- half_normal_os_mean(k)
  S0 <- sort_abs(matrix(rnorm(N * k), N, k))
  CONST[[paste0("qns_", k)]] <- 1 / mean(qns_raw(S0))
  C0 <- matrix(rnorm(N * k), N, k); st0 <- all_stats(C0)
  null_seq <- lapply(setNames(SEQ_METHODS, SEQ_METHODS), function(m) seq_null(k, m, N))
  ## Monte Carlo error of the CAL_IER critical values (10 batches)
  b <- rep(1:10, length.out = N)
  mcerr <- sapply(c("LEN89", "JP92", "DONG93", "MAD0", "QNS", "BM86"), function(m)
    sd(tapply(seq_len(N), b, function(ix) quantile(as.vector(st0$stat[[m]][ix, ]), 1 - ALPHA))) / sqrt(10))
  grid <- lapply(ALPHA_GRID, function(a) calib_rule(st0, null_seq, a, "IER"))
  names(grid) <- as.character(ALPHA_GRID)
  CAL[[paste0("k", k)]] <<- list(k = k, st0 = st0, null_seq = null_seq, mcerr = mcerr,
                                 OTS = ots_rule(k, null_seq),
                                 CAL_IER = calib_rule(st0, null_seq, ALPHA, "IER"),
                                 CAL_EER = calib_rule(st0, null_seq, ALPHA, "EER"),
                                 GRID = grid)
  invisible(TRUE)
}
## simulated p-values (pooled null of |c|/tau; BM86: posterior probability)
sim_pvalues <- function(st, k){
  cal <- CAL[[paste0("k", k)]]; res <- list()
  for (m in setdiff(METHODS, SEQ_METHODS)){
    t0 <- sort(as.vector(cal$st0$stat[[m]])); mx <- sort(apply(cal$st0$stat[[m]], 1, max))
    tt <- as.vector(st$stat[[m]])
    res[[m]] <- data.frame(p_ind = 1 - findInterval(tt, t0, left.open = TRUE) / length(t0),
                           p_sim = 1 - findInterval(tt, mx, left.open = TRUE) / length(mx))
  }
  res
}

## ---------------------------------------------------------------------
## 3c. Metrics (with Monte Carlo standard errors)
## ---------------------------------------------------------------------
## truth: TRUE = active main effect; excl: columns excluded from both the
## active and the inert set (factors involved ONLY through an interaction)
metrics <- function(D, truth, excl = rep(FALSE, length(truth))){
  keep <- !excl; D <- D[, keep, drop = FALSE]; truth <- truth[keep]
  k <- length(truth); A <- sum(truth); Z <- k - A; R <- nrow(D)
  TP <- rowSums(D[, truth, drop = FALSE]); FP <- rowSums(D[, !truth, drop = FALSE]); S <- TP + FP
  fdp <- ifelse(S > 0, FP / pmax(S, 1), 0); exact <- (TP == A) & (FP == 0)
  se <- function(x) sd(x) / sqrt(R)
  c(IER = if (Z > 0) mean(FP / Z) else NA, IER_se = if (Z > 0) se(FP / Z) else NA,
    EER = if (Z > 0) mean(FP >= 1) else NA, EER_se = if (Z > 0) se(FP >= 1) else NA,
    FDR = mean(fdp), FDR_se = se(fdp),
    Power = if (A > 0) mean(TP / A) else NA, Power_se = if (A > 0) se(TP / A) else NA,
    ExactID = mean(exact), ExactID_se = se(exact),
    Accuracy = mean((TP + Z - FP) / k), EFP = mean(FP))
}

## ---------------------------------------------------------------------
## 3d. Box & Meyer (1993) factor-level analysis via BsMD (scenario S3)
## ---------------------------------------------------------------------
bm93_post <- function(X, Y){        # Y: n x R ; returns R x k factor posteriors
  t(apply(Y, 2, function(y) BsMD::BsProb(X = X, y = y, mFac = BM93_MFAC, mInt = BM93_MINT,
                                         p = BM93_P, g = BM93_G, nMod = 1)$sprob[-1]))
}
BM93_CUT <- list()

## =====================================================================
## 4. CALIBRATION + CHECKS AGAINST PUBLISHED VALUES
## =====================================================================
cat("\n== Calibrating every method under the null (N_NULL =", N_NULL, ") ==\n")
for (k in KS) { cat("  k =", k, "\n"); calibrate_k(k) }
if (HAVE_BSMD){
  cat("  calibrating BM93 per design\n")
  for (d in names(DESIGNS)){
    set.seed(BASE_SEED + 9300 + DESIGN_K[d]); X <- DESIGNS[[d]]
    P0 <- bm93_post(X, matrix(rnorm(nrow(X) * N_NULL_BM93), nrow(X)))
    BM93_CUT[[d]] <- quantile(as.vector(P0), 1 - ALPHA, names = FALSE)
  }
}

checks <- list()
add_check <- function(what, source, published, ours, tol){
  checks[[length(checks) + 1]] <<- data.frame(check = what, source = source, published = published,
      this_code = round(ours, 4), tolerance = tol, PASS = abs(published - ours) <= tol)
}
for (d in names(DESIGNS)) add_check(paste("orthogonal design", d), "X'X = nI, balanced columns",
                                    1, as.numeric(check_orthogonal(DESIGNS[[d]])), 0)
BHH_y <- c(56, 93, 67, 60, 77, 65, 95, 49, 44, 63, 63, 61)
bhh_pub <- c(5.8, 21.2, -1.5, 7.2, -10.5, -2.2, 2.2, -8.8, 7.2, -0.5, -9.8)
bhh_eff <- get_effects(X12, BHH_y)$effect
for (j in 1:11) add_check(paste0("BHH Table 7.5 effect ", colnames(X12)[j]), "BHH (2005) p.287",
                          bhh_pub[j], bhh_eff[j], 0.05)
for (m in c(7, 15, 31)){
  add_check(paste0("Lenth t(.975, m/3), m=", m), "Lenth (1989) Table 1",
            c(`7` = 3.76, `15` = 2.57, `31` = 2.22)[as.character(m)], qt(.975, m / 3), 0.005)
  add_check(paste0("Lenth t(gamma, m/3), m=", m), "Lenth (1989) Table 1",
            c(`7` = 9.01, `15` = 5.22, `31` = 4.22)[as.character(m)], qt((1 + .95^(1 / m)) / 2, m / 3), 0.005)
}
BM_y <- list(
  I   = c(.23,.30,.52,.54,.70,.76,1.00,.96,.32,.39,.61,.66,.89,.97,1.07,1.21),
  II  = c(43.7,40.2,42.4,44.7,42.4,45.9,42.2,40.6,42.4,45.5,43.6,40.6,44.0,40.2,42.5,46.5),
  III = c(14.0,16.8,15.0,15.4,27.6,24.0,27.4,22.6,22.3,17.1,21.5,17.5,15.9,21.9,16.7,20.3),
  IV  = c(.08,.04,.53,.43,.31,.09,.12,.36,.79,.68,.73,.08,.77,.38,.49,.23))
eII <- get_effects(X16, BM_y$II)$effect; eIV <- get_effects(X16, BM_y$IV)$effect
add_check("Lenth s0, Example II", "Lenth (1989) Sec. 3", 0.45, 1.5 * row_median(sort_abs(eII)), 0.005)
add_check("Lenth PSE, Example II", "Lenth (1989) Sec. 3", 0.225, tau_lenth(sort_abs(eII)), 0.001)
add_check("Lenth PSE, Example IV", "Lenth (1989) Sec. 3", 0.114, tau_lenth(sort_abs(eIV)), 0.001)
aw <- uniroot(function(t) pnorm(t) - 0.5 * pnorm(JP_W * t) - 0.25, c(0.3, qnorm(.75)), tol = 1e-12)$root
add_check("J&P a_w (w = 3.5)", "Juan & Pena (1992) Table I", 0.6578, aw, 0.00005)
for (n in c(7, 15, 31)) add_check(paste0("J&P z_c, n=", n), "Juan & Pena (1992) p.1394",
    c(`7` = 2.68, `15` = 2.93, `31` = 3.15)[as.character(n)], qnorm((1 + .95^(1 / n)) / 2), 0.005)
for (n in c(15, 31, 63)) add_check(paste0("Daniel order statistic, n=", n), "Daniel (1959) p.325",
    c(`15` = 11, `31` = 22, `63` = 44)[as.character(n)], daniel_index(n), 0)
set.seed(BASE_SEED + 31)
Sb <- sort_abs(matrix(rnorm(31 * N_NULL / 4), ncol = 31)); tb <- Sb[, 31] / Sb[, daniel_index(31)]
for (a in c(.05, .10, .20)) add_check(paste0("Birnbaum k(31,", a, ")"), "Birnbaum (1959) Table II (2500 sims)",
    c(`0.05` = 3.36, `0.1` = 3.06, `0.2` = 2.75)[as.character(a)], quantile(tb, 1 - a, names = FALSE), 0.06)
set.seed(BASE_SEED + 93)
for (n in c(7, 15)){
  Sd <- sort_abs(matrix(rnorm(n * N_NULL / 4), ncol = n))
  add_check(paste0("Dong Tab.2 mean s0, n=", n), "Dong (1993) Table 2",
            c(`7` = 1.054, `15` = 1.031)[as.character(n)], mean(1.5 * row_median(Sd)), 0.01)
  add_check(paste0("Dong Tab.2 mean s1, n=", n), "Dong (1993) Table 2",
            c(`7` = 0.890, `15` = 0.915)[as.character(n)], mean(dong_parts(Sd)$tau), 0.01)
  ## tolerance = 2 x Lenth's Monte Carlo SE (SD .43/sqrt(2000), .30/sqrt(1000)) + rounding
  add_check(paste0("Lenth Tab.4 mean PSE (null), m=", n), "Lenth (1989) Table 4",
            c(`7` = 1.05, `15` = 1.01)[as.character(n)], mean(tau_lenth(Sd)), 0.025)
  add_check(paste0("Lenth Tab.4 mean s0 (null), m=", n), "Lenth (1989) Table 4",
            c(`7` = 1.06, `15` = 1.03)[as.character(n)], mean(1.5 * row_median(Sd)), 0.025)
}
add_check("Lenth IER(.05) critical value, k=15", "Ye & Hamada (2000)", 2.156, CAL$k15$CAL_IER$crit$LEN89, 0.02)
add_check("Lenth IER(.05) critical value, k=7", "Ye & Hamada (2000)", 2.30, CAL$k7$CAL_IER$crit$LEN89, 0.02)
## Box & Meyer (1986) Table 4. Examples II-IV from the data; Example I is
## reproduced exactly only from the two-decimal effects printed in Table 4
## (Box & Meyer evidently computed from rounded effects): column 1 = .31 from
## printed effects, .24 from the data (BsMD::BsProb also gives .240).
bm_pub <- list(
  I   = c(.31,1.00,.03,1.00,.02,.03,.02,.99,.05,.03,.03,.09,.03,.03,.03),
  II  = c(.03,.03,.05,.03,.08,.02,.07,.08,.02,.09,.03,.03,.07,1.00,1.00),
  III = c(.05,.03,.05,1.00,.11,.03,.03,.28,.06,.02,.03,1.00,1.00,.02,.05),
  IV  = c(.14,.03,.02,.04,.03,.04,.09,.34,.10,.27,.05,.03,.02,.06,.03))
bm_printed_effects_I <- c(.06,.25,-.01,.50,.00,-.02,.00,.14,.03,-.01,.02,.04,.02,-.01,.02)
for (ex in names(bm_pub)){
  e <- if (ex == "I") bm_printed_effects_I else get_effects(X16, BM_y[[ex]])$effect
  pe <- as.vector(bm_post_enum(e)); pg <- as.vector(bm_post_grid(e))
  tol <- if (ex == "IV") 0.0151 else 0.0051
  for (j in 1:15) add_check(paste0("Box-Meyer Ex.", ex, " column ", j, if (ex == "I") " (printed effects)" else ""),
                            "Box & Meyer (1986) Table 4", bm_pub[[ex]][j], pe[j], tol)
  add_check(paste0("BM grid vs enumeration, Ex.", ex), "internal (max abs diff)", 0, max(abs(pe - pg)), 1e-6)
}
eI <- get_effects(X16, BM_y$I)$effect
add_check("Box-Meyer Ex.I column 1 from the DATA (BsMD gives .240)", "BsMD::BsProb (Meyer's code)",
          0.240, bm_post_enum(eI)[1], 0.0005)
set.seed(BASE_SEED + 7)
for (k in c(7, 11, 15)){
  Cr <- matrix(rnorm(300 * k), 300, k); Cr[, 1:2] <- Cr[, 1:2] + sample(c(0, 4, 12), 600, TRUE)
  add_check(paste0("BM grid vs enumeration, random k=", k), "internal (max abs diff)", 0,
            max(abs(bm_post_grid(Cr) - bm_post_enum(Cr))), 1e-6)
}
for (k in KS) add_check(paste0("sum E|z|_(i), k=", k), "identity", k * sqrt(2 / pi),
                        sum(CONST[[paste0("mu_", k)]]), 1e-6)
## Qn-sym is asymptotically Qn: raw statistic * 2.2191 -> 1 as k grows (Section 3.6, Proposition 2)
set.seed(BASE_SEED + 2219)
qbig <- replicate(if (QUICK) 10 else 40, { z <- abs(rnorm(1001)); y <- c(z, -z)
  d <- abs(outer(y, y, "-"))[upper.tri(diag(2002))]; sort.int(d, partial = choose(1002, 2))[choose(1002, 2)] })
add_check("Qn-sym raw x 2.2191 at k = 1001 (asymptotic consistency)", "Rousseeuw & Croux (1993) constant",
          1, mean(qbig) * 2.2191, 0.01)
set.seed(BASE_SEED + 999)
for (k in KS){
  stv <- all_stats(matrix(rnorm(k * 50000), ncol = k))
  for (m in METHODS){
    Di <- decide(stv, m, CAL[[paste0("k", k)]]$CAL_IER); De <- decide(stv, m, CAL[[paste0("k", k)]]$CAL_EER)
    add_check(paste0("null IER after CAL_IER, ", m, ", k=", k), "target ALPHA (fresh sample)", ALPHA, mean(Di), 0.004)
    add_check(paste0("null EER after CAL_EER, ", m, ", k=", k), "target ALPHA (fresh sample)", ALPHA, mean(apply(De, 1, any)), 0.006)
  }
}
set.seed(BASE_SEED + 1998)
st15 <- all_stats(matrix(rnorm(15 * 50000), ncol = 15)); ots15 <- CAL$k15$OTS
hb <- list(LEN89 = c(.0290, .245), BM86 = c(.0262, .252), JP92 = c(.0294, .201))
for (m in names(hb)){
  D <- decide(st15, m, ots15)
  add_check(paste0("H&B Table 1 null IER, ", m), "Hamada & Balakrishnan (1998) Table 1", hb[[m]][1], mean(D), 0.004)
  add_check(paste0("H&B Table 1 null EER, ", m), "Hamada & Balakrishnan (1998) Table 1", hb[[m]][2], mean(apply(D, 1, any)), 0.015)
}
Dd <- st15$stat$DONG93 > qt(0.975, df = pmax(st15$dong_m, 1))
add_check("H&B Table 1 null IER, DONG93 (t.975,m)", "Hamada & Balakrishnan (1998) Table 1", .0418, mean(Dd), 0.004)
add_check("H&B Table 1 null EER, DONG93 (t.975,m)", "Hamada & Balakrishnan (1998) Table 1", .431, mean(apply(Dd, 1, any)), 0.015)
## step-down Lenth: published rule controls the null EER at .05 (YHW 2001)
add_check("SDLEN01 OTS null EER, k=15", "Ye, Hamada & Wu (2001): EER = .05", 0.05,
          mean(apply(decide(st15, "SDLEN01", ots15), 1, any)), 0.006)
rm(st15)
CHECKS <- do.call(rbind, checks); rownames(CHECKS) <- NULL
write.csv(CHECKS, file.path(tab_dir, "table_S1_unit_checks.csv"), row.names = FALSE)
cat("\nChecks passed:", sum(CHECKS$PASS), "of", nrow(CHECKS), "\n")
if (any(!CHECKS$PASS)) { print(CHECKS[!CHECKS$PASS, ]); warning("Some checks FAILED - see table_S1_unit_checks.csv") }

## critical values, their Monte Carlo error, constants
const_rows <- list()
for (k in KS) for (rule in c("OTS", "CAL_IER", "CAL_EER")) for (m in METHODS){
  cr <- CAL[[paste0("k", k)]][[rule]]$crit[[m]]
  if (is.null(cr)) next
  if (m == "DONG93" && rule == "OTS") cr <- cr[k]
  const_rows[[length(const_rows) + 1]] <- data.frame(k = k, rule = rule, method = m,
    critical_values = paste(round(cr, 4), collapse = "; "),
    step_alpha = if (!is.null(attr(cr, "step_alpha"))) round(attr(cr, "step_alpha"), 4) else NA,
    mc_se = if (rule == "CAL_IER" && m %in% names(CAL[[paste0("k", k)]]$mcerr))
              round(CAL[[paste0("k", k)]]$mcerr[[m]], 4) else NA)
}
write.csv(do.call(rbind, const_rows), file.path(tab_dir, "table_S2_critical_values.csv"), row.names = FALSE)

## ---- 4b. Qn-sym: finite-sample breakdown and efficiency of all scale estimators ----
bd <- do.call(rbind, lapply(KS, function(k){
  K <- choose(k + 1, 2)
  expl <- min(which(sapply(1:k, function(m) (k - m) * (2 * (k - m) - 1) < K)))   # m outliers -> infinity
  impl <- min(which(sapply(1:k, function(m) m * (2 * m - 1) >= K)))              # m effects -> 0
  data.frame(k = k, QNS_explosion = expl, QNS_implosion = impl,
             median_explosion = (k + 1) / 2, median_implosion = (k + 1) / 2,
             QNS_bdp = round(min(expl, impl) / k, 3), median_bdp = round(((k + 1) / 2) / k, 3))
}))
write.csv(bd, file.path(tab_dir, "table_S3_breakdown.csv"), row.names = FALSE)
eff <- do.call(rbind, lapply(KS, function(k){
  T0 <- CAL[[paste0("k", k)]]$st0$tau
  sc <- sweep(T0, 2, colMeans(T0), "/")            # each estimator rescaled to be unbiased
  mse <- colMeans((sc - 1)^2)
  data.frame(k = k, estimator = names(mse), rel_efficiency = round(mse["RMS"] / mse, 3),
             null_mean_unscaled = round(colMeans(T0), 3), null_rmse_unscaled = round(sqrt(colMeans((T0 - 1)^2)), 3))
}))
write.csv(eff, file.path(tab_dir, "table_S4_efficiency.csv"), row.names = FALSE)
write.csv(data.frame(k = KS, daniel_index = sapply(KS, daniel_index),
                     max_declarable_DAN_ZAHN = sapply(KS, function(k) k - daniel_index(k)),
                     max_declarable_SDLEN = floor(KS / 2),
                     qns_constant = sapply(KS, function(k) CONST[[paste0("qns_", k)]])),
          file.path(tab_dir, "table_S5_constants.csv"), row.names = FALSE)

## =====================================================================
## 5. WORKED EXAMPLE: BHH (2005) 12-run PB reactor data (Table 7.5)
## =====================================================================
cat("\n== Worked example: BHH 12-run PB reactor data ==\n")
ex_tab <- get_effects(X12, BHH_y)
stopifnot(all(abs(ex_tab$coefficient - ex_tab$effect / 2) < 1e-12))
st_ex <- all_stats(matrix(ex_tab$effect, nrow = 1), bm_engine = "enum")
pv <- sim_pvalues(st_ex, 11); ex_out <- ex_tab
for (m in c("LEN89", "JP92", "DONG93", "MAD0", "QNS")){
  ex_out[[paste0(m, "_stat")]] <- round(as.vector(st_ex$stat[[m]]), 3)
  ex_out[[paste0(m, "_p")]] <- round(pv[[m]]$p_ind, 4)
}
ex_out$BM86_post <- round(as.vector(st_ex$stat$BM86), 3)
if (HAVE_BSMD) ex_out$BM93_factor_post <- round(as.vector(bm93_post(X12, matrix(BHH_y))), 3)
write.csv(ex_out, file.path(tab_dir, "table_1_BHH_effects_pvalues.csv"), row.names = FALSE)
dec_rows <- list()
for (rule in c("OTS", "CAL_IER", "CAL_EER")) for (m in METHODS){
  D <- decide(st_ex, m, CAL$k11[[rule]])
  dec_rows[[length(dec_rows) + 1]] <- data.frame(rule = rule, method = m,
     tau_hat = if (m %in% colnames(st_ex$tau)) round(st_ex$tau[1, m], 3) else NA,
     declared = if (any(is.na(D))) "no rule" else paste(colnames(X12)[which(D[1, ])], collapse = ","))
}
if (HAVE_BSMD) dec_rows[[length(dec_rows) + 1]] <- data.frame(rule = c("OTS", "CAL_IER"), method = "BM93", tau_hat = NA,
  declared = c(paste(colnames(X12)[ex_out$BM93_factor_post > 0.5], collapse = ","),
               paste(colnames(X12)[ex_out$BM93_factor_post > BM93_CUT[["n12_PB"]]], collapse = ",")))
EX_DEC <- do.call(rbind, dec_rows); EX_DEC$declared[EX_DEC$declared == ""] <- "none"
write.csv(EX_DEC, file.path(tab_dir, "table_2_BHH_decisions.csv"), row.names = FALSE)
print(EX_DEC[EX_DEC$rule == "CAL_IER", ])
seq_tab <- do.call(rbind, lapply(SEQ_METHODS, function(m){
  q <- ncol(st_ex$seq[[m]]); ordj <- st_ex$ord[1, 1:q]
  data.frame(method = m, step = 1:q, factor = colnames(X12)[ordj], ratio = round(as.vector(st_ex$seq[[m]]), 3),
             crit_OTS = round(CAL$k11$OTS$crit[[m]], 3), crit_CAL_IER = round(CAL$k11$CAL_IER$crit[[m]], 3),
             p_step = sapply(1:q, function(j) mean(CAL$k11$null_seq[[m]][, j] >= st_ex$seq[[m]][1, j])))
}))
seq_tab$p_stepdown <- ave(seq_tab$p_step, seq_tab$method, FUN = cummax)
write.csv(seq_tab, file.path(tab_dir, "table_3_BHH_stepdown.csv"), row.names = FALSE)
set.seed(BASE_SEED + 5)
flip_res <- t(replicate(200, { s <- sample(c(-1, 1), 11, TRUE); cc <- ex_tab$effect * s
  c(qn_signed = qn_signed_raw(cc), qn_reflected = qns_raw(matrix(abs(cc), 1))) }))
write.csv(data.frame(statistic = c("Qn on signed contrasts (raw)", "Qn on reflected batch (raw)"),
                     min_over_relabellings = apply(flip_res, 2, min), max_over_relabellings = apply(flip_res, 2, max)),
          file.path(tab_dir, "table_S6_label_invariance.csv"), row.names = FALSE)

## =====================================================================
## 6. SIMULATION STUDY
## =====================================================================
## y = X(Delta/2) + interaction terms + e, e ~ N(0, n/4): every effect
## estimate has SE tau = 1 and an active main effect has E(c_j) = Delta.
DELTAS <- c(1, 2, 3, 4, 5, 6, 8); GAMMAS <- c(0, 2, 4, 6, 8)
## columns used in the interaction scenarios: A = col 1, B = col 2; the
## other active factors are the first columns NOT fully aliased with AB
other_cols <- function(X){ ab <- X[, 1] * X[, 2]
  setdiff(which(abs(as.vector(crossprod(X, ab))) < nrow(X)), 1:2) }
scen <- list()
add_scen <- function(...) scen[[length(scen) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
for (d in names(DESIGNS)){
  k <- unname(DESIGN_K[d]); pmax_d <- if (k <= 7) 3 else 4
  add_scen(scenario = "S1_main_equal", design = d, p = 0, delta = 0, sizes = "equal", errors = "normal", inter = "none", gamma = 0)
  for (p in 1:pmax_d) for (dl in DELTAS)
    add_scen(scenario = "S1_main_equal", design = d, p = p, delta = dl, sizes = "equal", errors = "normal", inter = "none", gamma = 0)
  for (it in c("a_strong", "b_weak", "c_none", "d_two")) for (g in GAMMAS)
    add_scen(scenario = "S3_interaction", design = d, p = NA, delta = 4, sizes = "equal", errors = "normal", inter = it, gamma = g)
}
for (p in 2:4) for (dl in DELTAS)
  add_scen(scenario = "S2_main_unequal", design = "n12_PB", p = p, delta = dl, sizes = "unequal", errors = "normal", inter = "none", gamma = 0)
add_scen(scenario = "S4_t5_errors", design = "n12_PB", p = 0, delta = 0, sizes = "equal", errors = "t5", inter = "none", gamma = 0)
for (p in 1:3) for (dl in DELTAS)
  add_scen(scenario = "S4_t5_errors", design = "n12_PB", p = p, delta = dl, sizes = "equal", errors = "t5", inter = "none", gamma = 0)
SCEN <- do.call(rbind, scen); SCEN$id <- seq_len(nrow(SCEN))
write.csv(SCEN, file.path(tab_dir, "table_S7_scenarios.csv"), row.names = FALSE)
cat("\n== Simulation:", nrow(SCEN), "conditions x", N_SIM, "replicates ==\n")

## true model of a condition: main-effect means, interaction columns, truth, excl
true_model <- function(sc, X){
  k <- ncol(X); dv <- rep(0, k); excl <- rep(FALSE, k); icol <- list(); g <- sc$gamma
  if (sc$scenario != "S3_interaction"){
    if (sc$p > 0) dv[1:sc$p] <- if (sc$sizes == "equal") rep(sc$delta, sc$p) else sc$delta * seq(1, 0.5, length.out = sc$p)
  } else {
    oc <- other_cols(X)
    if (sc$inter == "a_strong") { dv[1:2] <- sc$delta; icol <- list(c(1, 2)) }               # A, B, AB
    if (sc$inter == "b_weak")   { dv[1] <- sc$delta; icol <- list(c(1, 2)); excl[2] <- g > 0 } # A, AB
    if (sc$inter == "c_none")   { dv[oc[1:2]] <- sc$delta; icol <- list(c(1, 2)); excl[1:2] <- g > 0 } # C, D, AB
    if (sc$inter == "d_two")    { dv[c(1, 2, oc[1])] <- sc$delta; icol <- list(c(1, 2), c(1, oc[1])) } # A,B,C, AB, AC
  }
  mu <- as.vector(X %*% (dv / 2))
  for (ij in icol) mu <- mu + (g / 2) * X[, ij[1]] * X[, ij[2]]
  list(mu = mu, truth = dv != 0, excl = excl)
}

run_condition <- function(sc){
  set.seed(BASE_SEED + 1000 * sc$id)
  X <- DESIGNS[[sc$design]]; n <- nrow(X); k <- ncol(X); sigma <- sqrt(n) / 2
  tm <- true_model(sc, X)
  E <- if (sc$errors == "normal") matrix(rnorm(n * N_SIM, 0, sigma), n) else
         matrix(rt(n * N_SIM, df = 5) * sqrt(3 / 5) * sigma, n)
  Y <- tm$mu + E; C <- effects_mat(X, Y); st <- all_stats(C)
  cal <- CAL_LIGHT[[paste0("k", k)]]
  base <- data.frame(id = sc$id, scenario = sc$scenario, design = sc$design, k = k, p = sc$p, delta = sc$delta,
                     sizes = sc$sizes, errors = sc$errors, inter = sc$inter, gamma = sc$gamma)
  main <- list(); oc <- list(); pw <- list()
  for (rule in c("OTS", "CAL_IER", "CAL_EER")) for (m in METHODS){
    D <- decide(st, m, cal[[rule]]); if (any(is.na(D))) next
    main[[length(main) + 1]] <- cbind(base, rule = rule, method = m, t(metrics(D, tm$truth, tm$excl)),
      tau_mean = if (m %in% colnames(st$tau)) mean(st$tau[, m]) else NA,
      tau_rmse = if (m %in% colnames(st$tau)) sqrt(mean((st$tau[, m] - 1)^2)) else NA)
    if (rule == "CAL_IER" && any(tm$truth)) pw[[m]] <- rowMeans(D[, tm$truth & !tm$excl, drop = FALSE])
  }
  if (sc$scenario == "S3_interaction" && HAVE_BSMD){
    P93 <- bm93_post(X, Y)
    for (rule in c("OTS", "CAL_IER")){
      D <- P93 > if (rule == "OTS") 0.5 else BM93_CUT[[sc$design]]
      main[[length(main) + 1]] <- cbind(base, rule = rule, method = "BM93", t(metrics(D, tm$truth, tm$excl)),
                                        tau_mean = NA, tau_rmse = NA)
    }
  }
  ## operating curves: every method at every null IER level
  if (sc$scenario %in% c("S1_main_equal", "S3_interaction"))
    for (a in names(cal$GRID)) for (m in METHODS){
      mt <- metrics(decide(st, m, cal$GRID[[a]]), tm$truth, tm$excl)
      oc[[length(oc) + 1]] <- cbind(base, null_IER = as.numeric(a), method = m,
                                    t(mt[c("Power", "IER", "EER", "FDR", "EFP", "ExactID")]))
    }
  ## paired comparison with the most powerful method (same data sets)
  pair <- NULL
  if (length(pw)){
    pm <- sapply(pw, mean); best <- names(which.max(pm))
    pair <- do.call(rbind, lapply(names(pw), function(m){ d <- pw[[m]] - pw[[best]]
      data.frame(base, method = m, best = best, diff = mean(d), diff_se = sd(d) / sqrt(length(d)),
                 not_worse = mean(d) >= -2 * sd(d) / sqrt(length(d))) }))
  }
  list(main = do.call(rbind, main), oc = if (length(oc)) do.call(rbind, oc) else NULL, pair = pair)
}
CAL_LIGHT <- lapply(CAL, function(z) z[c("k", "OTS", "CAL_IER", "CAL_EER", "GRID")])
CONST_LIST <- as.list(CONST)
t_start <- Sys.time()
if (N_CORES > 1){
  cl <- parallel::makeCluster(N_CORES)
  parallel::clusterExport(cl, setdiff(ls(globalenv()), c("CAL", "cl")), envir = globalenv())
  parallel::clusterEvalQ(cl, { CONST <- list2env(CONST_LIST); NULL })
  RES <- parallel::parLapplyLB(cl, split(SCEN, SCEN$id), run_condition)
  parallel::stopCluster(cl)
} else RES <- lapply(split(SCEN, SCEN$id), run_condition)
SIM  <- do.call(rbind, lapply(RES, `[[`, "main")); rownames(SIM) <- NULL
OC   <- do.call(rbind, lapply(RES, `[[`, "oc"));   rownames(OC) <- NULL
PAIR <- do.call(rbind, lapply(RES, `[[`, "pair")); rownames(PAIR) <- NULL
cat("Simulation time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
write.csv(SIM, file.path(tab_dir, "simulation_results_all.csv"), row.names = FALSE)
write.csv(OC, file.path(tab_dir, "operating_curves_all.csv"), row.names = FALSE)
write.csv(PAIR, file.path(tab_dir, "paired_vs_best_all.csv"), row.names = FALSE)

## ---- 6b. Box & Meyer (1986): sensitivity to (alpha, K), k = 11 ----
cat("\n== Box-Meyer sensitivity ==\n")
BM_GRID <- list(c(0.1, 10), c(0.2, 5), c(0.2, 10), c(0.2, 15), c(0.3, 10))
set.seed(BASE_SEED + 8686)
C0 <- matrix(rnorm(11 * (N_NULL / 4)), ncol = 11)
bm_sens <- do.call(rbind, lapply(BM_GRID, function(ak){
  P0 <- bm_post_grid(C0, ak[1], ak[2]); cut <- quantile(as.vector(P0), 1 - ALPHA, names = FALSE)
  do.call(rbind, lapply(1:4, function(p) do.call(rbind, lapply(c(2, 4, 6), function(dl){
    C <- matrix(rnorm(11 * N_SIM), ncol = 11); C[, 1:p] <- C[, 1:p] + dl
    P <- bm_post_grid(C, ak[1], ak[2]); tr <- c(rep(TRUE, p), rep(FALSE, 11 - p))
    rbind(data.frame(alpha = ak[1], K = ak[2], p = p, delta = dl, rule = "OTS (post > .5)", t(metrics(P > 0.5, tr))),
          data.frame(alpha = ak[1], K = ak[2], p = p, delta = dl, rule = "CAL_IER", t(metrics(P > cut, tr))))
  }))))
}))
write.csv(bm_sens, file.path(tab_dir, "table_S8_BoxMeyer_sensitivity.csv"), row.names = FALSE)

## =====================================================================
## 7. TABLES AND FIGURES
## =====================================================================
ALLM <- c(METHODS, "BM93")
COL <- c(DAN59 = "#000000", ZAHN75 = "#8C8C8C", LEN89 = "#E69F00", SDLEN01 = "#A6761D", JP92 = "#56B4E9",
         DONG93 = "#009E73", BM86 = "#0072B2", MAD0 = "#CC79A7", QNS = "#D55E00", BM93 = "#6A3D9A")
LTY <- c(DAN59 = 1, ZAHN75 = 2, LEN89 = 1, SDLEN01 = 3, JP92 = 4, DONG93 = 2, BM86 = 1, MAD0 = 3, QNS = 1, BM93 = 5)
PCH <- c(DAN59 = 15, ZAHN75 = 0, LEN89 = 17, SDLEN01 = 25, JP92 = 2, DONG93 = 16, BM86 = 18, MAD0 = 4, QNS = 8, BM93 = 11)
LWD <- c(DAN59 = 1.5, ZAHN75 = 1.5, LEN89 = 1.5, SDLEN01 = 1.5, JP92 = 1.5, DONG93 = 1.5, BM86 = 1.5, MAD0 = 1.5, QNS = 2.6, BM93 = 2)
DLAB <- c("n8_2^(7-4)" = "8-run regular 2^(7-4)", "n12_PB" = "12-run Plackett-Burman",
          "n16_2^(15-11)" = "16-run regular 2^(15-11)", "n20_PB" = "20-run Plackett-Burman")
ILAB <- c(a_strong = "A, B active + AB", b_weak = "A active + AB", c_none = "C, D active + AB (A, B inert)",
          d_two = "A, B, C active + AB, AC")
mlab <- function(m) METHOD_LABEL[m]
panel_plot <- function(dat, metric, ylab, file, panels, panel_var, xvar = "delta",
                       xlab = expression("Active effect size " * Delta * " (SE units)"),
                       ylim = NULL, titles = NULL, ms = METHODS, hline = NULL, ncol_max = 4){
  ms <- ms[ms %in% unique(dat$method)]
  nr <- ceiling(length(panels) / ncol_max); nc <- min(ncol_max, length(panels))
  png(file, width = 420 * nc + 90, height = 440 * nr + 170, res = 140)
  par(mfrow = c(nr, nc), mar = c(4.2, 4.2, 2.4, 0.6), oma = c(4.6, 0, 0, 0))
  for (i in seq_along(panels)){
    sub <- dat[dat[[panel_var]] == panels[i], ]
    yl <- if (is.null(ylim)) range(c(0, sub[[metric]]), na.rm = TRUE) else ylim
    plot(NULL, xlim = range(sub[[xvar]]), ylim = yl, xlab = xlab, ylab = ylab,
         main = if (is.null(titles)) panels[i] else titles[i], cex.main = 0.9, cex.lab = 0.85)
    if (!is.null(hline)) abline(h = hline, col = "grey70", lty = 3)
    for (m in ms){ d <- sub[sub$method == m, ]; d <- d[order(d[[xvar]]), ]
      lines(d[[xvar]], d[[metric]], col = COL[m], lty = LTY[m], lwd = LWD[m])
      points(d[[xvar]], d[[metric]], col = COL[m], pch = PCH[m], cex = 0.65) }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = mlab(ms), col = COL[ms], lty = LTY[ms], pch = PCH[ms], lwd = LWD[ms],
         ncol = min(5, length(ms)), cex = 0.75, bty = "n", xpd = TRUE)
  invisible(dev.off())
}
ptitle <- function(ps, k = 11) paste0("k = ", k, ", p = ", ps, " active")

## ---- Fig 1: worked example ----
png(file.path(fig_dir, "fig1_BHH_halfnormal_and_BoxMeyer.png"), width = 1650, height = 640, res = 140)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.6, 1))
ae <- abs(ex_tab$effect); o <- order(ae); mu11 <- CONST$mu_11
plot(mu11, ae[o], pch = 19, xlab = "Expected half-normal order statistic", ylab = "|effect estimate|",
     main = "(a) Half-normal plot, BHH 12-run PB data", cex.main = 0.95, xlim = c(0, 2.3), ylim = c(0, 23))
text(mu11, ae[o], ex_tab$factor[o], pos = 2, cex = 0.8)
abline(0, st_ex$tau[1, "ZAHN75"], lty = 2); abline(0, st_ex$tau[1, "LEN89"], lty = 3, col = COL["LEN89"])
legend("topleft", c("Zahn slope", "Lenth PSE slope"), lty = c(2, 3), col = c("black", COL["LEN89"]), bty = "n", cex = 0.8)
text(2.3, 1.5, "Full 2^5 analysis (BHH p.287):\nB, D, E, BD, DE active", adj = c(1, 0), cex = 0.75)
bp <- rbind(as.vector(st_ex$stat$BM86), if (HAVE_BSMD) ex_out$BM93_factor_post else NULL)
barplot(bp, beside = TRUE, names.arg = ex_tab$factor, ylim = c(0, 1), col = c("grey60", "#6A3D9A")[1:nrow(bp)],
        ylab = "Posterior probability", main = "(b) Box & Meyer posteriors", cex.main = 0.95)
abline(h = 0.5, lty = 2)
legend("topright", c("1986, effect level (alpha=.2, K=10)", "1993, factor level (p=.25, g=2.5)")[1:nrow(bp)],
       fill = c("grey60", "#6A3D9A")[1:nrow(bp)], bty = "n", cex = 0.75)
invisible(dev.off())

## ---- Table 4 / Fig 2: null model, published rules ----
null_ots <- SIM[SIM$scenario == "S1_main_equal" & SIM$p == 0 & SIM$rule == "OTS",
                c("k", "method", "IER", "IER_se", "EER", "EER_se")]
set.seed(BASE_SEED + 4242)
extra <- do.call(rbind, lapply(KS, function(k){
  stn <- all_stats(matrix(rnorm(N_SIM * k), ncol = k)); ex <- ots_extra(k)
  D1 <- stn$stat$LEN89 > ex$LEN89_SME; D2 <- stn$stat$DONG93 > ex$DONG93_IND[pmax(stn$dong_m, 1)]
  data.frame(k = k, method = c("LEN89 (SME rule)", "DONG93 (individual t(.975, m))"),
             IER = c(mean(D1), mean(D2)), IER_se = NA,
             EER = c(mean(apply(D1, 1, any)), mean(apply(D2, 1, any))), EER_se = NA)
}))
NULL_OTS <- rbind(null_ots, extra); NULL_OTS <- NULL_OTS[order(NULL_OTS$k, NULL_OTS$method), ]
write.csv(NULL_OTS, file.path(tab_dir, "table_4_null_offtheshelf.csv"), row.names = FALSE)
ms_ots <- c("DAN59", "ZAHN75", "LEN89", "SDLEN01", "JP92", "DONG93", "BM86")
png(file.path(fig_dir, "fig2_null_error_offtheshelf.png"), width = 1650, height = 600, res = 140)
par(mfrow = c(1, 2), mar = c(7.8, 4.2, 2.4, 0.6))
for (met in c("IER", "EER")){
  M <- sapply(KS, function(k) sapply(ms_ots, function(m) null_ots[null_ots$k == k & null_ots$method == m, met]))
  barplot(t(M), beside = TRUE, names.arg = mlab(ms_ots), las = 2, cex.names = 0.72,
          col = c("grey85", "grey60", "grey40", "grey15"), ylab = met, main = paste("Null model, published rules:", met))
  abline(h = ALPHA, lty = 2, col = "red")
  if (met == "IER") legend("topleft", legend = paste("k =", KS), fill = c("grey85", "grey60", "grey40", "grey15"), bty = "n", cex = 0.75)
}
invisible(dev.off())

## ---- Figs 3-5: k = 11 main effects, CAL_IER (depends on k only) ----
s1 <- SIM[SIM$scenario == "S1_main_equal" & SIM$p > 0, ]
k11 <- s1[s1$k == 11 & s1$rule == "CAL_IER", ]
panel_plot(k11, "Power", "Power", file.path(fig_dir, "fig3_k11_Power.png"), 1:4, "p", titles = ptitle(1:4), ylim = c(0, 1))
panel_plot(k11, "IER", "IER (non-null)", file.path(fig_dir, "fig4_k11_IER.png"), 1:4, "p", titles = ptitle(1:4), hline = ALPHA)
panel_plot(k11, "ExactID", "P(exactly the true set)", file.path(fig_dir, "fig5_k11_ExactID.png"), 1:4, "p", titles = ptitle(1:4), ylim = c(0, 1))
panel_plot(k11, "FDR", "FDR", file.path(fig_dir, "figS1_k11_FDR.png"), 1:4, "p", titles = ptitle(1:4))
panel_plot(k11, "EER", "EER", file.path(fig_dir, "figS2_k11_EER.png"), 1:4, "p", titles = ptitle(1:4))
panel_plot(s1[s1$k == 11 & s1$rule == "CAL_EER", ], "Power", "Power (CAL_EER)", file.path(fig_dir, "figS3_k11_Power_calEER.png"),
           1:4, "p", titles = ptitle(1:4), ylim = c(0, 1))
panel_plot(s1[s1$k == 11 & s1$rule == "OTS", ], "Power", "Power (published rules)", file.path(fig_dir, "figS4_k11_Power_OTS.png"),
           1:4, "p", titles = ptitle(1:4), ylim = c(0, 1), ms = ms_ots)

## ---- Fig 6: operating curves, k = 11, Delta = 4 ----
oc11 <- OC[OC$scenario == "S1_main_equal" & OC$k == 11 & OC$delta == 4, ]
png(file.path(fig_dir, "fig6_operating_curves_k11.png"), width = 1770, height = 640, res = 140)
par(mfrow = c(1, 4), mar = c(4.2, 4.2, 2.4, 0.6), oma = c(4.6, 0, 0, 0))
for (p in 1:4){
  d <- oc11[oc11$p == p, ]
  plot(NULL, xlim = c(0, max(d$EFP)), ylim = c(0, 1), xlab = "Mean number of false positives",
       ylab = "Power", main = ptitle(p), cex.main = 0.9)
  for (m in METHODS){ dm <- d[d$method == m, ]; dm <- dm[order(dm$null_IER), ]
    lines(dm$EFP, dm$Power, col = COL[m], lty = LTY[m], lwd = LWD[m])
    i5 <- which(abs(dm$null_IER - 0.05) < 1e-9); points(dm$EFP[i5], dm$Power[i5], col = COL[m], pch = PCH[m]) }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = mlab(METHODS), col = COL[METHODS], lty = LTY[METHODS], pch = PCH[METHODS],
       lwd = LWD[METHODS], ncol = 5, cex = 0.75, bty = "n", xpd = TRUE)
invisible(dev.off())

## ---- Fig 7: power by k ----
panel_plot(s1[s1$p == 2 & s1$rule == "CAL_IER", ], "Power", "Power (p = 2)", file.path(fig_dir, "fig7_power_by_k_p2.png"),
           KS, "k", titles = paste("k =", KS), ylim = c(0, 1))
panel_plot(s1[s1$p == 3 & s1$rule == "CAL_IER", ], "Power", "Power (p = 3)", file.path(fig_dir, "figS5_power_by_k_p3.png"),
           KS, "k", titles = paste("k =", KS), ylim = c(0, 1))

## ---- Figs 8-9: interaction scenarios ----
s3 <- SIM[SIM$scenario == "S3_interaction", ]; s3$panel <- paste(s3$inter, s3$design)
for (rule in c("CAL_IER", "OTS")){
  d3 <- s3[s3$rule == rule, ]
  pan <- as.vector(t(outer(names(ILAB), names(DESIGNS), paste)))
  tt  <- as.vector(t(outer(ILAB, DLAB[names(DESIGNS)], function(a, b) paste0(b, ":\n", a))))
  ms3 <- if (rule == "OTS") c(ms_ots, "BM93") else ALLM
  tag <- if (rule == "CAL_IER") c("fig8", "fig9") else c("figS6", "figS7")
  panel_plot(d3, "Power", "Power for active main effects", file.path(fig_dir, paste0(tag[1], "_interaction_Power_", rule, ".png")),
             pan, "panel", xvar = "gamma", xlab = expression("Interaction size " * gamma), titles = tt, ms = ms3, ylim = c(0, 1))
  panel_plot(d3, "IER", "IER (factors with no effect)", file.path(fig_dir, paste0(tag[2], "_interaction_IER_", rule, ".png")),
             pan, "panel", xvar = "gamma", xlab = expression("Interaction size " * gamma), titles = tt, ms = ms3, hline = ALPHA)
}
write.csv(s3[, c("design", "inter", "gamma", "rule", "method", "Power", "Power_se", "IER", "IER_se", "EER", "FDR", "ExactID")],
          file.path(tab_dir, "table_7_interaction.csv"), row.names = FALSE)

## ---- unequal effects, t5 errors ----
panel_plot(SIM[SIM$scenario == "S2_main_unequal" & SIM$rule == "CAL_IER", ], "Power", "Power (unequal sizes)",
           file.path(fig_dir, "figS8_k11_unequal_Power.png"), 2:4, "p", titles = ptitle(2:4), ylim = c(0, 1),
           xlab = expression("Largest active effect " * Delta))
s4 <- SIM[SIM$scenario == "S4_t5_errors" & SIM$p > 0 & SIM$rule == "CAL_IER", ]
panel_plot(s4, "Power", "Power (t5 errors)", file.path(fig_dir, "figS9_k11_t5_Power.png"), 1:3, "p", titles = ptitle(1:3), ylim = c(0, 1))
write.csv(SIM[SIM$scenario == "S4_t5_errors" & SIM$p == 0, c("rule", "method", "IER", "IER_se", "EER", "EER_se")],
          file.path(tab_dir, "table_S9_t5_null.csv"), row.names = FALSE)

## ---- scale-estimator bias (Fig 10) ----
sc_m <- c("DAN59", "ZAHN75", "LEN89", "JP92", "DONG93", "MAD0", "QNS")
bias <- SIM[SIM$scenario == "S1_main_equal" & SIM$k == 11 & SIM$rule == "CAL_IER" & SIM$delta %in% c(0, 6) &
            SIM$method %in% sc_m, c("p", "delta", "method", "tau_mean", "tau_rmse")]
write.csv(bias, file.path(tab_dir, "table_6_scale_bias_k11.csv"), row.names = FALSE)
png(file.path(fig_dir, "fig10_scale_bias_k11.png"), width = 1100, height = 760, res = 140)
par(mar = c(4.2, 4.2, 2.4, 0.6))
plot(NULL, xlim = c(0, 4), ylim = range(bias$tau_mean), xlab = "Number of active effects p (Delta = 6)",
     ylab = expression("Mean of " * hat(tau) * " (true " * tau * " = 1)"), main = "Bias of the scale estimators, k = 11")
abline(h = 1, col = "grey60", lty = 3)
for (m in sc_m){ d <- bias[bias$method == m, ]; d <- d[order(d$p), ]
  lines(d$p, d$tau_mean, col = COL[m], lty = LTY[m], lwd = LWD[m]); points(d$p, d$tau_mean, col = COL[m], pch = PCH[m]) }
legend("topleft", legend = mlab(sc_m), col = COL[sc_m], lty = LTY[sc_m], pch = PCH[sc_m], lwd = LWD[sc_m], cex = 0.75, bty = "n")
invisible(dev.off())

## ---- Table 5: paired comparison with the best method ----
pr <- PAIR[PAIR$scenario == "S1_main_equal", ]
pr$pgrp <- ifelse(pr$p <= 2, "p = 1-2", ifelse(pr$p == 3, "p = 3", "p = 4"))
T5 <- aggregate(not_worse ~ k + pgrp + method, data = pr, FUN = mean)
write.csv(T5[order(T5$k, T5$pgrp, -T5$not_worse), ], file.path(tab_dir, "table_5_share_not_worse_than_best.csv"), row.names = FALSE)
## compact grid for k = 11 at Delta = 4
g4 <- k11[k11$delta == 4, c("p", "method", "Power", "Power_se", "EER", "ExactID", "ExactID_se")]
write.csv(g4[order(g4$p, -g4$Power), ], file.path(tab_dir, "table_S10_k11_delta4.csv"), row.names = FALSE)

cat("\nAll done. Tables in", tab_dir, "and figures in", fig_dir, "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
