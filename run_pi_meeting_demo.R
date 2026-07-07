# ============================================================
#  PI MEETING DEMO: Olivia's numerical lambda + Kato validation
#  Real resolution (n=300, matching the thesis pipeline).
#  Total runtime: ~15 minutes on a normal RStudio Cloud instance.
#  Just click "Source" (or Ctrl+Shift+Enter / Cmd+Shift+Enter) and wait.
#  Progress messages will print in the Console as it goes.
# ============================================================

# ---- 0. Settings ----
# n_mesh = 300 matches the real thesis/manuscript resolution.
# n_years/burn_yrs are SHORTER than the production 70yr/40yr-burn setup,
# purely to keep this runnable in ~15 min for tomorrow's meeting.
# (We confirmed separately that runtime here is essentially independent
#  of n_mesh -- it's dominated by the 365-day-per-year loop -- so using
#  the real n=300 costs nothing extra. Only n_years drives runtime.)
n_mesh   <- 300
n_years  <- 4
burn_yrs <- 1

# ---- 1. Load parameters and build the size grid ----
source("Code/Salmonparams.R")
min.size <- log(minmass)
max.size <- log(maxmass)
n <- n_mesh
b <- seq(min.size, max.size, length.out = n + 1)
y <- 0.5 * (b[1:n] + b[2:(n + 1)])
h <- y[2] - y[1]

source("Code/MIPM_functions.R")   # the bug-fixed version

burnin2remove <- 1:(burn_yrs * 12)

create_fixed_temp_df <- function(fixdtemp, n_yrs = n_years) {
  meantemp <- rep(fixdtemp, 365 * n_yrs)
  df <- as.data.frame(meantemp)
  df$no.years.for.sim <- rep(1:n_yrs, each = 365)
  df
}

# ---- 2. Olivia's pipeline, wrapped: build one annual projection matrix ----
build_annual_matrix <- function(fxd_Temp, at0, z0) {
  fixed_temp_data <- create_fixed_temp_df(fxd_Temp)
  allmonths <- get_month_matrices(y = y, h = h, n = n,
                                   Temp_vector_data = fixed_temp_data,
                                   spawningmonth = 11, spawningday = 1,
                                   at0 = at0, z0 = z0)
  allmonths_noburnin <- allmonths[-burnin2remove]
  Reduce("%*%", allmonths_noburnin)
}
lambda_of <- function(A) Re(eigen(A, only.values = TRUE)$values[1])

# ---- 3. Baseline: Olivia's brute-force lambda at the MT parameter estimates ----
cat("[1/6] Building baseline matrix...\n")
t0 <- Sys.time()
A0 <- build_annual_matrix(10, at0 = MT_at0, z0 = MT_z0)
lambda0 <- lambda_of(A0)
cat("      done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    "sec. Baseline lambda =", lambda0, "\n\n")

# ---- 4. Brute-force sweep: lambda vs at0 (5 points) ----
cat("[2/6] Sweeping at0 (5 points, brute force)...\n")
at0_seq <- MT_at0 * seq(0.9, 1.1, length.out = 5)
lambda_seq_at0 <- numeric(5)
for (i in seq_along(at0_seq)) {
  t0 <- Sys.time()
  lambda_seq_at0[i] <- lambda_of(build_annual_matrix(10, at0 = at0_seq[i], z0 = MT_z0))
  cat("      point", i, "/5 done (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec)\n")
}
cat("\n")

# ---- 5. Brute-force sweep: lambda vs z0 (5 points) ----
cat("[3/6] Sweeping z0 (5 points, brute force)...\n")
z0_seq <- MT_z0 * seq(0.9, 1.1, length.out = 5)
lambda_seq_z0 <- numeric(5)
for (i in seq_along(z0_seq)) {
  t0 <- Sys.time()
  lambda_seq_z0[i] <- lambda_of(build_annual_matrix(10, at0 = MT_at0, z0 = z0_seq[i]))
  cat("      point", i, "/5 done (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec)\n")
}
cat("\n")

# ---- 6. Kato sensitivity matrix from baseline eigenvectors (Griffith Eq.4) ----
cat("[4/6] Computing Kato sensitivity matrix from baseline eigenvectors...\n")
kato_eigen <- function(A) {
  eR <- eigen(A); eL <- eigen(t(A))
  idxR <- which.max(Re(eR$values)); idxL <- which.max(Re(eL$values))
  lambda <- Re(eR$values[idxR])
  w <- Re(eR$vectors[, idxR])   # stable size distribution
  v <- Re(eL$vectors[, idxL])   # reproductive value
  if (sum(v * w) < 0) v <- -v
  list(lambda = lambda, v = v, w = w)
}
ke <- kato_eigen(A0)
SV <- outer(ke$v, ke$w) / sum(ke$v * ke$w)   # Kato / Griffith Eq.4
cat("      done.\n\n")

# ---- 7. Kato prediction vs brute force: d(lambda)/d(at0) ----
cat("[5/6] d(lambda)/d(at0): Kato prediction vs brute force...\n")
delta_at0 <- MT_at0 * 1e-3
A_plus  <- build_annual_matrix(10, at0 = MT_at0 + delta_at0, z0 = MT_z0)
A_minus <- build_annual_matrix(10, at0 = MT_at0 - delta_at0, z0 = MT_z0)
dA_dat0 <- (A_plus - A_minus) / (2 * delta_at0)
kato_dlam_dat0       <- sum(SV * dA_dat0)
bruteforce_dlam_dat0 <- (lambda_of(A_plus) - lambda_of(A_minus)) / (2 * delta_at0)
relerr_at0 <- abs(kato_dlam_dat0 - bruteforce_dlam_dat0) / abs(bruteforce_dlam_dat0)
cat(sprintf("      Kato prediction : %.6f\n", kato_dlam_dat0))
cat(sprintf("      Brute force     : %.6f\n", bruteforce_dlam_dat0))
cat(sprintf("      Relative error  : %.3e\n\n", relerr_at0))

# ---- 8. Kato prediction vs brute force: d(lambda)/d(z0) ----
cat("[6/6] d(lambda)/d(z0): Kato prediction vs brute force...\n")
delta_z0 <- MT_z0 * 1e-3
B_plus  <- build_annual_matrix(10, at0 = MT_at0, z0 = MT_z0 + delta_z0)
B_minus <- build_annual_matrix(10, at0 = MT_at0, z0 = MT_z0 - delta_z0)
dA_dz0 <- (B_plus - B_minus) / (2 * delta_z0)
kato_dlam_dz0       <- sum(SV * dA_dz0)
bruteforce_dlam_dz0 <- (lambda_of(B_plus) - lambda_of(B_minus)) / (2 * delta_z0)
relerr_z0 <- abs(kato_dlam_dz0 - bruteforce_dlam_dz0) / abs(bruteforce_dlam_dz0)
cat(sprintf("      Kato prediction : %.6f\n", kato_dlam_dz0))
cat(sprintf("      Brute force     : %.6f\n", bruteforce_dlam_dz0))
cat(sprintf("      Relative error  : %.3e\n\n", relerr_z0))

# ---- 9. PLOTS ----
# Each plot is defined as a function, then:
#   (a) called once directly -- this is what will appear live in RStudio's
#       Plots pane when you run this script there, and
#   (b) called again wrapped in png(...)/dev.off() -- this saves an
#       identical copy as a file you can attach to slides/email.

draw_plot1 <- function() {
  plot(at0_seq, lambda_seq_at0, type = "b", pch = 19, col = "steelblue", lwd = 2,
       xlab = "at0", ylab = "lambda",
       main = "Brute-force lambda vs at0\nwith Kato-predicted local slope")
  abline(a = lambda0 - kato_dlam_dat0 * MT_at0, b = kato_dlam_dat0,
         col = "firebrick", lwd = 2, lty = 2)
  points(MT_at0, lambda0, pch = 8, col = "firebrick", cex = 1.5)
  legend("topleft", legend = c("Brute-force lambda", "Kato-predicted tangent"),
         col = c("steelblue", "firebrick"), lty = c(1, 2), pch = c(19, NA), bty = "n")
}
draw_plot1()
png("plot1_lambda_vs_at0.png", width = 1000, height = 750, res = 130); draw_plot1(); dev.off()

draw_plot2 <- function() {
  plot(z0_seq, lambda_seq_z0, type = "b", pch = 19, col = "steelblue", lwd = 2,
       xlab = "z0", ylab = "lambda",
       main = "Brute-force lambda vs z0\nwith Kato-predicted local slope")
  abline(a = lambda0 - kato_dlam_dz0 * MT_z0, b = kato_dlam_dz0,
         col = "firebrick", lwd = 2, lty = 2)
  points(MT_z0, lambda0, pch = 8, col = "firebrick", cex = 1.5)
  legend("topright", legend = c("Brute-force lambda", "Kato-predicted tangent"),
         col = c("steelblue", "firebrick"), lty = c(1, 2), pch = c(19, NA), bty = "n")
}
draw_plot2()
png("plot2_lambda_vs_z0.png", width = 1000, height = 750, res = 130); draw_plot2(); dev.off()

draw_plot3 <- function() {
  barplot(c(relerr_at0, relerr_z0), names.arg = c("d(lambda)/d(at0)", "d(lambda)/d(z0)"),
          log = "y", col = "steelblue",
          ylab = "relative error (log scale, Kato vs brute force)",
          main = "Kato validation: relative error")
}
draw_plot3()
png("plot3_relative_error.png", width = 900, height = 650, res = 130); draw_plot3(); dev.off()

# ---- 10. Final summary printed to the Console ----
cat("\n========== SUMMARY (read this out / screenshot for the meeting) ==========\n")
cat(sprintf("Baseline lambda                 : %.6f\n", lambda0))
cat(sprintf("d(lambda)/d(at0)   Kato         : %.6f\n", kato_dlam_dat0))
cat(sprintf("d(lambda)/d(at0)   Brute force  : %.6f\n", bruteforce_dlam_dat0))
cat(sprintf("Relative error                   : %.3e\n", relerr_at0))
cat(sprintf("d(lambda)/d(z0)    Kato         : %.6f\n", kato_dlam_dz0))
cat(sprintf("d(lambda)/d(z0)    Brute force  : %.6f\n", bruteforce_dlam_dz0))
cat(sprintf("Relative error                   : %.3e\n", relerr_z0))
cat("Plots saved in your Files pane: plot1_lambda_vs_at0.png, plot2_lambda_vs_z0.png, plot3_relative_error.png\n")
