my_wilcoxon <- function(data_1, data_2) {
  m      <- length(data_1)
  n      <- length(data_2)
  group  <- c(rep(1, m), rep(2, n))
  pooled <- c(data_1, data_2)
  ranks  <- rank(pooled)         
  T1     <- sum(ranks[group == 1])
  return(T1)
}

# Data
Group1 <- c(1.9168, 3.5102, 3.4567, 2.7052, 1.3366,
            2.7654, 2.9780, 2.8394, 2.1455, 1.7992)

Group2 <- c(2.9204, 3.1500, 3.3405, 2.7374, 4.0499,
            3.8547, 4.1359, 4.0159, 3.4328, 2.0894,
            3.7789, 4.2853, 2.8724, 3.0208, 3.1968)

# Compute T1
T1 <- my_wilcoxon(Group1, Group2)
cat("T1 =", T1, "\n")
# T1 = 83

# Theoretical mean and variance under H0
m <- length(Group1)
n <- length(Group2)
ET1  <- m * (m + n + 1) / 2
VarT1 <- m * n * (m + n + 1) / 12
cat("E[T1]   =", ET1,   "\n")   # 130
cat("Var[T1] =", VarT1, "\n")   # 325

# Verify via wilcox.test (returns U = T1 - m(m+1)/2)
wt <- wilcox.test(Group1, Group2, alternative = "less")
cat("Mann-Whitney U =", wt$statistic, "\n")   # W = 28
cat("T1 from U      =", wt$statistic + m*(m+1)/2, "\n")  # 28 + 55 = 83

# One-tailed p-value for H_A: mu1 - mu2 < 0
cat("p-value =", wt$p.value, "\n")
