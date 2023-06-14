library(reticulate)

# Install an ARM build of conda and use it
# reticulate::install_miniconda()

# Install packages
# py_install("numpy", pip = TRUE); py_install("gurobipy", pip = TRUE); py_install("scipy", pip = TRUE)

# Set the path to the Python executable file
reticulate::use_python("/Users/xwan0362/Library/r-miniconda-arm64/bin/python3.10", required = T)

# Check the version of Python
py_config()

# Evaluate the chosen script
reticulate::source_python("Python/subset.py")

y <- c(10, 6, 4, 2, 4, 1, 6) |> matrix(nrow = 7)
S <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1), c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
W <- diag(nrow=7)

fit <- miqp(y = fc, S = S, W = W, l0 = 0.5)


