## neos solver

 A valid email address is required for all NEOS jobs.

```R
job_neos <- ROI_solve(op, "neos", email = "putyouremailhere")
```

## gurobi solver

Warning: One year free license.

- Step 1: Register a Gurobi account as an academic user at https://pages.gurobi.com/registration.
- Step 2: Request for a license at  https://www.gurobi.com/downloads/end-user-license-agreement-academic/ and you will obtain a Gurobi key. Copy and paste the Gurobi key to the terminal.
- Step 3: Go to https://www.gurobi.com/downloads/gurobi-optimizer-eula/ , read and accept the End User License Agreement. Then download the current version of Gurobi optimizer. Install the optimizer with default settings.
- In RStudio, install the `slam` package using `install.packages('slam')` and gurobi plugin using `remotes:::install_github("fl0sch/ROI.plugin.gurobi")`.

## cplex solver

Plugin installed but didn't work...

Warning: This version allows models with up to 1,000 variables and 1,000 constraints to be run.

- Step 1: install ILOG CPLEX Optimization Studio (Register first, download link: https://ibm.ent.box.com/s/wjuh81fmorssmwwoe4eg2spfl9xrakvn).

- Step 2: Find "Makefile" using the path "/Applications/CPLEX_Studio_Community2211/cplex/examples/arm64_osx/static_pic/Makefile" (Refer to https://github.com/cran/Rcplex/blob/master/inst/INSTALL).

- Step 3: Makefile might have

```
CLNFLAGS = -l{CPLEX_LIB_PATH}/static_pic -m64 -lm -lpthread -framework CoreFoundation -framework IOKit
CFLAGS = $(COPT)  -I${CPLEX_INCLUDE_PATH}
```

Set/replace the ${CPLEX_INCLUDE_PATH} and ${CPLEX_LIB_PATH} with the corresponding paths to the 'include' and 'lib' directories of your CPLEX installation on your system.

So, we have

```
CLNFLAGS  = -l/Applications/CPLEX_Studio_Community2211/cplex/lib/arm64_osx/static_pic -m64 -lm -lpthread -framework CoreFoundation -framework IOKit

CFLAGS  = $(COPT)  -I/Applications/CPLEX_Studio_Community2211/cplex/include
```

- Step 4: Open the terminal, and run

```
R CMD INSTALL --configure-args="PKG_CFLAGS='-m64 -fPIC'
PKG_CPPFLAGS=-I/Applications/CPLEX_Studio_Community2211/cplex/include
PKG_LIBS='-L/Applications/CPLEX_Studio_Community2211/cplex/lib/arm64_osx/static_pic
-lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit'" Rcplex
```

Then, the Rcplex is installed!

- Step 5: In RStudio, install cplex plugin using `install.packages("ROI.plugin.cplex")`.