# ROI package

## neos solver

 A **valid email address** is required for all NEOS jobs.

```R
job_neos <- ROI_solve(op, "neos", email = "putyouremailhere")
```

## gurobi solver

- Step 1: Register a Gurobi account as an academic user at https://pages.gurobi.com/registration.
- Step 2: Request for a license at  https://www.gurobi.com/downloads/end-user-license-agreement-academic/ and you will obtain a Gurobi key. Copy and paste the Gurobi key to the terminal.
  - **Warning: One year free license. One license can only support one computer to download this optimizer.**

- Step 3: Go to https://www.gurobi.com/downloads/gurobi-optimizer-eula/ , read and accept the End User License Agreement. Then download the current version of Gurobi optimizer. Install the optimizer with default settings.
- In RStudio, install the `slam` package using `install.packages('slam')` and gurobi plugin using `remotes:::install_github("fl0sch/ROI.plugin.gurobi")`.

## mosek solver

Not applicable to our OP.

- Step 1: Install **MOSEK**.

  - Download **MOSEK** Optimizer from https://mosek.com/downloads/.

  - Run the command

    ```
    python3 <MSKHOME>/mosek/10.0/tools/platform/<PLATFORM>/bin/install.py
    ```
     <MSKHOME> is the directory where MOSEK was installed and <PLATFORM> is osx64x86 or osxaarch64 depending on the version of MOSEK installed.
    
  - Add the path `<MSKHOME>/mosek/10.0/tools/platform/<PLATFORM>/bin` to the OS variable `PATH`.

    ```markdown
    - Open up Terminal.
    - Run the command `sudo nano /etc/paths`.
    - Go to the bottom of the file, and enter the path `<MSKHOME>/mosek/10.0/tools/platform/<PLATFORM>/bin`.
    - Save it and exit.
    ```

- Step2: Get the Academic License (valid for 365 days, but can be renewed multiple times) from https://www.mosek.com/products/academic-licenses/ and set up the license. It should be saved to a file called `/Users/YOUR_USER_NAME/mosek/mosek.lic` (OSX).

- Step 3: Make sure `xcode` installed and install the `Rmosek` R package.

  ```R
  source("<RMOSEKDIR>/builder.R") # <RMOSEKDIR> is /Users/xwan0362/mosek/10.0/tools/platform/osxaarch64/rmosek
  attachbuilder()
  install.rmosek()
  ```

- Step 4: Test the Installation using `require("Rmosek")`. Open a terminal, change folder to `/Users/xwan0362/mosek/10.0/tools/examples/rmosek` and use R to run a selected example, for instance `R -f gp1.R`.

- Step 5: Install mosek plugin in R using `remotes:::install_github("fl0sch/ROI.plugin.mosek")`.

## cplex solver

**Plugin installed. It can handle simple examples given by the ROI package but will terminate the R session when dealing with our OP...**

- Step 1: install ILOG CPLEX Optimization Studio (Academic edition VS No-cost edition).
  - The **no-cost edition** is restricted to problems up to 1,000 variables and 1,000 constraints. https://www.ibm.com/account/reg/us-en/signup?formid=urx-20028

  - There is **an unlimited version** of IBM ILOG CPLEX Optimization Studio for students and faculty members. https://www.ibm.com/academic/home However, the Academic Initiative/SkillsBuild website has experienced technical difficulties for several days...

- Step 2: Find "Makefile" using the path "/Applications/CPLEX_Studio_Community2211/cplex/examples/arm64_osx/static_pic/Makefile" (Refer to https://github.com/cran/Rcplex/blob/master/inst/INSTALL).

- Step 3: Makefile might have

```
CLNFLAGS = -l{CPLEXLIB} -m64 -lm -lpthread -framework CoreFoundation -framework IOKit
CFLAGS = $(COPT)  -I${CPLEXINCDIR}
```

Set/replace the `${CPLEXINCDIR}` and `${CPLEXLIB}` with the corresponding paths to the 'include' and 'lib' directories of your CPLEX installation on your system.

So, we have

```
CLNFLAGS  = -l/Applications/CPLEX_Studio_Community2211/cplex/lib/arm64_osx/static_pic -m64 -lm -lpthread -framework CoreFoundation -framework IOKit

CFLAGS  = $(COPT)  -I/Applications/CPLEX_Studio_Community2211/cplex/include
```

- Step 4: Download Rcplex_0.3-5.tar.gz from R CRAN and open the terminal, run

```
R CMD INSTALL --configure-args="PKG_CFLAGS='-m64 -fPIC' \
PKG_CPPFLAGS=-I/Applications/CPLEX_Studio_Community2211/cplex/include \
PKG_LIBS='-L/Applications/CPLEX_Studio_Community2211/cplex/lib/arm64_osx/static_pic
-lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit'" Rcplex_0.3-5.tar.gz
```

Then, the Rcplex is installed!

- Step 5: In RStudio, install cplex plugin using `install.packages("ROI.plugin.cplex")`.



# MiniZinc

## MiniZinc 

### Introduction

- MiniZinc --> Compiler
- FlatZinc --> A solver input language that is understood by a wide range of solvers
- MiniZin IDE --> Intergrated Development Enviroment
- Solvers --> Gurobi, CPLEX, Gecode ...
- Bundled binary packages --> They contain the compiler and IDE, as well as the following solvers: Gecode, Chuffed, COIN-OR CBC, and a interfaces to the Gurobi and CPLEX solvers (the Gurobi and CPLEX libraries are not included).

### Installation

- Download and install **bundled binary packages**, available from http://www.minizinc.org/software.html.

- In order to use the MiniZinc tools from a terminal, add the path to the MiniZinc installation to the PATH environment variable. For macOS:

  ```
  $ export PATH=/Applications/MiniZincIDE.app/Contents/Resources:$PATH
  ```

- Configuring existing solvers. For example, [Section 3.2.5.2 of the MiniZinc Handbook](https://www.minizinc.org/doc-2.5.5/en/minizinc_ide.html#sec-ide-add-solvers) shows the configuratio for Gurobi and CPLEX solvers.

### Problems

1. Some solver interfaces to MiniZinc currently don't support **quadratic constraints**. We can multiply decision variables by constants, and can add these terms together, but we cannot multiply two decision variables. Gurobi supports the variables multiplication and we need to define `QuadrFloat=true` to specify that the solver supports quadratic constraints. But:

2. OP takes too much time. `float` and `int` variables should have as tight domains as possible to improve solving. Some solvers don't like **unbounded variables** at all and might hang or give an error if one is encountered.

3. Matrix multiplication may not be supported. No examples in the MiniZinc handbook for MIP, quadratic constraints, objective functions or constraints formed with matrix.



## R interface to MiniZinc

**installed but didn't work...**

This package is compatible with MiniZinc 2.5.5 and all the features haven't been tested for all of MiniZinc.

- Step 1: Build libminizinc
  - Download the libminizinc 2.5.5 release from https://github.com/MiniZinc/libminizinc/releases.
  - Extract the downloaded tar.gz or zip (for Windows) and name the extracted folder `libminizinc`.
  - `cd libminizinc/`
  -  Add `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` in the 3rd line of the file `CMakeLists.txt`.
  - `sudo cmake CMakeLists.txt`
  - `sudo make`
  - `sudo make install`
- Step 2: Get Solver Binaries
  - Download the MiniZinc bundled binary package from https://www.minizinc.org/software.html (MiniZincIDE).
  - MiniZinc binary solvers are save in the folder `/Applications/MiniZincIDE.app/Contents/Resources/bin`.
  - Copy the solver binaries to another folder and just provide the path to that folder with `--with-bin`.

- Step 3: Install rminizinc using `install.packages("rminizinc", configure.args="--with-mzn=/path/to/libminizinc --with-bin=/path/to/bin")`.

