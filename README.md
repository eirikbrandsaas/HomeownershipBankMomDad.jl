# HomeownershipBankMomDad
[![Build Status](https://github.com/eirikbrandsaas/HomeownershipBankMomDad.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/eirikbrandsaas/HomeownershipBankMomDad.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/eirikbrandsaas/HomeownershipBankMomDad.jl/graph/badge.svg?token=4IVYUDK98Q)](https://codecov.io/gh/eirikbrandsaas/HomeownershipBankMomDad.jl)

This repo contains all material for the paper "Illiquid Homeownership and the Bank of Mom and Dad".
## Replication Instructions
1. Fork the repo to your local machine.
    1. Change the global `basepath` in `src/Stata/runall.do` to refer the directory of this file
    2. The directory called `HomeownershipBankMomDad.jl/` will have the correct folder structure and all you will have to do is manually is to download PSID data and run one Stata program and one Julia program to recreate all results in the paper.
2. The file `docs/paper.tex` is compilable immediately (pdflatex paper.tex -> bibtex paper -> pdflatex paper.tex -> pdflatex paper.tex)
    * Each table and figure is tracked with source control. This makes it very easy to see the impact of any changes to the code or data on results.
3. Download the [Panel Study of Income Dynamics (PSID) packaged data](https://simba.isr.umich.edu/Zips/ZipMain.aspx). Note that you have to manually change the paths in the PSID-provided transfer supplement, see below.
    1. You will need 
       1. `famYYYY.zip` from 1999-2017, 
       2. `wlthYYYY.zip` from 1999-2007,
       3. `ind2019er.zip`,
       4. `pid19.zip` (Parental Identification) and,
       5. `RT13.zip` (2013 Family Roster and Transfers)
    2. Place the *zipped* folders into `data/PSID`
    3. Change the `[paths]` of the two do files in `RT13/` to `data/PSID/RT13`:
       1. In `RT13FAM.do` change `using [path]\RT13FAM.txt, clear ` to `using data/PSID/RT13/RT13FAM.txt, clear`, line 31.
       2. In `RT13PARCHD.do` change `using [path]\RT13PARCHD.txt, clear` to `using data/PSID/RT13/RT13PARCHD.txt, clear`, line 38.
4. Run Stata code by running `src/Stata/runall.do`
   1. Automatically installs necessary packages
   2. Defines necessary globals
   3. Downloads some data from FRED. (If you havent used FRED through Stata before you must set up acces key. Follow on-screen instructions)
   4. Merges and cleans PSID data using the `PSIDtools` package. Will be slow (~2-5 minutes) first time when it unzips etc.
   5. Runs a series of files that creates tables, plots, estimates calibrated parameters, and calculates moments used in the structural estimation
   6. Finally it downloads SHED data and plots SHED and AHS plots used in the introduction.
      1. `SHED.do` will automatically download and manipulate the 2013-2017 "Stata Files (ZIP)" available at [SHED's Board website](https://www.federalreserve.gov/consumerscommunities/shed_data.htm)
      2. `AHS.do` loads the data in `data/AHS/`, which consist of CSV from the Table generator (National->Year->Value, Purchase Price, and Source of Downpayment -> get table) for 2011-2019 and manually inputted data for 1991-2009, which I took from the AHS Summary Table (with `.pdf`s in the directory as well). All data is from Table 3-14 "Major Source of Down Payment". 
      3. `SHED_AHS_plots.do` plots the figure in the introduction.
5. Run Julia code:
   1. Start Julia in `--project` mode:
      1. If using VSCode, I recommend that you start VSCode in the directory by launching `code PathTo/HomeownershipBankMomDad.jl/` which will automatically launch Julia with the right settings
      2. If using the terminal/Julia REPL, navigate to the directory in the terminal and launch Julia with `julia --project`
   2. Run the file `run_all.jl`, which generates all tables and figures the estimation procedure. The program is divided in two where the first repeats the structural estimation. The program is written so that you can skip the estimation (which is very slow) and recreate all plots.
      1. The first section runs the structural estimation. This can be skipped. **If** you chose to re-estimate the model note that the other code does not manually update the parameters.
         1. To use new parameters, change the inputted parameters in `function benchpar()` in `run_all_utils.jl` manually
      2. If you don't wish to re-estimate, you can run the second part of the program.
