# Nonignorable Missing Data

To start the simulation, run `src/simulation-study/main.r`. Parallelization is based on process forking, so it is only available on UNIX-like systems. To run this code on a Windows machine, use sequential execution: replace all instances of `mclapply()` with `lapply()`. 

The R code in `src/simulation-study/tsonaka_etal/` was written and provided by Roula Tsonaka.