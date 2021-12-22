# Sufficient Dimension Reduction for Mixed Predictors

## Description

### This repository contains code to reproduce results in E. Bura, L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, "Sufficient Reductions in Regression with Mixed Predictors" (submitted).

The main code for the proposed method are written in Matlab. However, some comparisons with other methods require running scripts in R. Indications to reproduce results reported in the manuscript are detailed below.

## Organization

The code is organized in several folders:

```Data-Analysis/```: Contains the scripts to reproduce the examples presented in Section 6: Data Analyses.

```Figures/```: Contains the scripts to reproduce Figures 1 and 2 presented in Section 5: Simulation Studies.

```Internal/```: Contains the procedures and auxiliary unctions to actually implement the proposed methods.

```Main-Functions/```: Contains the functions which are called to apply the proposed methods.

```Other-Used-Packages/```: Contains other tools used to implement and run the code.

```Simulations/```: Contains the scripts to reproduce the results reported in Section 5: Simulation Studies.

```Test-Dimension/```: Contains the scripts to reproduce results in Table 3 from Section 5.


## Usage

To start using the code, set Matlab's working directory to the main folder of this package and run

```
> setpaths.m
```

This will add paths to all the internal folders so that all the functions and datasets become available.


## Reproducing experiments with synthetic data in Section 5: Simulation Studies

To repreoduce Figures 1 and 2, run the following scripts from ```Figures/```:

```
> script-Figure-1.R
> script-Figure-2.R
```

If you want to reproduce the results used in Figure 1, Figure 2, and Table 2, do the following:

In ```Simulations/ContinuousPredictors/```, run:
```
> script-to-simulate-ContinuousPredictors-d1.m
> script-to-simulate-ContinuousPredictors-d2.m
```

In ```Simulations/BinaryPredictors/```, run:
```
> script-to-simulate-BinaryPredictors-d1.m
> script-to-simulate-BinaryPredictors-d2.m
```

In ```Simulations/MixedPredictors/```, run:
```
> script-to-simulate-MixedPredictors-d1.m
> script-to-simulate-MixedPredictors-d2.m
```

To reproduce results reported in Table 3, do the following:

In ```Test-Dimension/ContinuousPredictors/```, run:
```
> test-dimension-ContinuousPredictors-d1.m
> test-dimension-ContinuousPredictors-d2.m
```

In ```Test-Dimension/BinaryPredictors/```, run:
```
> test-dimension-BinaryPredictors-d1.m
> test-dimension-BinaryPredictors-d2.m
```

In ```Test-Dimension/MixedPredictors/```, run:
```
> test-dimension-MixedPredictors-d1.m
> test-dimension-MixedPredictors-d2.m
> test-dimension-Suboptimal.m
```

To reproduce results reported in Table 6, do the following:

In ```Simulations/ContinuousPredictors/```, run:
```
> script-to-simulate-ContinuousPredictors-d1-NonNormal.m
> script-to-simulate-ContinuousPredictors-d2-NonNormal.m
```

In ```Simulations/MixedPredictors/```, run:
```
> script-to-simulate-MixedPredictors-d1-NonNormal.m
> script-to-simulate-MixedPredictors-d2-NonNormal.m
```

## Reproducing experiments with real data

To reproduce results reported in Table 4, run the following scripts from ```Data-Analysis/Example1-Krzanowski-DataSets/```:

```
> run-krzanowski-PCA.R
> run-krzanowski-PFC.R
```

To reproduce Figures 3, 4 and 5, run the following script from ```Data-analysis/Example2-Governance-Index/```:

```
> script-for-Figures-3-4-5.m
```

If you want to obtain the Composite Governance indices from different methodologies, used in Figures 3, run the following scripts from ```Data-Analysis/Example2-Governance-Index/```:
```
> script-for-PCAmixindex.R
> script-for-PFC-Optimal-Suboptimal-Indices.m
```

To reproduce results reported in Table 5, run the following script from ```Data-Analysis/Example2-Governance-Index/```:

```
> script-for-Table5.R
```











