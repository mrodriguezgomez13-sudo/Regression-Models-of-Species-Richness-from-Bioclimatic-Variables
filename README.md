# Statistical Modelling of Global Species Richness from Bioclimatic Variables

**Author:** Mònica Rodríguez Gómez  
**Dataset:** WorldClim bioclimatic variables (n = 120)  
**Tools:** R, ggplot2, GGally, car, MASS, leaps, pls, glmnet

---

## Overview

This project investigates the relationship between global species richness (`richness`) and six bioclimatic predictors from the WorldClim database. The analysis combines exploratory data analysis, multiple linear regression, model diagnostics, multicollinearity assessment, variable selection, and predictive modelling with PCR, PLS, Ridge, and LASSO regression.

The goal is to identify which climatic variables are the main determinants of species richness and to compare the predictive performance of different regularization and dimension-reduction techniques.


---

## Data

The dataset contains **120 observations** and **8 variables**:

- `richness`: number of species (response variable)
- `bio1`: annual mean temperature
- `bio4`: temperature seasonality
- `bio5`: max temperature of the warmest month
- `bio6`: min temperature of the coldest month
- `bio12`: annual precipitation
- `bio15`: precipitation seasonality
- `biome`: categorical variable (Grassland, Forest, Desert, etc.)

No missing values were found. The `biome` variable was converted to a factor for exploratory analysis.

---

## Methods

### 1. Exploratory Data Analysis
- Distribution of the response variable
- Scatterplot matrix by biome
- Correlation analysis

### 2. Multiple Linear Regression
- Full model with six bioclimatic predictors
- Multicollinearity assessment via VIF and condition number (κ)
- Model diagnostics: residuals, Q-Q plot, Shapiro-Wilk test, Cook's distance, component + residual plots

### 3. Variable Selection
- Stepwise AIC selection
- Best-subset selection based on adjusted R²
- Nested model comparison via ANOVA F-test

### 4. Predictive Modelling
- Train/test split (80/20) with standardization based on training data
- Principal Component Regression (PCR)
- Partial Least Squares (PLS)
- Ridge regression
- LASSO regression
- Model comparison using test-set RMSE

---

## Key Results

### Multicollinearity
- `bio1` (VIF = 18.75) and `bio5` (VIF = 16.08) show **severe multicollinearity**.
- `bio6` (VIF = 6.54) shows high multicollinearity.
- Condition number **κ = 147.48**, confirming a severe near-linear dependency among predictors.

### Variable Selection
Both stepwise AIC and best-subset selection converge on the same model:

**`richness ~ bio1 + bio4 + bio12`**

A nested F-test confirms that `bio5`, `bio6`, and `bio15` do not add significant information (p = 0.82).

### Model Diagnostics
- Residuals show slight deviations from normality (Shapiro-Wilk p = 0.00024).
- No observation exceeds Cook's distance of 1, so no influential outliers are problematic.
- `bio1`, `bio5`, and `bio6` show slight curvature in component + residual plots, suggesting possible transformations.

### Predictive Performance

| Model | RMSE |
|-------|------|
| PCR   | 9.553 |
| PLS   | 9.539 |
| Ridge | 9.369 |
| LASSO | 9.539 |

- **Ridge** achieves the lowest test RMSE (9.369).
- PLS and LASSO perform slightly better than PCR.
- Differences between methods are small, indicating that the predictive information is stable across techniques.

---

## Conclusions

- Annual mean temperature (`bio1`), temperature seasonality (`bio4`), and annual precipitation (`bio12`) are the main determinants of global species richness in this dataset.
- Severe multicollinearity among temperature variables inflates standard errors and destabilizes individual coefficient estimates.
- Ridge regression provides the best predictive performance by stabilizing coefficients without eliminating variables.
- The analysis highlights the importance of combining model diagnostics, variable selection, and regularized regression when predictors are highly collinear.

---

## Requirements

- R (>= 4.0)
- Packages: `ggplot2`, `GGally`, `car`, `MASS`, `leaps`, `pls`, `glmnet`

```r
install.packages(c(
  "ggplot2", "GGally", "car", "MASS",
  "leaps", "pls", "glmnet"
))
