# ============================================================
# Statistical Modelling of Global Species Richness
# from Bioclimatic Variables
# ------------------------------------------------------------
# Author: Mònica Rodríguez Gómez
# Dataset: WorldClim bioclimatic variables (n = 120)
# ============================================================

# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

# Load required packages
library(ggplot2)
library(GGally)
library(car)
library(MASS)
library(leaps)
library(pls)
library(glmnet)

# Set seed for reproducibility
set.seed(123)

# ------------------------------------------------------------
# 1. Data Loading and Exploration
# ------------------------------------------------------------

# Load the dataset
datos <- read.csv("data/dataset_WorldClim.csv", header = TRUE)

# Inspect the data
head(datos)
dim(datos)        # 120 observations, 8 variables
str(datos)
sum(is.na(datos)) # 0 missing values

# Convert biome to factor
datos$biome <- as.factor(datos$biome)

# ------------------------------------------------------------
# 2. Exploratory Data Analysis
# ------------------------------------------------------------

# 2.1 Distribution of the response variable
ggplot(datos, aes(x = richness)) +
  geom_histogram(fill = "lightblue", color = "white") +
  labs(
    title = "Distribution of Species Richness",
    x = "Number of species", y = "Frequency"
  ) +
  theme_minimal()

# 2.2 Scatterplot matrix by biome
p <- ggpairs(
  datos, columns = 1:7,
  aes(color = biome, alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 2)),
  lower = list(continuous = wrap("points", size = 0.7)),
  diag  = list(continuous = wrap("densityDiag"))
) +
  labs(title = "Scatterplot Matrix by Biome") +
  theme_bw()

p + theme(axis.text = element_text(size = 6))

# ------------------------------------------------------------
# 3. Multiple Linear Regression
# ------------------------------------------------------------

# Full model with six bioclimatic predictors
m_1 <- lm(richness ~ bio1 + bio4 + bio5 + bio6 + bio12 + bio15, data = datos)
summary(m_1)

# 3.1 Multicollinearity assessment
vif(m_1)

kappa_m1 <- model.matrix(m_1)[, -1]
kappa(kappa_m1, exact = TRUE)

# ------------------------------------------------------------
# 4. Model Diagnostics
# ------------------------------------------------------------

# 4.1 Diagnostic plots
par(mfrow = c(2, 2))
plot(m_1)

# 4.2 Normality of residuals
shapiro.test(residuals(m_1))

# 4.3 Component + residual plots
crPlots(m_1)

# 4.4 Influential observations (Cook's distance)
cooksd <- cooks.distance(m_1)
which(cooksd > 4 / (nrow(datos) - 6 - 1))
which(cooksd > 1)

# ------------------------------------------------------------
# 5. Variable Selection
# ------------------------------------------------------------

# 5.1 Stepwise AIC selection
m1_aic <- stepAIC(m_1, direction = "both")
formula(m1_aic)

# 5.2 Best-subset selection based on adjusted R²
seleccion <- regsubsets(
  richness ~ bio1 + bio4 + bio5 + bio6 + bio12 + bio15,
  data = datos, nbest = 1, nvmax = 6
)
resumen_sel <- summary(seleccion)

which.max(resumen_sel$adjr2)
coef(seleccion, which.max(resumen_sel$adjr2))

plot(
  resumen_sel$adjr2, type = "b",
  xlab = "Number of predictors",
  ylab = expression(R^2 - adjusted),
  main = "Subset Selection"
)

# 5.3 Nested model comparison (F-test)
m_1_reducido <- lm(richness ~ bio1 + bio4 + bio12, data = datos)
anova(m_1_reducido, m_1)
summary(m_1_reducido)

# ------------------------------------------------------------
# 6. Predictive Modelling
# ------------------------------------------------------------

# 6.1 Train/test split
test_obs <- seq(1, nrow(datos), by = 5)
test  <- datos[test_obs, ]
train <- datos[-test_obs, ]

x_train <- train[, c("bio1", "bio4", "bio5", "bio6", "bio12", "bio15")]
y_train <- train$richness

x_test <- test[, c("bio1", "bio4", "bio5", "bio6", "bio12", "bio15")]
y_test <- test$richness

# Standardize using training-set statistics
media <- apply(x_train, 2, mean)
desv  <- apply(x_train, 2, sd)

x_train_sc <- scale(x_train, center = media, scale = desv)
x_test_sc  <- scale(x_test,  center = media, scale = desv)

# 6.2 Principal Component Regression (PCR)
m1_pcr <- pcr(y_train ~ x_train_sc, scale = FALSE, validation = "CV")
summary(m1_pcr)
validationplot(m1_pcr, val.type = "MSEP")

pred_pcr_m1 <- predict(m1_pcr, newdata = x_test_sc, ncomp = 4)
rmse_pcr_m1 <- sqrt(mean((y_test - pred_pcr_m1)^2))
rmse_pcr_m1

# 6.3 Partial Least Squares (PLS)
pls_m1 <- plsr(y_train ~ x_train_sc, scale = FALSE, validation = "CV")
summary(pls_m1)

pred_pls <- predict(pls_m1, newdata = x_test_sc, ncomp = 3)
rmse_pls <- sqrt(mean((y_test - pred_pls)^2))
rmse_pls

# 6.4 Ridge regression
x_train_mat <- as.matrix(x_train_sc)
x_test_mat  <- as.matrix(x_test_sc)

set.seed(123)
ridge_m1 <- cv.glmnet(x_train_mat, y_train, alpha = 0, standardize = TRUE)

lambda_m1_ridge <- ridge_m1$lambda.min
lambda_m1_ridge

ridge_lambd <- glmnet(
  x_train_mat, y_train,
  alpha = 0, lambda = lambda_m1_ridge
)
coef(ridge_lambd)

plot(glmnet(x_train_mat, y_train, alpha = 0))
abline(v = -log(lambda_m1_ridge), col = "red", lty = 2)

pred_ridge <- predict(ridge_m1, newx = x_test_mat, s = "lambda.min")
rmse_ridge <- sqrt(mean((y_test - pred_ridge)^2))
rmse_ridge

# 6.5 LASSO regression
lasso_m1 <- cv.glmnet(x_train_mat, y_train, alpha = 1, standardize = TRUE)
lambda_m1_lasso <- lasso_m1$lambda.min
lambda_m1_lasso

m1_lasso <- glmnet(
  x_train_mat, y_train,
  alpha = 1, lambda = lambda_m1_lasso
)

pred_lasso <- predict(lasso_m1, newx = x_test_mat, s = "lambda.min")
rmse_lasso <- sqrt(mean((y_test - pred_lasso)^2))
rmse_lasso

# 6.6 Coefficient comparison (Ridge vs LASSO)
coef_ridge <- coef(ridge_lambd)
coef_lasso <- coef(m1_lasso)

tabla_coef <- data.frame(
  Variable = rownames(coef_ridge),
  Ridge = round(as.vector(coef_ridge), 4),
  LASSO = round(as.vector(coef_lasso), 4)
)
tabla_coef

# 6.7 Model comparison
tabla_rmse <- data.frame(
  Model = c("PCR", "PLS", "Ridge", "LASSO"),
  RMSE  = c(rmse_pcr_m1, rmse_pls, rmse_ridge, rmse_lasso)
)
tabla_rmse

# ------------------------------------------------------------
# 7. Session Information
# ------------------------------------------------------------

sessionInfo()