## using mice to impute the missing data
```{r, echo=FALSE}
set.seed(2025)
imputed_training_data <- mice(df_train, m=5, method='pmm', maxit=50)
imputed_training_data_c <- complete(imputed_training_data, 1)

imputed_test_data <- mice(df_test_dl, m=5, method='pmm', maxit=50)
imputed_test_data_c <- complete(imputed_test_data, 1)
```

```{r}
#F2F17, F2F4, F3F1  deleted
fs<-list(Comprehension=c("F1F1",  "F1F2",  "F1F3",  "F1F4",  "F1F5",  "F1F6",  "F1F7",  "F1F8",  "F1F9",  "F1F10", "F1F11",                              "F1F12", "F1F13", "F1F14"),
         Evaluation=c('F2F1',"F2F2",  "F2F3",  "F2F5",  "F2F6",  "F2F7",  "F2F8",  "F2F9",  "F2F10", "F2F11", "F2F12",                                "F2F13", "F2F14", "F2F15", "F2F16", "F2F18",'F2F19','F2F20'),
         Integration=c("F3F2",'F3F3' ,"F3F4",  "F3F5",  "F3F6",  "F3F7", "F3F8","F3F9"),
         Communication=c('F4F1','F4F2','F4F3','F4F4','F4F5','F4F6','F4F7','F4F8'),
         Statistics=c("F5F1",  "F5F2",  "F5F3",  "F5F4",  "F5F5",  "F5F6",  "F5F7", "F5F8",  "F5F9",  "F5F10",                                       "F5F11","F5F12", "F5F13", "F5F14", "F5F15", "F5F16", "F5F17", "F5F18")
)
```

## ordinal data transformation?
```{r}
ords <- imputed_training_data_merged3[, names(imputed_training_data_merged3)%in%unlist(fs)]
ords <- lapply(ords, as.ordered)
ords <- do.call(data.frame, ords)

# Print the resulting data frame
print(ords)
```

## Objective-Function
```{r}
library(stuart)
library(lavaan)


#ordinal data
objective.normal <- function(rmsea.scaled, srmr, cfi.scaled) {
  out1 = 0.5-(0.5/(1 + exp(- 100 * (rmsea.scaled-.05))))
  out2 = 0.5-(0.5/(1 + exp(- 100 * (srmr-.05))))
  out3 = (1/(1 + exp(- 100 * (cfi.scaled-.95))))
  out = (out1 + out2 + out3)/3
  return(out)                                
}



obj_default3 <- function(chisq, df, pvalue, rmsea.robust, srmr, cfi.robust, crel, nu)
{(1-rmsea.robust)+(1-srmr)+(1+cfi.robust)} # best jet


```

## Genetic Algorithm
### Kfold-Crossvalidation
```{r}
############ k-Folding #########################

kfold_sel <- kfold(
  'gene', 
  k = 3, 
  data = imputed_training_data_c,
  factor.structure = fs,  
  max.invariance = "strict", 
  capacity = list(4, 4, 4, 4, 4), 
  seed = 2028,
  seeded.search = TRUE,
  auxiliary = NULL,
  software = "lavaan",
  cores = 4,
  objective = obj_default3,
  ignore.errors = TRUE,
  burnin = 5,
  generations = 500,
  individuals = 34,
  selection = "tournament",
  selection.pressure = NULL,
  elitism = NULL,
  reproduction = 0.5,
  mutation = 0.05,
  mating.index = 0,
  mating.size = 0.25,
  mating.criterion = "similarity",
  immigration = 0,
  convergence.criterion = "geno.between",
  tolerance = NULL,
  reinit.n = 1,
  reinit.criterion = "geno.between",
  reinit.tolerance = NULL,
  reinit.prop = 0.75,
  schedule = "run",
  analysis.options = NULL,
  suppress.model = FALSE,
  remove.details = TRUE
)
#########################################################
```

### Summary Results No.1
```{r}
summary(kfold_sel)

inspect(kfold_sel$final, 'est')
inspect(kfold_sel$final, 'fit')

lavaan::summary(kfold_sel$final, standardized = TRUE) 
lavaan::inspect(kfold_sel$final, 'est')$A$nu
lavaan::inspect(kfold_sel$final, 'est')$B$nu
```
SUMMARY OF ANALYSIS:
  
  Number of Folds: 3 
Analysis Type: gene 
Estimation Software: lavaan 
Models Estimated: 78720 
Time Required: 8852.842 seconds

Crossvalidation Results with STRICT Invariance:
  
  Average Jaccard Similarity: Comprehension: 0.206, Evaluation: 0.000, Integration: 0.206, Communication: 0.511, Statistics: 0.270

Constructed Subtests: (k = 3)
Comprehension: F1F3 F1F4 F1F5 F1F14
Evaluation: F2F1 F2F7 F2F14 F2F20
Integration: F3F1 F3F4 F3F7 F3F8
Communication: F4F3 F4F6 F4F7 F4F8
Statistics: F5F1 F5F7 F5F9 F5F17

### Kfold-Crossvalidation with multiple iterations
```{r}
# Load the necessary libraries
library(stuart)

# Number of iterations for averaging
num_iterations <- 3

# Initialize a list to store the results of each iteration
results_list <- vector("list", num_iterations)

# Loop over the number of iterations
for (i in 1:num_iterations) {
  # Perform k-fold cross-validation and store the result
  results_list[[i]] <- kfold(
    'gene', 
    k = 3, 
    data = imputed_training_data_c,
    factor.structure = fs,
    max.invariance = "strict", 
    capacity = list(4, 4, 4, 4, 4), 
    seed = 2028 + i, # Change seed each time for variability
    seeded.search = TRUE,
    auxiliary = NULL,
    software = "lavaan",
    cores = 4,
    objective = obj_default3,
    ignore.errors = TRUE,
    burnin = 5,
    generations = 500,
    individuals = 34,
    selection = "tournament",
    selection.pressure = NULL,
    elitism = NULL,
    reproduction = 0.5,
    mutation = 0.05,
    mating.index = 0,
    mating.size = 0.25,
    mating.criterion = "similarity",
    immigration = 0,
    convergence.criterion = "geno.between",
    tolerance = NULL,
    reinit.n = 1,
    reinit.criterion = "geno.between",
    reinit.tolerance = NULL,
    reinit.prop = 0.75,
    schedule = "run",
    analysis.options = NULL,
    suppress.model = FALSE,
    remove.details = TRUE
  )
}

# Initialize a list to store averaged fit indices
fit_indices <- list(
  SRMR = numeric(num_iterations),
  RMSEA = numeric(num_iterations),
  CFI = numeric(num_iterations)
)

# Extract fit indices from each result
for (i in 1:num_iterations) {
  # Ensure the results are not NULL
  if (!is.null(results_list[[i]])) {
    # Use the inspect function to get the fit measures
    fit_measures <- lavaan::inspect(results_list[[i]]$final, "fit")
    
    # Store the fit indices if they exist
    if (!is.null(fit_measures)) {
      fit_indices$SRMR[i] <- fit_measures["srmr"]
      fit_indices$RMSEA[i] <- fit_measures["rmsea"]
      fit_indices$CFI[i] <- fit_measures["cfi"]
    }
  }
}

# Compute the average of each fit index, ensuring NA values are handled
average_fit_indices <- sapply(fit_indices, function(x) mean(x, na.rm = TRUE))

# Print the average fit indices
print("Average Fit Indices from Results:")
print(average_fit_indices)

# Inspect the complete output for available information
inspect(results_list[[1]]$final, 'fit')
inspect(results_list[[1]]$final, 'est')

lavaan::standardizedSolution(results_list[[i]]$final)

```
Comprehension: F1F2 F1F5 F1F6 F1F14
Evaluation: F2F6 F2F10 F2F12 F2F13
Integration: F3F2 F3F6 F3F7 F3F9
Communication: F4F1 F4F6 F4F7 F4F8
Statistics: F5F1 F5F6 F5F7 F5F18


### Validation via MG-CFA for MI
```{r}
library(lavaan)

cfa_model_m <- '
# Measurement model
Comprehension =~ F1F2 + F1F5 + F1F6 + F1F14 
Evaluation =~ F2F6 + F2F10 + F2F12 + F2F13 
Integration =~ F3F2 + F3F6 + F3F7 + F3F9 
Communication =~ F4F1 + F4F6 + F4F7 + F4F8 
Statistics =~ F5F1 + F5F6 + F5F7 + F5F18 

# Factor variances
Comprehension ~~ Comprehension
Evaluation ~~ Evaluation
Integration ~~ Integration
Communication ~~ Communication
Statistics ~~ Statistics

# Factor correlations
Comprehension ~~ Evaluation
Comprehension ~~ Integration
Comprehension ~~ Communication
Comprehension ~~ Statistics
Evaluation ~~ Integration
Evaluation ~~ Communication
Evaluation ~~ Statistics
Integration ~~ Communication
Integration ~~ Statistics
Communication ~~ Statistics
'
fit_cfa_model_m <- cfa(cfa_model_m, joined_mgcfa_data, estimator = "MLR", group = "split" )
summary(fit_cfa_model_m, fit.measures = TRUE, standardized = TRUE)

inspect(fit_cfa_model_m, 'fit')
```

#### configural to metric
```{r}
fit_cfa_model_metric <- cfa(cfa_model_m, joined_mgcfa_data, estimator = "MLR", group = "split", group.equal= "loadings")

comp_fit_metric <- compareFit(fit_cfa_model_m, fit_cfa_model_metric)
summary(comp_fit_metric)

```
#### configural to metric
Cfi delta: .011 #check?
rmsea delta: .000 #check
srmr delta: .009 #check

#### metric to scalar
```{r}
fit_cfa_model_scalar <- cfa(cfa_model_m, joined_mgcfa_data, estimator = "MLR", group = "split",
                            group.equal = c('loadings','intercepts'))

comp_fit_scalar <- compareFit(fit_cfa_model_metric,fit_cfa_model_scalar)
summary(comp_fit_scalar)

```
#### metric to scalar
Cfi delta: .009 #check
rmsea delta: .000 #check
srmr delta: .002 #check

#### scalar to strict
```{r}
fit_cfa_model_strict <- cfa(cfa_model_m, joined_mgcfa_data, estimator = "MLR", group = "split",
                            group.equal = c('loadings','intercepts','residuals'))

comp_fit_strict <- compareFit(fit_cfa_model_scalar, fit_cfa_model_strict)
summary(comp_fit_strict)

```
#### scalar to strict
Cfi delta: .019 #fail
rmsea delta: .001 #check
srmr delta: .004 #check

### Objective Function with NU (Metric and MLR)
```{r}
obj_default_nu <- function(chisq, df, pvalue, cfi.robust, rmsea.robust, srmr, crel, nu) {
  # Ensure nu is an atomic vector
  if (!is.atomic(nu)) {
    nu <- unlist(nu)
  }
  
  # Calculate the standard deviation of nu
  nu_sd <- sd(nu)
  
  # Compute the objective function
  1/(1 + exp(-10 * (crel - 0.6))) + 
    1/(1 + exp(-10 * (cfi.robust - 0.6))) +
    0.5 * (1 - (1/(1 + exp(-100 * (rmsea.robust - 0.05))))) + 
    0.5 * (1 - (1/(1 + exp(-100 * (srmr - 0.06))))) + 
    1 / (1 + exp(-10 * (nu_sd - 0.1))) # Adjust the threshold as needed
}
```

### Kfold-Crossvalidation with NU (Metric and MLR) - looped
```{r}
library(stuart)

# Number of iterations for averaging
num_iterations <- 3

# Initialize a list to store the results of each iteration
results_list_nu <- vector("list", num_iterations)

# Loop over the number of iterations
for (i in 1:num_iterations) {
  # Perform k-fold cross-validation and store the result
  results_list_nu[[i]] <- kfold(
    'gene', 
    k = 3, 
    data = imputed_training_data_c, 
    factor.structure = fs,  # Pass the variable directly
    max.invariance = "strict", 
    capacity = list(4, 4, 4, 4, 4), 
    seed = 2024,
    seeded.search = TRUE,
    auxiliary = NULL,
    software = "lavaan",
    cores = 4,
    objective = obj_default_nu,
    ignore.errors = TRUE,
    burnin = 5,
    generations = 500,
    individuals = 34,
    selection = "tournament",
    selection.pressure = NULL,
    elitism = NULL,
    reproduction = 0.5,
    mutation = 0.05,
    mating.index = 0,
    mating.size = 0.25,
    mating.criterion = "similarity",
    immigration = 0,
    convergence.criterion = "geno.between",
    tolerance = NULL,
    reinit.n = 1,
    reinit.criterion = "geno.between",
    reinit.tolerance = NULL,
    reinit.prop = 0.75,
    schedule = "run",
    analysis.options = NULL,
    suppress.model = FALSE,
    remove.details = TRUE
  )
}

# Initialize a list to store averaged fit indices
fit_indices_nu <- list(
  SRMR = numeric(num_iterations),
  RMSEA = numeric(num_iterations),
  CFI = numeric(num_iterations)
)

# Extract fit indices from each result
for (i in 1:num_iterations) {
  # Ensure the results are not NULL
  if (!is.null(results_list_nu[[i]])) {
    # Use the inspect function to get the fit measures
    fit_measures <- lavaan::inspect(results_list_nu[[i]]$final, "fit")
    
    # Store the fit indices if they exist
    if (!is.null(fit_measures)) {
      fit_indices_nu$SRMR[i] <- fit_measures["srmr"]
      fit_indices_nu$RMSEA[i] <- fit_measures["rmsea"]
      fit_indices_nu$CFI[i] <- fit_measures["cfi"]
    }
  }
}

# Compute the average of each fit index, ensuring NA values are handled
average_fit_indices <- sapply(fit_indices_nu, function(x) mean(x, na.rm = TRUE))

# Print the average fit indices
print("Average Fit Indices from Results:")
print(average_fit_indices)

# Inspect the complete output for available information
inspect(results_list_nu[[1]]$final, 'fit')
inspect(results_list_nu[[1]]$final, 'est')

lavaan::standardizedSolution(results_list_nu[[i]]$final)
#########################################################
```

### MG-CFA with NU and Metric data 
```{r}
cfa_model_m1 <- '
# Measurement model
Comprehension =~ F1F5 + F1F6 + F1F8 + F1F9 
Evaluation =~ F2F3 + F2F7 + F2F8 + F2F20 
Integration =~ F3F1 + F3F2 + F3F7 + F3F9 
Communication =~ F4F3 + F4F4 + F4F6 + F4F7 
Statistics =~ F5F3 + F5F5 + F5F7 + F5F8 

# Factor variances
Comprehension ~~ Comprehension
Evaluation ~~ Evaluation
Integration ~~ Integration
Communication ~~ Communication
Statistics ~~ Statistics

# Factor correlations
Comprehension ~~ Evaluation
Comprehension ~~ Integration
Comprehension ~~ Communication
Comprehension ~~ Statistics
Evaluation ~~ Integration
Evaluation ~~ Communication
Evaluation ~~ Statistics
Integration ~~ Communication
Integration ~~ Statistics
Communication ~~ Statistics
'
fit_cfa_model_m1 <- cfa(cfa_model_m1, joined_mgcfa_data, estimator = "MLR", group = "split" )
summary(fit_cfa_model_m1, fit.measures = TRUE, standardized = TRUE)

inspect(fit_cfa_model_m1, 'fit')
```

#### configural to metric
```{r}
fit_cfa_model_metric1 <- cfa(cfa_model_m1, joined_mgcfa_data, estimator = "MLR", group = "split", group.equal= "loadings")

comp_fit_metric1 <- compareFit(fit_cfa_model_m1, fit_cfa_model_metric1)
summary(comp_fit_metric1)

```
#### metric to scalar
```{r}
fit_cfa_model_scalar1 <- cfa(cfa_model_m1, joined_mgcfa_data, estimator = "MLR", group = "split",
                             group.equal = c('loadings','intercepts'))

comp_fit_scalar1 <- compareFit(fit_cfa_model_metric1,fit_cfa_model_scalar1)
summary(comp_fit_scalar1)

```

#### scalar to strict
```{r}
fit_cfa_model_strict1 <- cfa(cfa_model_m1, joined_mgcfa_data, estimator = "MLR", group = "split",
                             group.equal = c('loadings','intercepts','residuals'))

comp_fit_strict1 <- compareFit(fit_cfa_model_scalar1, fit_cfa_model_strict1)
summary(comp_fit_strict1)

```

################################################

## K-fold Crossvalidation with Ordinal data Structure (WLSMV)
```{r, eval=FALSE}
library(stuart)
set.seed(125)
# Number of iterations for averaging
num_iterations <- 3

# Initialize a list to store the results of each iteration
results_list_ord <- vector("list", num_iterations)

# Loop over the number of iterations
for (i in 1:num_iterations) {
  # Perform k-fold cross-validation and store the result
  results_list_ord[[i]] <- kfold(
    'gene', 
    k = 3, 
    data = ords1,
    factor.structure = fs,  
    max.invariance = "strong", 
    capacity = list(4, 4, 4, 4, 4), 
    seed = 2026,
    seeded.search = TRUE,
    auxiliary = NULL,
    software = "lavaan",
    cores = 4,
    objective = objective.normal,
    ignore.errors = TRUE,
    burnin = 5,
    generations = 500,
    individuals = 34,
    selection = "tournament",
    selection.pressure = NULL,
    elitism = NULL,
    reproduction = 0.5,
    mutation = 0.05,
    mating.index = 0,
    mating.size = 0.25,
    mating.criterion = "similarity",
    immigration = 0,
    convergence.criterion = "geno.between",
    tolerance = NULL,
    reinit.n = 1,
    reinit.criterion = "geno.between",
    reinit.tolerance = NULL,
    reinit.prop = 0.75,
    schedule = "run",
    analysis.options = NULL,
    suppress.model = FALSE,
    remove.details = TRUE
  )
}

# Initialize a list to store averaged fit indices
fit_indices_ord <- list(
  SRMR = numeric(num_iterations),
  RMSEA = numeric(num_iterations),
  CFI = numeric(num_iterations)
)

# Extract fit indices from each result
for (i in 1:num_iterations) {
  # Ensure the results are not NULL
  if (!is.null(results_list_ord[[i]])) {
    # Use the inspect function to get the fit measures
    fit_measures <- lavaan::inspect(results_list_ord[[i]]$final, "fit")
    
    # Store the fit indices if they exist
    if (!is.null(fit_measures)) {
      fit_indices_ord$SRMR[i] <- fit_measures["srmr"]
      fit_indices_ord$RMSEA[i] <- fit_measures["rmsea"]
      fit_indices_ord$CFI[i] <- fit_measures["cfi"]
    }
  }
}

# Compute the average of each fit index, ensuring NA values are handled
average_fit_indices_ord <- sapply(fit_indices_ord, function(x) mean(x, na.rm = TRUE))

# Print the average fit indices
print("Average Fit Indices from Results:")
print(average_fit_indices_ord)

# Inspect the complete output for available information
inspect(results_list_ord[[1]]$final, 'fit')
inspect(results_list_ord[[1]]$final, 'est')

lavaan::standardizedSolution(results_list_ord[[i]]$final)

```
```{r, eval=FALSE}
library(stuart)
set.seed(127)
# Number of iterations for averaging
num_iterations <- 5

# Initialize a list to store the results of each iteration
results_list_ord1 <- vector("list", num_iterations)

# Loop over the number of iterations
for (i in 1:num_iterations) {
  # Perform k-fold cross-validation and store the result
  results_list_ord1[[i]] <- kfold(
    'gene', 
    k = 3, 
    data = ords1,
    factor.structure = fs,  
    max.invariance = "strong", 
    capacity = list(4, 4, 4, 4, 4), 
    seed = 2026,
    seeded.search = TRUE,
    auxiliary = NULL,
    software = "lavaan",
    cores = 6,
    objective = objective.normal,
    ignore.errors = TRUE,
    burnin = 10,
    generations = 1000,
    individuals = 50,
    selection = "tournament",
    selection.pressure = NULL,
    elitism = NULL,
    reproduction = 0.5,
    mutation = 0.05,
    mating.index = 0,
    mating.size = 0.25,
    mating.criterion = "similarity",
    immigration = 0,
    convergence.criterion = "geno.between",
    tolerance = 1e-6, # stricter tolerance
    reinit.n = 2, # adjusted reinit
    reinit.criterion = "geno.between",
    reinit.tolerance = 1e-6, # stricter reinit tolerance
    reinit.prop = 0.75,
    schedule = "run",
    analysis.options = NULL,
    suppress.model = FALSE,
    remove.details = TRUE
  )
}

# Initialize a list to store averaged fit indices
fit_indices_ord1 <- list(
  SRMR = numeric(num_iterations),
  RMSEA = numeric(num_iterations),
  CFI = numeric(num_iterations)
)

# Extract fit indices from each result
for (i in 1:num_iterations) {
  # Ensure the results are not NULL
  if (!is.null(results_list_ord1[[i]])) {
    # Use the inspect function to get the fit measures
    fit_measures_ord1 <- lavaan::inspect(results_list_ord1[[i]]$final, "fit")
    
    # Store the fit indices if they exist
    if (!is.null(fit_measures)) {
      fit_indices_ord1$SRMR[i] <- fit_measures["srmr"]
      fit_indices_ord1$RMSEA[i] <- fit_measures["rmsea"]
      fit_indices_ord1$CFI[i] <- fit_measures["cfi"]
    }
  }
}

# Compute the average of each fit index, ensuring NA values are handled
average_fit_indices_ord1 <- sapply(fit_indices_ord1, function(x) mean(x, na.rm = TRUE))

# Print the average fit indices
print("Average Fit Indices from Results:")
print(average_fit_indices_ord1)

# Inspect the complete output for available information
inspect(results_list_ord1[[1]]$final, 'fit')
inspect(results_list_ord1[[1]]$final, 'est')

lavaan::standardizedSolution(results_list_ord1[[i]]$final)

```

```{r}
# Measurement model before fuck up
#Comprehension =~ F1F4 + F1F6 + F1F10 + F1F11 
#Evaluation =~ F2F1 + F2F6 + F2F13 + F2F16 
#Integration =~ F3F5 + F3F6 + F3F7 + F3F8 
#Communication =~ F4F3 + F4F6 + F4F7 + F4F8 
#Statistics =~ F5F6 + F5F8 + F5F9 + F5F16 

cfa_model_ord <- '
# Measurement model
Comprehension =~ F1F5 + F1F6 + F1F8 + F1F9 
Evaluation =~ F2F10 + F2F14 + F2F15 + F2F20 
Integration =~ F3F2 + F3F7 + F3F8 + F3F8 
Communication =~ F4F3 + F4F4 + F4F6 + F4F7 
Statistics =~ F5F1 + F5F5 + F5F8 + F5F15 

# Factor variances
Comprehension ~~ Comprehension
Evaluation ~~ Evaluation
Integration ~~ Integration
Communication ~~ Communication
Statistics ~~ Statistics

# Factor correlations
Comprehension ~~ Evaluation
Comprehension ~~ Integration
Comprehension ~~ Communication
Comprehension ~~ Statistics
Evaluation ~~ Integration
Evaluation ~~ Communication
Evaluation ~~ Statistics
Integration ~~ Communication
Integration ~~ Statistics
Communication ~~ Statistics
'
fit_cfa_model_ords <- cfa(cfa_model_ord, ords, ordered = TRUE, estimator = "WLSMV")
summary(fit_cfa_model_ords, fit.measures = TRUE, standardized = TRUE)

inspect(fit_cfa_model_ords,'fit')
```

### MG-CFA with Ordinal data Structure (WLSMV)
```{r}
fit_cfa_model_ords <- cfa(cfa_model_ord, ords_mi, ordered = TRUE, estimator = "WLSMV", group = "split")
summary(fit_cfa_model_ords, fit.measures = TRUE, standardized = TRUE)

inspect(fit_cfa_model_ords,'fit')
```

#### configural to metric
```{r}
library(lavaan)
library(semTools)
fit_cfa_model_ords_m <- cfa(cfa_model_ord, ords_mi, estimator = "WLSMV", group = "split", group.equal= "loadings")

comp_fit_ord <- compareFit(fit_cfa_model_ords, fit_cfa_model_ords_m)
summary(comp_fit_ord)
inspect(fit_cfa_model_ords_m ,'fit')
```