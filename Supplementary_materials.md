---
title: "Coming out from the shadows: facultative slave-making ants reveal their chemical identity during colony development"
subtitle: "Behavioral Ecology and Sociobiology"
author:
  - Tomasz Włodarczyk^[Faculty of Biology, University of Białystok, Poland, tomwlo@gmail.com, ORCID 0000-0003-1554-9699]
  - Thomas Schmitt^[Department of Animal Ecology and Tropical Biology, University of Würzburg, Germany, thomas.schmitt@uni-wuerzburg.de, ORCID 0000-0002-6719-8635]
bibliography: book.bib
params:
  save_globals: TRUE
clean: false
output:
  bookdown::pdf_document2:
    keep_tex: true
---

# Introduction
This document presents graphical results generated as output from analytic workflows written in *R* language and available in [the public GitHub repository - https://github.com/TomVuod/sangchem/tree/supp](https://github.com/TomVuod/sangchem/tree/supp).   





<!--chapter:end:index.Rmd-->

# Changes in CHC amount over time

This section contains a series of linear (mixed) models fitted to examine how the proportion of *F. sanguinea* workers in a colony is associated with the amount of CHC, analyzed as a whole or in subsets. Included are summary reports, results of the Shapiro-Wilk test of normality of conditional residuals, as well as diagnostic plots and tests generated with the use of `DHARMa` R package[@R-DHARMa]. The response variable often needed to be transformed to meet model assumptions. Random terms with no variance have been dropped. 

In model specification $(1|\text{random_factor})$ denotes random intercept (separate for each level of the random factor) and $((1|\text{random_factor_1:random_factor_2}))$ denotes random intercepts generated from the combinations of two factors. $\beta_{0}$ and $\epsilon$ denote intercept and error term, respectively. 



## Colony growth and body size

n the analyses below, normalized CHC amounts were used as the response variable. Because body shape in both species deviates from an isometric growth pattern (see Section \@ref(surface-area)), normalizing CHC amounts by the square of head width may introduce bias. To address this, we first examined the relationship between body size (head width) and the proportion of F. sanguinea workers in a colony, and then incorporated the square of head width as an explanatory variable in the linear model.

$$
\text{head_width} =  \beta_0 + \text{sanguinea_proportion} + (1|\text{colony_ID}) + \epsilon,
$$


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: head_width ~ sang_prop + (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: -103.8
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.10138 -0.61076  0.06058  0.43106  2.28980 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.003155 0.05617 
##  Residual             0.009110 0.09545 
## Number of obs: 68, groups:  colony, 16
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  1.21806    0.02641 38.92702  46.121  < 2e-16 ***
## sang_prop    0.12514    0.04344 62.09975   2.881  0.00544 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.715
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.9861, p-value = 0.6514
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-3-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-3-2} \end{center}

## Change in the total CHC mass

In this series of models the total CHC mass (normalized by assuming head width of 1.3 mm) was regressed on the proportion of *Formica sanguinea* workers in the colony. 

### Mature *F. sanguinea*

$$
\log(\text{CHC_mass_sanguinea}) =  \beta_0 + \text{sanguinea_proportion} +\newline (1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_sanguinea}$ denotes the total normalized CHC mass on the body of adult *F. sanguinea* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass)) ~ sang_prop + (1 | colony) + (1 | colony:census_date)
##    Data: model_input
## 
## REML criterion at convergence: 92.1
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.83723 -0.39249  0.01222  0.53571  1.73459 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.12975  0.3602  
##  colony             (Intercept) 0.02626  0.1620  
##  Residual                       0.10333  0.3215  
## Number of obs: 68, groups:  colony:census_date, 44; colony, 16
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)   1.2035     0.1295 36.0905   9.296 4.11e-11 ***
## sang_prop     0.7991     0.2392 36.9985   3.340  0.00192 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.786
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98623, p-value = 0.6591
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-4-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-4-2} \end{center}
The model shows no significant effect of body size on CHC amount in *F. sanguinea*, and this term will therefore be excluded from subsequent models for this species.

### *F. fusca*

$$
\log(\text{CHC_mass_fusca}) =  \beta_0 + \text{sanguinea_proportion} + (1|\text{colony_ID}) + \epsilon,
$$
where $\text{CHC_mass_fusca}$ denotes the total normalized CHC mass on the body of adult *F. fusca* worker. 


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass)) ~ sang_prop + (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 103.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.79685 -0.62513 -0.00685  0.63754  2.49220 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.05014  0.2239  
##  Residual             0.17722  0.4210  
## Number of obs: 78, groups:  colony, 20
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.81635    0.08805 31.60123   9.272 1.57e-10 ***
## sang_prop    0.17928    0.17517 70.00046   1.024     0.31    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.594
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.99196, p-value = 0.9107
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-5-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-5-2} \end{center}

### Callow *F. sanguinea*

$$
\log(\text{CHC_mass_callow}) =  \beta_0 + \text{sanguinea_proportion} + \newline(1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_callow}$ denotes the total normalized CHC mass on the body of callow *F. sanguinea* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass)) ~ sang_prop + (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 51.2
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.16539 -0.61246  0.07593  0.62418  1.49264 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.06331  0.2516  
##  Residual             0.22509  0.4744  
## Number of obs: 32, groups:  colony, 16
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)   
## (Intercept)  0.03373    0.17027 28.88675   0.198  0.84435   
## sang_prop    1.05723    0.30345 25.91578   3.484  0.00177 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.775
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.95692, p-value = 0.2258
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-6-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-6-2} \end{center}
Similarly to mature *F. sanguinea* ants, in callow workers there is negative but not significant relationship between body size and normalized CHC amount.

## Change of the CHC characteristic of *F. sanguinea*

In this section the subset of peaks were used to calculate normalized CHC mass. These were peaks identified as *F. sanguinea* markers.

### Mature *F. sanguinea* ants

$$
\log(\text{CHC_mass_mature}) =  \beta_0 + \text{sanguinea_proportion} +  \newline(1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_mature}$ denotes the normalized mass of *F. sanguinea* markers on the body of adult *F. sanguinea* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass + 1)) ~ sang_prop + (1 | colony:census_date) +      (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 30.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.8821 -0.3658 -0.0362  0.3439  2.1883 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.06558  0.2561  
##  colony             (Intercept) 0.00000  0.0000  
##  Residual                       0.03804  0.1950  
## Number of obs: 68, groups:  colony:census_date, 44; colony, 16
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)   0.2833     0.0805 40.1068   3.519  0.00109 ** 
## sang_prop     1.3131     0.1564 39.5702   8.394 2.54e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.822
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98702, p-value = 0.7046
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-7-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-7-2} \end{center}

### *F. fusca* ants

$$
\sqrt{(\text{CHC_mass_fusca})} =  \beta_0 + \text{sanguinea_proportion} +  \newline(1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_fusca}$ denotes the normalized mass of *F. sanguinea* markers on the body of adult *F. sanguinea* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(sqrt(normalized_mass)) ~ sang_prop + (1 | colony:census_date)
##    Data: model_input
## 
## REML criterion at convergence: -58.3
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.06253 -0.51814 -0.07238  0.45433  1.95590 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.009805 0.09902 
##  Residual                       0.016560 0.12869 
## Number of obs: 78, groups:  colony:census_date, 55
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.41024    0.02938 54.47670  13.965  < 2e-16 ***
## sang_prop    0.61852    0.06935 51.99602   8.919 4.63e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.729
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.99067, p-value = 0.8465
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-8-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-8-2} \end{center}

### Callow *F. sanguinea* ants

$$
\sqrt[3]{(\text{CHC_mass_mature})} = \beta_0 + \text{sanguinea_proportion} + \epsilon,
$$
where $\text{CHC_mass_mature}$ denotes the normalized mass of *F. sanguinea* markers on the body of callow *F. sanguinea* worker.


```
## 
## Call:
## lm(formula = I((normalized_mass)^(1/3)) ~ sang_prop, data = model_input)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.41432 -0.06803  0.02430  0.07590  0.35296 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.31592    0.04638   6.812 1.48e-07 ***
## sang_prop    0.71797    0.08719   8.235 3.42e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1473 on 30 degrees of freedom
## Multiple R-squared:  0.6933,	Adjusted R-squared:  0.6831 
## F-statistic: 67.81 on 1 and 30 DF,  p-value: 3.423e-09
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.96738, p-value = 0.4305
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-9-1} \end{center}

## Change of the CHC characteristic of *F. fusca*

Similarly to the analysis with the use of markers of *F. sanguinea* workers, the procedure was repeated using the mass of *F. fusca* markers as a response variable. 

### Mature *F. sanguinea* ants

$$
\log(\text{CHC_mass_mature}) =  \beta_0 + \text{sanguinea_proportion} + (1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_mature}$ denotes the normalized mass of *F. fusca* markers on the body of adult *F. sanguinea* worker.



```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass)) ~ sang_prop + (1 | colony:census_date) +      (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 131.2
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.82877 -0.45873  0.05316  0.55616  1.98511 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.1471   0.3835  
##  colony             (Intercept) 0.2359   0.4857  
##  Residual                       0.1804   0.4247  
## Number of obs: 68, groups:  colony:census_date, 44; colony, 16
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)  -0.4566     0.1947 33.8202  -2.345 0.025027 *  
## sang_prop    -1.1016     0.2951 33.2726  -3.733 0.000707 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.656
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98971, p-value = 0.8512
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-10-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-10-2} \end{center}

### *F. fusca* ants

$$
\log(\text{CHC_mass_mature}+0.01) =  \beta_0 + \text{sanguinea_proportion} + \newline(1|\text{colony_ID}) + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{CHC_mass_mature}$ denotes the normalized mass of *F. fusca* markers on the body of *F. fusca* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass + 0.01)) ~ sang_prop + (1 | colony) + (1 |      colony:census_date)
##    Data: model_input
## 
## REML criterion at convergence: 159.8
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.34120 -0.60468 -0.05522  0.59876  2.14890 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.0000   0.0000  
##  colony             (Intercept) 0.2325   0.4821  
##  Residual                       0.3213   0.5668  
## Number of obs: 78, groups:  colony:census_date, 55; colony, 20
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)  -0.8547     0.1462 27.2550  -5.844 3.08e-06 ***
## sang_prop    -0.8184     0.2436 65.7807  -3.360   0.0013 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.479
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.99138, p-value = 0.8835
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-11-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-11-2} \end{center}

### Callow *F. sanguinea* ants

$$
\log(\text{CHC_mass_mature}) = \beta_0 + \text{sanguinea_proportion} + (1|\text{colony_ID}) + \epsilon,
$$
where $\log(\text{CHC_mass_mature})$ denotes the normalized mass of *F. fusca* markers on the body of callow *F. sanguinea* worker.


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(log(normalized_mass)) ~ sang_prop + (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 59.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.88270 -0.59925 -0.01403  0.65154  1.88156 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.2208   0.4699  
##  Residual             0.2225   0.4717  
## Number of obs: 32, groups:  colony, 16
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept) -1.66659    0.20293 28.05257  -8.213 6.04e-09 ***
## sang_prop   -0.07903    0.32299 19.69963  -0.245    0.809    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.685
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98634, p-value = 0.9487
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-12-1} \end{center}

<!--chapter:end:CHCs_over_time.Rmd-->

# Comparison of the amount and proportions of CHC between species and age categories

This section presents the results of Wilcoxon tests comparing features of CHC profiles between callow and mature *F. sanguinea* ants as well as *F. fusca* slaves. Samples from the same colony were averaged before calculation of the final statistics to account for their non-independence.  



## Difference in the total CHC amount between *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-14-1} \end{center}

## Difference in the total CHC amount between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-15-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. sanguinea* between mature *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-16-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-17-1} \end{center}

## Difference in the proportion of CHC characteristic of callow *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-18-1} \end{center}


## Difference in the amount of CHC characteristic of callow *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-19-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. fusca* between mature *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-20-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. fusca* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-21-1} \end{center}


## Difference in the proportion of CHC characteristic of *F. fusca* between callow *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-22-1} \end{center}

## Difference in the mass of *n*-alkanes between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-23-1} \end{center}

## Difference in the proprtion of *n*-alkanes between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-24-1} \end{center}

<!--chapter:end:Non-parametric_tests.Rmd-->

# Change in the CHC profile of separated callow *F. sanguinea* ants


## Change in the total CHC amount over time

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(sqrt(normalized_mass)) ~ mean_delta + (1 | colony)
##    Data: separation_data[-74, ]
## 
## REML criterion at convergence: 9.8
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.55998 -0.69772 -0.07088  0.67280  2.51792 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.01765  0.1329  
##  Residual             0.04716  0.2172  
## Number of obs: 79, groups:  colony, 11
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)  1.367210   0.055393 16.298104  24.682 2.43e-14 ***
## mean_delta   0.002934   0.001933 68.480649   1.518    0.134    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## mean_delta -0.515
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.9908, p-value = 0.8482
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-26-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-26-2} \end{center}

## Change in the amount of compounds characteristic of callow *F. sanguinea*

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: log(normalized_mass_part) ~ mean_delta + (1 | colony)
##    Data: separation_data
## 
## REML criterion at convergence: 92.1
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.3303 -0.7669  0.1012  0.6201  2.5706 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.1129   0.3360  
##  Residual             0.1238   0.3518  
## Number of obs: 80, groups:  colony, 11
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)
## (Intercept) -0.003285   0.118687 13.016273  -0.028    0.978
## mean_delta  -0.004592   0.003054 68.334089  -1.503    0.137
## 
## Correlation of Fixed Effects:
##            (Intr)
## mean_delta -0.389
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98886, p-value = 0.7216
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-27-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-27-2} \end{center}

## Change in the amount of compounds characteristic of *F. sanguinea*

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(sqrt(normalized_mass_part)) ~ poly(mean_delta, 1) + (1 | colony)
##    Data: separation_data[-c(69, 74), ]
## 
## REML criterion at convergence: -67.6
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.27083 -0.64402  0.06067  0.73259  2.42696 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.002444 0.04944 
##  Residual             0.021031 0.14502 
## Number of obs: 78, groups:  colony, 11
## 
## Fixed effects:
##                     Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          0.60656    0.02252  7.88363  26.930 4.81e-09 ***
## poly(mean_delta, 1)  1.08605    0.14693 69.21997   7.392 2.54e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## ply(mn_d,1) 0.020
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.98936, p-value = 0.7684
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-28-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-28-2} \end{center}

## Change in the amount of compounds characteristic of *F. fusca*

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: normalized_mass_part ~ mean_delta + (1 | colony)
##    Data: separation_data
## 
## REML criterion at convergence: -25.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.5589 -0.4676 -0.0943  0.2358  6.8831 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.02205  0.1485  
##  Residual             0.02765  0.1663  
## Number of obs: 80, groups:  colony, 11
## 
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)  2.240e-01  5.346e-02  1.383e+01    4.19 0.000932 ***
## mean_delta  -4.298e-05  1.443e-03  6.870e+01   -0.03 0.976326    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## mean_delta -0.408
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.59858, p-value = 1.852e-13
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-29-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-29-2} \end{center}

<!--chapter:end:Separation_experiment.Rmd-->

# Impact of slaves on the development of the *F. sanguinea* CHC profile

We analyzed the impact of *F. fusca* slaves on the CHC profile of callow *F. sanguinea* ants. In this experiment, *F. sanguinea* ants were isolated in pairs before eclosion, and their CHC profiles were analyzed at one of four time points marking their adult age. We then tested whether callow *F. sanguinea* ants were chemically more similar to their slave relatives from free-living colonies. Unrelated *F. fusca* ants served as a background group."

## Analysis based on all samples


\begin{figure}

{\centering \includegraphics[width=0.8\linewidth]{Supplementary_materials_files/figure-latex/unnamed-chunk-31-1} 

}

\caption{Distribution of lifespans of ants used in the experiment. The precise age of individuals could not be determined because ants in pairs were not marked. Therefore, all possible combinations were considered to illustrate the potential range.}(\#fig:unnamed-chunk-31)
\end{figure}





```
## Empirical p-value for the chemical distance difference for the ants aged 1-3 days:0.00000
```

```
## Number of colonies used in the treatment of 1-3 days of separation: 9
```



```
## Empirical p-value for the chemical distance difference for the ants aged 8-10 days:0.00000
```

```
## Number of colonies used in the treatment of 8-10 days of separation: 9
```


```
## Empirical p-value for the chemical distance difference for the ants aged 17-20 days:0.00079
```

```
## Number of colonies used in the treatment of 17-20 days of separation: 9
```


```
## Empirical p-value for the chemical distance difference for the ants aged 35-40 days:0.00703
```

```
## Number of colonies used in the treatment of 35-40 days of separation: 8
```

\begin{table}

\caption{(\#tab:unnamed-chunk-40)Chemical distance of the isolated \textit{F. sanguinea} ant to \textit{F. fusca} ants realted and unrelated to slaves.}
\centering
\begin{tabular}[t]{l|>{\raggedleft\arraybackslash}p{4cm}|>{\raggedleft\arraybackslash}p{4cm}}
\hline
Adult age & Chemical distance to slaves' relatives & Chemical distance to random colonies\\
\hline
1-3 days & 0.3332 & 0.4718\\
\hline
8-10 days & 0.3518 & 0.4802\\
\hline
17-20 days & 0.4521 & 0.5291\\
\hline
35-40 days & 0.5729 & 0.6100\\
\hline
\end{tabular}
\end{table}


## Analysis restricted to the individuals hatched from cocoons

Since naked pupae might have acquired CHC through physical contact with adult ants after pupation, we subset our data to include only those samples which derived from the ants hatched from cocoons. In this case silky envelop during pupa stage ruled out possibility of the CHC transfer from the environment. 





```
## Empirical p-value for the chemical distance difference for the ants aged 1-3 days:0.00000
```

```
## Number of colonies used in the treatment of 1-3 days of separation: 8
```


```
## Empirical p-value for the chemical distance difference for the ants aged 8-10 days:0.00001
```

```
## Number of colonies used in the treatment of 8-10 days of separation: 8
```



```
## Empirical p-value for the chemical distance difference for the ants aged 17-20 days:0.00438
```

```
## Number of colonies used in the treatment of 17-20 days of separation: 8
```




```
## Empirical p-value for the chemical distance difference for the ants aged 35-40 days:0.03298
```

```
## Number of colonies used in the treatment of 35-40 days of separation: 7
```

\begin{table}

\caption{(\#tab:unnamed-chunk-51)Chemical distance of the isolated \textit{F. sanguinea} ant to \textit{F. fusca} ants realted and unrelated to slaves. Isolated ant spun a silky envelope before pupal stage, which prevented CHC transfer from the enviroment.}
\centering
\begin{tabular}[t]{l|>{\raggedleft\arraybackslash}p{4cm}|>{\raggedleft\arraybackslash}p{4cm}}
\hline
Adult age & Chemical distance to slaves' relatives & Chemical distance to random colonies\\
\hline
1-3 days & 0.3473 & 0.4869\\
\hline
8-10 days & 0.3740 & 0.4932\\
\hline
17-20 days & 0.4643 & 0.5407\\
\hline
35-40 days & 0.6175 & 0.6458\\
\hline
\end{tabular}
\end{table}


<!--chapter:end:Slave_impact.Rmd-->

# Effect of dummy ants

In these experiments, *F. sanguinea* ants released from their cocoon envelopes were maintained in Petri dishes along with glass beads coated with CHC extracted from either *F. fusca* or *F. sanguinea* ants. We examined whether the contact with dummy nestmates influenced the development of CHC profile of young *F. sanguinea* workers. In particular, we aimed to determine there there is evidence for active chemical mimicry in *F. sanguinea*. 



## Distance to the CHC profile of the dummy ants treated by *F. fusca* CHC

We computed the chemical distance to the CHC profile of ants whose CHC were used to cover glass beads serving as dummy ants. In the control variant, we used the ants subjected to clean glass beads. 


```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  to_test$mean_dist.x and to_test$mean_dist.y
## V = 18.5, p-value = 0.6781
## alternative hypothesis: true location shift is not equal to 0
```


\begin{table}

\caption{(\#tab:unnamed-chunk-55)Chemical distance of separated \textit{F. sanguinea} ants to the CHC profile of \textit{F. fusca} ants from colonies that served as a source of CHC to caot the glass beads. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{l|l|>{\raggedleft\arraybackslash}p{2cm}|l|>{\raggedleft\arraybackslash}p{2cm}}
\hline
Colony ID & Treamtent & Averaged distance & Contrast & Averaged distance\\
\hline
SD18-3 & 12-15 days + F. fusca hydrocarbons & 0.396 & 12-15 days (control) & 0.509\\
\hline
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.616 & 12-15 days (control) & 0.513\\
\hline
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.842 & 12-15 days (control) & 0.832\\
\hline
SD19-11 & 12-15 days + F. fusca hydrocarbons & 0.467 & 12-15 days (control) & 0.469\\
\hline
SD19-3 & 12-15 days + F. fusca hydrocarbons & 0.500 & 12-15 days (control) & 0.436\\
\hline
SD19-4 & 12-15 days + F. fusca hydrocarbons & 0.305 & 12-15 days (control) & 0.424\\
\hline
SD19-6 & 12-15 days + F. fusca hydrocarbons & 0.674 & 12-15 days (control) & 0.699\\
\hline
SD19-8 & 12-15 days + F. fusca hydrocarbons & 0.595 & 12-15 days (control) & 0.587\\
\hline
SD20-2 & 12-15 days + F. fusca hydrocarbons & 0.699 & 12-15 days (control) & 0.709\\
\hline
\end{tabular}
\end{table}
## Distance to the CHC profile of the dummy ants treated by *F. sanguinea* CHC

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  to_test$mean_dist.x and to_test$mean_dist.y
## V = 12, p-value = 0.2361
## alternative hypothesis: true location shift is not equal to 0
```

\begin{table}

\caption{(\#tab:unnamed-chunk-56)Chemical distance of separated \textit{F. sanguinea} ants to the CHC profile of \textit{F. sanguinea} ants from colonies that served as a source of CHC to caot the glass beads. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{l|l|>{\raggedleft\arraybackslash}p{2cm}|l|>{\raggedleft\arraybackslash}p{2cm}}
\hline
Colony ID & Treamtent & Averaged distance & Contrast & Averaged distance\\
\hline
SD18-2 & 12-15 days + F. sanguinea hydrocarbons & 0.571 & 12-15 days (control) & 0.373\\
\hline
SD18-3 & 12-15 days + F. sanguinea hydrocarbons & 0.457 & 12-15 days (control) & 0.599\\
\hline
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.525 & 12-15 days (control) & 0.580\\
\hline
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.250 & 12-15 days (control) & 0.232\\
\hline
SD19-11 & 12-15 days + F. sanguinea hydrocarbons & 0.595 & 12-15 days (control) & 0.607\\
\hline
SD19-3 & 12-15 days + F. sanguinea hydrocarbons & 0.419 & 12-15 days (control) & 0.511\\
\hline
SD19-4 & 12-15 days + F. sanguinea hydrocarbons & 0.648 & 12-15 days (control) & 0.668\\
\hline
SD19-6 & 12-15 days + F. sanguinea hydrocarbons & 0.553 & 12-15 days (control) & 0.591\\
\hline
SD19-8 & 12-15 days + F. sanguinea hydrocarbons & 0.850 & 12-15 days (control) & 0.850\\
\hline
SD20-2 & 12-15 days + F. sanguinea hydrocarbons & 0.646 & 12-15 days (control) & 0.654\\
\hline
\end{tabular}
\end{table}



## Fraction of *n*-docosane in the CHC profile

For the control treatment, glass beads were left uncoated with CHC. In the experimental treatments, the CHC used to coat the glass beads were supplemented with *n*-docosane, which served as an internal standard. We examined differences in the proportion of *n*-docosane in the CHC profiles of the tested ants, as these could indicate the acquisition of chemicals from the dummy ants.

### *F. fusca*


```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  to_test$C22_prop.x and to_test$C22_prop.y
## V = 55, p-value = 0.005857
## alternative hypothesis: true location shift is not equal to 0
```

\begin{table}

\caption{(\#tab:unnamed-chunk-57)Proportion of \textit{n}-docosane in CHC extracted from \textit{F. sanguinea} ants maintained with the glass beads coated with the CHC of \textit{F. fusca} ants and contaminated with \textit{n}-docosane. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{l|l|>{\raggedleft\arraybackslash}p{2cm}|l|>{\raggedleft\arraybackslash}p{2cm}}
\hline
Colony ID & Treamtent & Averaged C22 fraction & Contrast & Averaged C22 fraction\\
\hline
SD18-3 & 12-15 days + F. fusca hydrocarbons & 0.018 & 12-15 days (control) & 0.000\\
\hline
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.001 & 12-15 days (control) & 0.000\\
\hline
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
\hline
SD19-11 & 12-15 days + F. fusca hydrocarbons & 0.011 & 12-15 days (control) & 0.002\\
\hline
SD19-3 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.001\\
\hline
SD19-4 & 12-15 days + F. fusca hydrocarbons & 0.017 & 12-15 days (control) & 0.000\\
\hline
SD19-6 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.001\\
\hline
SD19-8 & 12-15 days + F. fusca hydrocarbons & 0.003 & 12-15 days (control) & 0.002\\
\hline
SD20-2 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.005\\
\hline
W17-1 & 12-15 days + F. fusca hydrocarbons & 0.021 & 12-15 days (control) & 0.001\\
\hline
\end{tabular}
\end{table}

### *F. sanguinea*


```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  to_test$C22_prop.x and to_test$C22_prop.y
## V = 54, p-value = 0.007093
## alternative hypothesis: true location shift is not equal to 0
```


\begin{longtable}[t]{l|l|>{\raggedleft\arraybackslash}p{2cm}|l|>{\raggedleft\arraybackslash}p{2cm}}
\caption{(\#tab:unnamed-chunk-58)Proportion of *n*-docosane in CHC extracted from \textit{F. sanguinea} ants maintained with the glass beads coated with the CHC of \textit{F. sanguinea} ants and contaminated with \textit{n}-docosane. In control variant, glass bead were left clean.}\\
\hline
Colony ID & Treamtent & Averaged C22 fraction & Contrast & Averaged C22 fraction\\
\hline
SD18-2 & 12-15 days + F. sanguinea hydrocarbons & 0.005 & 12-15 days (control) & 0.002\\
\hline
SD18-3 & 12-15 days + F. sanguinea hydrocarbons & 0.014 & 12-15 days (control) & 0.000\\
\hline
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.000\\
\hline
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
\hline
SD19-11 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.002\\
\hline
SD19-3 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.001\\
\hline
SD19-4 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
\hline
SD19-6 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.001\\
\hline
SD19-8 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.002\\
\hline
SD20-2 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.005\\
\hline
W17-1 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.001\\
\hline
\end{longtable}

<!--chapter:end:Dummy_ants.Rmd-->

# Body surface area of *F. sanguinea* ants and their slaves {#surface-area}

In the publication, CHC amounts are reported relative to the square of head width. To assess whether this measure reliably reflects body size, we measured the planar projection areas of different body parts of *F. sanguinea* and *F. fusca*. The plots below suggest that head width slightly underestimates body surface area in *F. fusca* compared to *F. sanguinea*. As a result, the relative concentration of CHC in *F. fusca* may be overestimated. However, this bias acts in the opposite direction of our findings—*F. sanguinea* exhibits greater CHC concentration on its cuticle. Therefore, the shape-related bias between species is conservative, as it reduces the likelihood of rejecting the null hypothesis  

\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-59-1} 

}

\caption{Ratios of the planar projection areas of three body parts (head, thorax dorsal view, thorax lateral view) to the square of head width.}(\#fig:unnamed-chunk-59)
\end{figure}

<!--chapter:end:Body_surface.Rmd-->

