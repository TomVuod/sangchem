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

# Species discrimination and prediction
To classify the CHC samples collected from mixed colonies, we need a measure which positions a sample between the *F. sanguinea* and *F. fusca* CHC profiles. Therefore, using `mixOmics` R package [@R-mixOmics], we performed a discriminant analysis training the model on the chemical data collected from the filed colonies in the earlier study [@WLODARCZYK201798]. The model was subsequently used to predict identity of new samples.

## Species markers
To figure out the importance of each input variable in the prediction of the species identity made by discriminant model, we calculated the vector of weights using the following formula [@Tenenhaus, @Rohart].
$$
\textbf{W}(\textbf{P}^\intercal\textbf{W})^{-1}\textbf{c}
$$
$$
\textbf{P}=\textbf{X}\textbf{V}
$$
$$
\textbf{c}=\textbf{V}^\intercal\textbf{y},
$$
where\newline
$\textbf{W}$ is the matrix of loading vectors with the number of rows corresponding to the number of input variables and the number of columns corresponding to the number of the latent components\newline
$\textbf{V}$ is the matrix of the coordinates of each observation on latent components\newline 
$\textbf{X}$ is the input data matrix\newline
$\textbf{y}$ is the vector of zero-centered one-hot encoded class labels of the samples for a selected class (here, *F. sanguinea* sample)

To project the weights assigned to principal components onto the original variables one has to reverse data transformation. Technically, this is achieved by multiplying with the pseudoinverse of the rotation matrix produced by PCA. We use pseudoinverse because the rotation matrix has been reduced by removing principal components with low variance.

$$
\textbf{R}^{+}\textbf{W}(\textbf{P}^\intercal\textbf{W})^{-1}\textbf{c},
$$

where $\textbf{R}^+$ denotes the pseudoinverse of rotation matrix from PCA.







\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-5-1} 

}

\caption{Variance of each principal compnent. Data was fed into the discriminat model to find features separating both species.}(\#fig:unnamed-chunk-5)
\end{figure}





\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-9-1} 

}

\caption{Distribution of the samples projected onto two first latent components before optimizing the number of model paramaters.}(\#fig:unnamed-chunk-9)
\end{figure}



\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-11-1} 

}

\caption{Results of the performance test of the discriminant model with different number of the latent components. For details see the documunetation of R mixOmics package.}(\#fig:unnamed-chunk-11)
\end{figure}



\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-13-1} 

}

\caption{Projection of the samples in onto two first discriminant analysis components after model tuning.}(\#fig:unnamed-chunk-13)
\end{figure}



\begin{longtable}[t]{rlrll}
\caption{(\#tab:unnamed-chunk-15)List of the peaks with the species identity score indicating their importance as a marker of \textit{F. fusca} or \textit{F. sanguinea} samples.}\\
\toprule
Peak ID & Compound & Marker score & fusca marker & sanguinea marker\\
\midrule
\endfirsthead
\caption[]{(\#tab:unnamed-chunk-15)List of the peaks with the species identity score indicating their importance as a marker of \textit{F. fusca} or \textit{F. sanguinea} samples. \textit{(continued)}}\\
\toprule
Peak ID & Compound & Marker score & fusca marker & sanguinea marker\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
15 & 4-MeC24 & -0.0255545 & TRUE & FALSE\\
60 & 13,17-diMeC29 & -0.0253947 & TRUE & FALSE\\
10 & 3-MeC23 & -0.0199077 & TRUE & FALSE\\
6 & 11-Me-C23 + 9-Me-C23 & -0.0187827 & TRUE & FALSE\\
32 & 14-; 10-MeC26 + 3,7,11-triMeC25 & -0.0181734 & TRUE & FALSE\\
\addlinespace
12 & 3,13-; 3,11-; 3,9-; 3,7-diMeC23 & -0.0176211 & TRUE & FALSE\\
14 & 6-MeC24 & -0.0174156 & TRUE & FALSE\\
7 & 7-Me-C23 & -0.0171928 & TRUE & FALSE\\
13 & 12-; 11-; 10-; 9-; 8-MeC24 & -0.0170398 & TRUE & FALSE\\
69 & C31 & -0.0169988 & TRUE & FALSE\\
\addlinespace
44 & 7-MeC27 & -0.0168845 & TRUE & FALSE\\
37 & 7,x-; 5,x-; 10,14-diMeC26 & -0.0168292 & TRUE & FALSE\\
31 & 3,13-diMeC25 & -0.0164987 & TRUE & FALSE\\
19 & 4,12-diMeC24 & -0.0160872 & TRUE & FALSE\\
11 & C24 & -0.0151770 & TRUE & FALSE\\
\addlinespace
22 & 13-; 11-; 9-MeC25 & -0.0150263 & FALSE & FALSE\\
49 & 7,11-diMeC27 & -0.0145556 & FALSE & FALSE\\
5 & C23 & -0.0145410 & FALSE & FALSE\\
33 & 6-MeC26 & -0.0145371 & FALSE & FALSE\\
27 & 7,11-diMeC25 + 3-MeC25 & -0.0141446 & FALSE & FALSE\\
\addlinespace
55 & C29-ene & -0.0120979 & FALSE & FALSE\\
40 & 4,12-diMeC26 & -0.0117432 & FALSE & FALSE\\
21 & 4,12,16-triMeC24 & -0.0115601 & FALSE & FALSE\\
65 & 3,7-; 3,9-; 3,11-diMeC29 & -0.0109379 & FALSE & FALSE\\
3 & C23-ene & -0.0105481 & FALSE & FALSE\\
\addlinespace
8 & 5-Me C23 & -0.0099939 & FALSE & FALSE\\
57 & 4,8,12-triMeC28 & -0.0097465 & FALSE & FALSE\\
42 & 4,8,14-MeC26 & -0.0092601 & FALSE & FALSE\\
26 & 11,15-; 9,13-diMeC25 & -0.0074051 & FALSE & FALSE\\
63 & C30-ene & -0.0073431 & FALSE & FALSE\\
\addlinespace
34 & 5-MeC26 & -0.0072263 & FALSE & FALSE\\
53 & 14-; 12-; 10-MeC28 & -0.0059975 & FALSE & FALSE\\
16 & 10,14-diMeC24 & -0.0049985 & FALSE & FALSE\\
56 & C29 & -0.0026256 & FALSE & FALSE\\
24 & 7-MeC25 & -0.0024319 & FALSE & FALSE\\
\addlinespace
58 & 15-; 13-; 11-; 9-; 7-MeC29 & -0.0017646 & FALSE & FALSE\\
23 & 9-MeC25 & -0.0006146 & FALSE & FALSE\\
20 & C25 & -0.0002985 & FALSE & FALSE\\
77 & C33-ene & 0.0000217 & FALSE & FALSE\\
43 & 13-MeC27 & 0.0000912 & FALSE & FALSE\\
\addlinespace
51 & C28 + 3,15-diMeC27 & 0.0002801 & FALSE & FALSE\\
47 & 9,17-; 9,15-diMeC27 & 0.0027638 & FALSE & FALSE\\
59 & 11,17-diMeC29 & 0.0032268 & FALSE & FALSE\\
66 & 14-; 13-; 12-; 11-MeC30 & 0.0037694 & FALSE & FALSE\\
48 & 9,13-diMeC27 & 0.0051902 & FALSE & FALSE\\
\addlinespace
28 & 5,13-diMeC25 & 0.0054108 & FALSE & FALSE\\
9 & 9,13-diMeC23 & 0.0059293 & FALSE & FALSE\\
4 & C23-ene & 0.0060814 & FALSE & FALSE\\
80 & 11,21-; 11,19-diMeC33 & 0.0091259 & FALSE & FALSE\\
36 & 10,16-; 10,14-diMeC26 & 0.0094887 & FALSE & FALSE\\
\addlinespace
52 & 3,9-; 3,7-diMeC27 & 0.0095619 & FALSE & FALSE\\
54 & 8,12,16-triMeC28 & 0.0097784 & FALSE & FALSE\\
18 & C25-ene & 0.0098478 & FALSE & FALSE\\
72 & 13,19-;13,17-diMeC31 & 0.0098522 & FALSE & FALSE\\
30 & C26 & 0.0103193 & FALSE & FALSE\\
\addlinespace
62 & 5,9-; 5,11-; 5,13-diMeC29 & 0.0103978 & FALSE & FALSE\\
75 & x-MeC32 & 0.0106378 & FALSE & TRUE\\
79 & x-MeC32 & 0.0131912 & FALSE & TRUE\\
46 & 11,15-diMeC27 & 0.0143588 & FALSE & TRUE\\
74 & (di)MeC31, mix & 0.0145589 & FALSE & TRUE\\
\addlinespace
70 & 15-; 13-; 11-; 9-MeC31 & 0.0150606 & FALSE & TRUE\\
45 & 5-MeC27 & 0.0154463 & FALSE & TRUE\\
35 & 4-MeC26 & 0.0156798 & FALSE & TRUE\\
76 & 13,17-diMe\_C32 + x,y-diMeC32 & 0.0157185 & FALSE & TRUE\\
17 & C25-ene & 0.0177937 & FALSE & TRUE\\
\addlinespace
71 & 11,19-; 11,17-; 11,15-diMeC31 & 0.0180838 & FALSE & TRUE\\
41 & C27 & 0.0196690 & FALSE & TRUE\\
38 & C27-ene + 6,12-diMeC26 & 0.0222253 & FALSE & TRUE\\
67 & (di)MeC30 (mix) & 0.0296362 & FALSE & TRUE\\
25 & 5-MeC25 & 0.0309633 & FALSE & TRUE\\
\addlinespace
50 & C28-ene & 0.0333263 & FALSE & TRUE\\*
\end{longtable}


## Species Identity Score

Ants in the mixed colonies exchanged their CHC to produce a hybrid chemical identity. We developed a method to assign each sample a value measuring similarity to one or the other species. It is based on a discriminant model that predicts the class of new samples using `predict` method from `mixOmics` R package. The input data are processed in the same way as the training set, i.e. central-log transformation is followed by projection into principal component space. In the latter step, the same rotation matrix and scaling vector is applied as those used for the training set are applied. This is achieved by using `predict` method on the object of `prcomp` class.



### Predicted species of different categories of ants as a function of the proportion of *F. sanguinea* ants in a colony

The Species Identity Score was computed for *F. sanguinea* and *F. fusca* samples from mixed colonies and used as a response variable in linear mixed models. As a fixed effect we used the proportion of *F. sanguinea* workers among all workers in a colony. Colony ID and sampling occasion were incorporated as random effects. 

In the following model specifications $(1|\text{random_factor})$ denotes random intercept (separate for each level of the random factor) and $((1|\text{random_factor_1:random_factor_2}))$ denotes random intercepts generated from the combinations of two factors. $\beta_{0}$ and $\epsilon$ denote intercept and error term, respectively. Included are summary reports, results of the Shapiro-Wilk test of nornalizty of conditional residuals, as well as diagnostic plots and tests generated with the use of `DHARMa` R package[@R-DHARMa]. The response variable often needed to be transformed to meet model assumptions. Random terms with no variance have been dropped. 



#### Mature *F. sanguinea*
$$
\sqrt{\text{Species_Identity_Score+1}} = \beta_0 + \text{sanguinea_proportion} + (1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{Species_Identity_Score}$ denotes the Species Identity Score of the mature *F. sanguinea* workers.

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: I(sqrt(predicted_species + 1)) ~ sang_prop + (1 | colony:census_date)
##    Data: model_input
## 
## REML criterion at convergence: 22.6
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.39156 -0.43663  0.06261  0.57841  1.64975 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.01917  0.1385  
##  Residual                       0.06029  0.2455  
## Number of obs: 62, groups:  colony:census_date, 42
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  0.96117    0.07083 42.22312  13.570  < 2e-16 ***
## sang_prop    1.00289    0.13239 39.12506   7.575 3.46e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.841
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.9752, p-value = 0.2413
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-18-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-18-2} \end{center}

#### *F. fusca*
$$
\text{Species_Identity_Score}= \beta_0 + \text{sanguinea_proportion} + (1|\text{colony})+(1|\text{colony_ID:sampling_occasion_ID}) + \epsilon,
$$
where $\text{Species_Identity_Score}$ denotes the Species Identity Score of the mature *F. fusca* workers.

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: predicted_species ~ sang_prop + (1 | colony) + (1 | colony:census_date)
##    Data: model_input
## 
## REML criterion at convergence: 240.5
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.73942 -0.41900 -0.07327  0.52549  2.41878 
## 
## Random effects:
##  Groups             Name        Variance Std.Dev.
##  colony:census_date (Intercept) 0.2073   0.4553  
##  colony             (Intercept) 0.1479   0.3846  
##  Residual                       0.9773   0.9886  
## Number of obs: 78, groups:  colony:census_date, 55; colony, 20
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)  -1.5743     0.2110 42.6506   -7.46 2.93e-09 ***
## sang_prop     3.3945     0.4606 36.1850    7.37 1.05e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.661
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.9847, p-value = 0.4753
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-19-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-19-2} \end{center}

#### Callow *F. sanguinea*
$$
\text{Species_Identity_Score}= \beta_0 + \text{sanguinea_proportion} + (1|\text{colony}) + \epsilon,
$$
where $\text{Species_Identity_Score}$ denotes the Species Identity Score of the callow *F. sanguinea* workers.

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
## Formula: predicted_species ~ sang_prop + (1 | colony)
##    Data: model_input
## 
## REML criterion at convergence: 76.3
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.25648 -0.51742 -0.08274  0.69855  1.88538 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  colony   (Intercept) 0.01975  0.1406  
##  Residual             0.62134  0.7883  
## Number of obs: 32, groups:  colony, 16
## 
## Fixed effects:
##             Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)  -0.9408     0.2530 28.3318  -3.719 0.000877 ***
## sang_prop     2.5681     0.4727 29.6945   5.433 7.07e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##           (Intr)
## sang_prop -0.821
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.97706, p-value = 0.7106
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-20-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-20-2} \end{center}

### Difference in Species Identity Score between mature *F. sanguinea* ants and *F. fusca* slaves

Since the linear model does not meet one of the diagnostics criteria, Wilcoxon paired test was applied.


```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  model_input$predicted_species and model_input$slave_species_index
## V = 1910, p-value = 4.507e-08
## alternative hypothesis: true location shift is not equal to 0
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-21-1} \end{center}

### Difference in Species Identity Score between callow *F. sanguinea* ants and *F. fusca* slaves
$$
\text{SIS_diff}= \beta_0 + \text{sanguinea_proportion} + (1|\text{colony}) + \epsilon,
$$
where $\text{SIS_diff}$ denotes the difference in Species Identity Score between the mature *F. sanguinea* and *F. fusca* sampled at the same time and from the same colony. Replicated samples were averaged before the calculations.

Since the linear model does not meet one of diagnostics criterion we will apply Wilcoxon paired test.

```
## 
## Call:
## lm(formula = index_diff ~ sang_prop, data = model_input)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.1745 -0.6534 -0.1962  0.7732  1.9023 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   0.3914     0.3492   1.121    0.272
## sang_prop    -0.4235     0.7201  -0.588    0.561
## 
## Residual standard error: 1.065 on 28 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.0122,	Adjusted R-squared:  -0.02308 
## F-statistic: 0.3458 on 1 and 28 DF,  p-value: 0.5612
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.97627, p-value = 0.7202
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-22-1} \end{center}

### Difference in Species Identity Score between callow and mature *F. sanguinea* ants.

```
## 
## Call:
## lm(formula = index_diff ~ sang_prop, data = model_input)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.3512 -0.4805 -0.1439  0.6242  2.1423 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept) -0.92662    0.33558  -2.761     0.01 *
## sang_prop    0.07754    0.61392   0.126     0.90  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.002 on 28 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.0005694,	Adjusted R-squared:  -0.03512 
## F-statistic: 0.01595 on 1 and 28 DF,  p-value: 0.9004
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  residuals(lm_model)
## W = 0.99264, p-value = 0.9987
```



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-23-1} \end{center}

<!--chapter:end:Species_identity.Rmd-->

# Finding peaks differentiating callow and mature *F. sanguinea* workers

The procedure identifying species markers was also applied to determine peaks characteristic of *F. sanguinea* callow ants. In discriminant analysis, samples were classified into one of the three groups: callow *F. sanguinea*,  adult *F. sanguinea*, or adult *F. fusca*.

## Marker indentification



\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-26-1} 

}

\caption{Distribution of the samples projected onto two first latent components before optimizing the number of model paramaters.}(\#fig:unnamed-chunk-26)
\end{figure}



\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-28-1} 

}

\caption{Projection of the samples in onto two first discriminant analysis components after model tuning.}(\#fig:unnamed-chunk-28)
\end{figure}



\begin{longtable}[t]{rlrlrl}
\caption{(\#tab:unnamed-chunk-30)List of the peaks with the score indicating their importance as a marker of mature \textit{F. fusca}, mature \textit{F. sanguinea}, or callow \textit{F. sanguinea} samples.}\\
\toprule
Peak ID & Compound & Callow score & Callow marker & Mature score & Mature marker\\
\midrule
\endfirsthead
\caption[]{(\#tab:unnamed-chunk-30)List of the peaks with the score indicating their importance as a marker of mature \textit{F. fusca}, mature \textit{F. sanguinea}, or callow \textit{F. sanguinea} samples. \textit{(continued)}}\\
\toprule
Peak ID & Compound & Callow score & Callow marker & Mature score & Mature marker\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
31 & 3,13-diMeC25 & 0.0416903 & TRUE & -0.0178266 & FALSE\\
24 & 7-MeC25 & 0.0386839 & TRUE & -0.0327207 & FALSE\\
28 & 5,13-diMeC25 & 0.0348565 & TRUE & -0.0317241 & FALSE\\
43 & 13-MeC27 & 0.0305096 & TRUE & -0.0234311 & FALSE\\
32 & 14-; 10-MeC26 + 3,7,11-triMeC25 & 0.0296012 & TRUE & -0.0373441 & FALSE\\
\addlinespace
68 & C31ene; 3-MeC30 & 0.0295308 & TRUE & -0.0189449 & FALSE\\
58 & 15-; 13-; 11-; 9-; 7-MeC29 & 0.0264321 & TRUE & -0.0072612 & FALSE\\
53 & 14-; 12-; 10-MeC28 & 0.0244925 & TRUE & -0.0172651 & FALSE\\
22 & 13-; 11-; 9-MeC25 & 0.0241580 & TRUE & -0.0457490 & FALSE\\
13 & 12-; 11-; 10-; 9-; 8-MeC24 & 0.0235367 & TRUE & -0.0429234 & FALSE\\
\addlinespace
10 & 3-MeC23 & 0.0233218 & TRUE & -0.0332150 & FALSE\\
12 & 3,13-; 3,11-; 3,9-; 3,7-diMeC23 & 0.0217166 & TRUE & -0.0258548 & FALSE\\
52 & 3,9-; 3,7-diMeC27 & 0.0207065 & FALSE & -0.0121810 & FALSE\\
27 & 7,11-diMeC25 + 3-MeC25 & 0.0199437 & FALSE & -0.0195216 & FALSE\\
38 & C27-ene + 6,12-diMeC26 & 0.0184220 & FALSE & 0.0157725 & FALSE\\
\addlinespace
19 & 4,12-diMeC24 & 0.0174162 & FALSE & -0.0198991 & FALSE\\
49 & 7,11-diMeC27 & 0.0161116 & FALSE & -0.0162587 & FALSE\\
26 & 11,15-; 9,13-diMeC25 & 0.0147629 & FALSE & -0.0161852 & FALSE\\
6 & 11-Me-C23 + 9-Me-C23 & 0.0145882 & FALSE & -0.0286965 & FALSE\\
23 & 9-MeC25 & 0.0127650 & FALSE & 0.0068551 & FALSE\\
\addlinespace
7 & 7-Me-C23 & 0.0126559 & FALSE & -0.0155144 & FALSE\\
47 & 9,17-; 9,15-diMeC27 & 0.0123952 & FALSE & -0.0148502 & FALSE\\
36 & 10,16-; 10,14-diMeC26 & 0.0111900 & FALSE & -0.0088143 & FALSE\\
17 & C25-ene & 0.0111205 & FALSE & 0.0258040 & TRUE\\
54 & 8,12,16-triMeC28 & 0.0073623 & FALSE & -0.0072445 & FALSE\\
\addlinespace
70 & 15-; 13-; 11-; 9-MeC31 & 0.0054970 & FALSE & 0.0234269 & TRUE\\
3 & C23-ene & 0.0044117 & FALSE & -0.0328012 & FALSE\\
5 & C23 & 0.0040366 & FALSE & -0.0384879 & FALSE\\
66 & 14-; 13-; 12-; 11-MeC30 & 0.0035717 & FALSE & 0.0179697 & FALSE\\
45 & 5-MeC27 & 0.0006637 & FALSE & 0.0242265 & TRUE\\
\addlinespace
72 & 13,19-;13,17-diMeC31 & 0.0002333 & FALSE & 0.0133734 & FALSE\\
60 & 13,17-diMeC29 & -0.0004252 & FALSE & 0.0173024 & FALSE\\
61 & 7,11-; 7,13-; 7,15-diMeC29 & -0.0024105 & FALSE & 0.0053526 & FALSE\\
44 & 7-MeC27 & -0.0043353 & FALSE & -0.0217227 & FALSE\\
75 & x-MeC32 & -0.0048799 & FALSE & 0.0234993 & TRUE\\
\addlinespace
77 & C33-ene & -0.0086063 & FALSE & 0.0217736 & FALSE\\
50 & C28-ene & -0.0108589 & FALSE & 0.0394054 & TRUE\\
59 & 11,17-diMeC29 & -0.0115280 & FALSE & 0.0183653 & FALSE\\
67 & (di)MeC30 (mix) & -0.0137091 & FALSE & 0.0378014 & TRUE\\
76 & 13,17-diMe\_C32 + x,y-diMeC32 & -0.0141628 & FALSE & 0.0414515 & TRUE\\
\addlinespace
37 & 7,x-; 5,x-; 10,14-diMeC26 & -0.0143395 & FALSE & -0.0207037 & FALSE\\
51 & C28 + 3,15-diMeC27 & -0.0144112 & FALSE & 0.0029394 & FALSE\\
79 & x-MeC32 & -0.0166740 & FALSE & 0.0323524 & TRUE\\
48 & 9,13-diMeC27 & -0.0187604 & FALSE & 0.0029709 & FALSE\\
46 & 11,15-diMeC27 & -0.0188128 & FALSE & 0.0345456 & TRUE\\
\addlinespace
11 & C24 & -0.0194179 & FALSE & -0.0262966 & FALSE\\
78 & x-MeC32 & -0.0200349 & FALSE & 0.0050400 & FALSE\\
56 & C29 & -0.0272272 & FALSE & -0.0143242 & FALSE\\
64 & C30 & -0.0274913 & FALSE & -0.0042224 & FALSE\\
80 & 11,21-; 11,19-diMeC33 & -0.0300790 & FALSE & 0.0357923 & TRUE\\
\addlinespace
71 & 11,19-; 11,17-; 11,15-diMeC31 & -0.0324225 & FALSE & 0.0406780 & TRUE\\
20 & C25 & -0.0341250 & FALSE & -0.0184365 & FALSE\\
30 & C26 & -0.0342711 & FALSE & 0.0024356 & FALSE\\
25 & 5-MeC25 & -0.0360842 & FALSE & 0.0417524 & TRUE\\
69 & C31 & -0.0421607 & FALSE & 0.0033436 & FALSE\\
\addlinespace
41 & C27 & -0.0464208 & FALSE & 0.0009641 & FALSE\\*
\end{longtable}


## Verifying the performace of the discrimination model

The performance of the discriminant model was evaluated by examining the accuracy of its predictions in a cross-validation test, measured as the area under the receiver operating characteristic curve (AUROC). This procedure was repeated $10^3$ times to account for variance resulting from the random splitting of samples into training and validation sets. In a parallel analysis, the model was trained on data with permuted age status of *F. sanguinea* to assess whether the true assignment of samples influenced model performance. Results from both training approaches were paired, and the difference in AUROC was calculated each for pair. The empirical p-value was defined as the proportion of models in which the AUROC after permutation was equal to or greater than the AUROC based on the original data.



The p-value is calculated as proportion of differences with value equal to or less than zero.




\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-34-1} 

}

\caption{Distribution of the differences in the performance of the discriminant models trained on true and permuted data. Age status of \textit{F. sanguinea} ants were shuffled in the null variant.}(\#fig:unnamed-chunk-34)
\end{figure}

All differences are positive, so the *p*-value is less than 0.001 since the null distribution consists of 1000 values.

<!--chapter:end:Callow_markers.Rmd-->

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-36-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-36-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-37-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-37-2} \end{center}
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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-38-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-38-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-39-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-39-2} \end{center}
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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-40-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-40-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-41-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-41-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-42-1} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-43-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-43-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-44-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-44-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-45-1} \end{center}

<!--chapter:end:CHCs_over_time.Rmd-->

# Comparison of the amount and proportions of CHC between species and age categories

This section presents the results of Wilcoxon tests comparing features of CHC profiles between callow and mature *F. sanguinea* ants as well as *F. fusca* slaves. Samples from the same colony were averaged before calculation of the final statistics to account for their non-independence.  



## Difference in the total CHC amount between *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-47-1} \end{center}

## Difference in the total CHC amount between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-48-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. sanguinea* between mature *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-49-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-50-1} \end{center}

## Difference in the proportion of CHC characteristic of callow *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-51-1} \end{center}


## Difference in the amount of CHC characteristic of callow *F. sanguinea* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-52-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. fusca* between mature *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-53-1} \end{center}

## Difference in the proportion of CHC characteristic of *F. fusca* between mature and callow *F. sanguinea* 

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-54-1} \end{center}


## Difference in the proportion of CHC characteristic of *F. fusca* between callow *F. sanguinea* and *F. fusca*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-55-1} \end{center}

## Difference in the mass of *n*-alkanes between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-56-1} \end{center}

## Difference in the proprtion of *n*-alkanes between callow and mature *F. sanguinea*

\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-57-1} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-59-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-59-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-60-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-60-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-61-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-61-2} \end{center}

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



\begin{center}\includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-62-1} \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-62-2} \end{center}

<!--chapter:end:Separation_experiment.Rmd-->

# Impact of slaves on the development of the *F. sanguinea* CHC profile

We analyzed the impact of *F. fusca* slaves on the CHC profile of callow *F. sanguinea* ants. In this experiment, *F. sanguinea* ants were isolated in pairs before eclosion, and their CHC profiles were analyzed at one of four time points marking their adult age. We then tested whether callow *F. sanguinea* ants were chemically more similar to their slave relatives from free-living colonies. Unrelated *F. fusca* ants served as a background group."

## Analysis based on all samples


\begin{figure}

{\centering \includegraphics[width=0.8\linewidth]{Supplementary_materials_files/figure-latex/unnamed-chunk-64-1} 

}

\caption{Distribution of lifespans of ants used in the experiment. The precise age of individuals could not be determined because ants in pairs were not marked. Therefore, all possible combinations were considered to illustrate the potential range.}(\#fig:unnamed-chunk-64)
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
\centering
\caption{(\#tab:unnamed-chunk-73)Chemical distance of the isolated \textit{F. sanguinea} ant to \textit{F. fusca} ants realted and unrelated to slaves.}
\centering
\begin{tabular}[t]{l>{\raggedleft\arraybackslash}p{4cm}>{\raggedleft\arraybackslash}p{4cm}}
\toprule
Adult age & Chemical distance to slaves' relatives & Chemical distance to random colonies\\
\midrule
1-3 days & 0.3332 & 0.4718\\
8-10 days & 0.3518 & 0.4802\\
17-20 days & 0.4521 & 0.5291\\
35-40 days & 0.5729 & 0.6100\\
\bottomrule
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
\centering
\caption{(\#tab:unnamed-chunk-84)Chemical distance of the isolated \textit{F. sanguinea} ant to \textit{F. fusca} ants realted and unrelated to slaves. Isolated ant spun a silky envelope before pupal stage, which prevented CHC transfer from the enviroment.}
\centering
\begin{tabular}[t]{l>{\raggedleft\arraybackslash}p{4cm}>{\raggedleft\arraybackslash}p{4cm}}
\toprule
Adult age & Chemical distance to slaves' relatives & Chemical distance to random colonies\\
\midrule
1-3 days & 0.3473 & 0.4869\\
8-10 days & 0.3740 & 0.4932\\
17-20 days & 0.4643 & 0.5407\\
35-40 days & 0.6175 & 0.6458\\
\bottomrule
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
\centering
\caption{(\#tab:unnamed-chunk-88)Chemical distance of separated \textit{F. sanguinea} ants to the CHC profile of \textit{F. fusca} ants from colonies that served as a source of CHC to caot the glass beads. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{ll>{\raggedleft\arraybackslash}p{2cm}l>{\raggedleft\arraybackslash}p{2cm}}
\toprule
Colony ID & Treamtent & Averaged distance & Contrast & Averaged distance\\
\midrule
SD18-3 & 12-15 days + F. fusca hydrocarbons & 0.396 & 12-15 days (control) & 0.509\\
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.616 & 12-15 days (control) & 0.513\\
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.842 & 12-15 days (control) & 0.832\\
SD19-11 & 12-15 days + F. fusca hydrocarbons & 0.467 & 12-15 days (control) & 0.469\\
SD19-3 & 12-15 days + F. fusca hydrocarbons & 0.500 & 12-15 days (control) & 0.436\\
\addlinespace
SD19-4 & 12-15 days + F. fusca hydrocarbons & 0.305 & 12-15 days (control) & 0.424\\
SD19-6 & 12-15 days + F. fusca hydrocarbons & 0.674 & 12-15 days (control) & 0.699\\
SD19-8 & 12-15 days + F. fusca hydrocarbons & 0.595 & 12-15 days (control) & 0.587\\
SD20-2 & 12-15 days + F. fusca hydrocarbons & 0.699 & 12-15 days (control) & 0.709\\
\bottomrule
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
\centering
\caption{(\#tab:unnamed-chunk-89)Chemical distance of separated \textit{F. sanguinea} ants to the CHC profile of \textit{F. sanguinea} ants from colonies that served as a source of CHC to caot the glass beads. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{ll>{\raggedleft\arraybackslash}p{2cm}l>{\raggedleft\arraybackslash}p{2cm}}
\toprule
Colony ID & Treamtent & Averaged distance & Contrast & Averaged distance\\
\midrule
SD18-2 & 12-15 days + F. sanguinea hydrocarbons & 0.571 & 12-15 days (control) & 0.373\\
SD18-3 & 12-15 days + F. sanguinea hydrocarbons & 0.457 & 12-15 days (control) & 0.599\\
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.525 & 12-15 days (control) & 0.580\\
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.250 & 12-15 days (control) & 0.232\\
SD19-11 & 12-15 days + F. sanguinea hydrocarbons & 0.595 & 12-15 days (control) & 0.607\\
\addlinespace
SD19-3 & 12-15 days + F. sanguinea hydrocarbons & 0.419 & 12-15 days (control) & 0.511\\
SD19-4 & 12-15 days + F. sanguinea hydrocarbons & 0.648 & 12-15 days (control) & 0.668\\
SD19-6 & 12-15 days + F. sanguinea hydrocarbons & 0.553 & 12-15 days (control) & 0.591\\
SD19-8 & 12-15 days + F. sanguinea hydrocarbons & 0.850 & 12-15 days (control) & 0.850\\
SD20-2 & 12-15 days + F. sanguinea hydrocarbons & 0.646 & 12-15 days (control) & 0.654\\
\bottomrule
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
\centering
\caption{(\#tab:unnamed-chunk-90)Proportion of \textit{n}-docosane in CHC extracted from \textit{F. sanguinea} ants maintained with the glass beads coated with the CHC of \textit{F. fusca} ants and contaminated with \textit{n}-docosane. In control variant, glass bead were left clean.}
\centering
\begin{tabular}[t]{ll>{\raggedleft\arraybackslash}p{2cm}l>{\raggedleft\arraybackslash}p{2cm}}
\toprule
Colony ID & Treamtent & Averaged C22 fraction & Contrast & Averaged C22 fraction\\
\midrule
SD18-3 & 12-15 days + F. fusca hydrocarbons & 0.018 & 12-15 days (control) & 0.000\\
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.001 & 12-15 days (control) & 0.000\\
SD18-5 & 12-15 days + F. fusca hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
SD19-11 & 12-15 days + F. fusca hydrocarbons & 0.011 & 12-15 days (control) & 0.002\\
SD19-3 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.001\\
\addlinespace
SD19-4 & 12-15 days + F. fusca hydrocarbons & 0.017 & 12-15 days (control) & 0.000\\
SD19-6 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.001\\
SD19-8 & 12-15 days + F. fusca hydrocarbons & 0.003 & 12-15 days (control) & 0.002\\
SD20-2 & 12-15 days + F. fusca hydrocarbons & 0.009 & 12-15 days (control) & 0.005\\
W17-1 & 12-15 days + F. fusca hydrocarbons & 0.021 & 12-15 days (control) & 0.001\\
\bottomrule
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


\begin{longtable}[t]{ll>{\raggedleft\arraybackslash}p{2cm}l>{\raggedleft\arraybackslash}p{2cm}}
\caption{(\#tab:unnamed-chunk-91)Proportion of \textit{n}-docosane in CHC extracted from \textit{F. sanguinea} ants maintained with the glass beads coated with the CHC of \textit{F. sanguinea} ants and contaminated with \textit{n}-docosane. In control variant, glass bead were left clean.}\\
\toprule
Colony ID & Treamtent & Averaged C22 fraction & Contrast & Averaged C22 fraction\\
\midrule
SD18-2 & 12-15 days + F. sanguinea hydrocarbons & 0.005 & 12-15 days (control) & 0.002\\
SD18-3 & 12-15 days + F. sanguinea hydrocarbons & 0.014 & 12-15 days (control) & 0.000\\
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.000\\
SD18-5 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
SD19-11 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.002\\
\addlinespace
SD19-3 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.001\\
SD19-4 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.000\\
SD19-6 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.001\\
SD19-8 & 12-15 days + F. sanguinea hydrocarbons & 0.002 & 12-15 days (control) & 0.002\\
SD20-2 & 12-15 days + F. sanguinea hydrocarbons & 0.004 & 12-15 days (control) & 0.005\\
\addlinespace
W17-1 & 12-15 days + F. sanguinea hydrocarbons & 0.003 & 12-15 days (control) & 0.001\\
\bottomrule
\end{longtable}

<!--chapter:end:Dummy_ants.Rmd-->

# Body surface area of *F. sanguinea* ants and their slaves {#surface-area}

The areas of planar projection of different body parts were regressed on the head width to obtain a model which was subsequently used to approximate the area for each sample. The regression was performed separately for each species. The following formula was applied:

$$
\hat{A} = \sum_{i \in I}\beta_{0,i}+\beta_{1,i}w+\beta_{2,i}w^2,
$$
where\newline
$\hat{A}$ denotes the estimated body area,\newline
$I$ denotes the set of body parts for which the planar projection area was measured\newline
$w$ denotes head width\newline
$\beta_{j,i}$ are coefficients from the polynomial regression\newline

These calculations do not reflect the absolute value of body surface area but were converted into scaling factors by using them a divisors of a standard area. In this way, the unit of measurement is canceled out, and absolute values are converted into dimensionless ratios. The standard area was calculated by substituting the $w$ = 1.15 mm into the formula above and using coefficients obtained from measurements of *F. sanguinea*. 


\begin{figure}

{\centering \includegraphics{Supplementary_materials_files/figure-latex/unnamed-chunk-93-1} 

}

\caption{Relation between head width and area of the planar projections of ant body parts.}(\#fig:unnamed-chunk-93)
\end{figure}

<!--chapter:end:Body_surface.Rmd-->

