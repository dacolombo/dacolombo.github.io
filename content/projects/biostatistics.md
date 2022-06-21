---
title: "DES treatment in advanced prostate cancer"
date: 2022-05-12
draft: false
cover: "/biostatistics/cover.png"
description: "In this group project me and three of my collegues evaluated the efficacy of DES treatment for advanced prostate cancer, using different models to evaluate the success of the therapy depending on dosage and on clinical parameters" 
keywords: ["Biostatistics","Prostate Cancer","Clinical Trial","Survival Analysis","Cox Model","Logistic Regression","Competing Risk Analysis"]
tags: ["Biostatistics","Cancer","Survival","Regression","CompetingRiskAnalysis"]
----

In this project, we used the Byar & Greene prostate cancer data available from [Vanderbilt Biostatistics Datasets](https://hbiostat.org/data/).
This is a dataset coming from a randomised clinical trial with the aim of comparing survival outcome according to different treatments. In particular, the patients are divided into four treatment groups, each one receiving 0.2 mg DES, 1.0 mg DES, 5.0 mg DES and placebo respectively.
The dataset used contains data about survival time and clinical parameters of 506 patients suffering from prostate cancer in stages 3 and 4.

# Introduction
The prostate is a gland in the male reproductive system which produce and contain the fluid that forms part of semen. The semen-secreting gland cells can mutate into cancer cells which led to the prostate cancer.
Prostate cancer cells usually require androgen hormones, such as testosterone, to grow. One possible approach for prostate cancer treatment is the androgen suppression therapy, which reduces the level of these hormones.
The drug used in the study is the Diethylstilbestrol (DES), a synthetic ethinyl estrogen. Its anti-tumor efficacy is attributable to several actions on the prostate cancer:
* Negative feedback on hypothalamic-pituitary complex
* Increasing SHBG, hence decreasing TST
* Suppress leydig cells function
* Telomerase inhibitor of cancer cells

However, from some studies it is known that DES has some adverse effects, in particular cardiotoxicity.


# Data preprocessing
The preprocessing step consisted in variable categorization (e.g. death cause, age, blood pressure), removal of null values and of unnecessary varaibles.
After that we were left with data about 474 patients.

# Exposure-event correlation
Before doing an evaluation of the survival and the efficacy of the therapy, we evaluated the correlation between DES administration and cardiovascular death, using the Chi-squared test and odds ratios.
The Chi-squared test was performed comparing vascular death with other events in all the treatment groups, resulting in a significative p-value indicating a relationship between the treatment and cardiovascular events.
As expected, we also found a correlation between DES and cardiovascular events with the OR indicating a higher risk of cardiovascular event related to the higher dosage of DES, i.e. 5.0 mg

{{<image src="/biostatistics/odds_ratio.png" position="center">}}

# Survival analysis
We then started looking at the survival data, firstly by taking a qualitative look using the Kaplan-Meier curves.
As can be seen in the stratified plot, there seems to be a higher survival for patients treated with 1.0 mg of DES with respect to the other treatment groups which seems to have similar outcomes.

{{<image src="/biostatistics/kaplan_meier.png" position="center">}}

# Cox Model
Through the Cox Model we could evaluate relevant covariates for survival. Once fitted, we plotted the model parameters through a forest plot and we found both:
* Risk factors: age group between 76-89, positive history of cardiovascular disease, presence of bone methastasis, heart strain ECG reading, tumor size.
* Protective factors: treatment with 1.0 mg DES

{{<image src="/biostatistics/cox_forest.png" position="center">}}

After fitting the Cox model, it was possible to predict the theoretical survival of a patient based on its health parameters. The selected patients for this evaluation had the following parameters:
* Patient 1 (median patient): no history of CVD, normal ECG, median tumor size, median weight index, no bone methastasis, lower age group.
* Patient 2 (severe cancer and cardiovascular conditions): history of CVD, heart strain ECG reading, high tumor size, high weight index, bone methastasis, higher age group.
* Patient 3 (severe cancer conditions): no history of CVD, normal ECG, high tumor size, median weight index, bone methastasis, lower age group.
* Patient 4 (severe cardiovascular conditions): history of CVD, heart strain ECG reading, low tumor size, median weight index, no bone methastasis, lower age group.

As expected, patients with severe conditions showed a lower survival, especially the ones having both bad cancer and cardiovascular conditions. What was found of interest is also that in all patients, the 1.0 mg DES treatment performed better, with the other treatments performing similarly.

{{<image src="/biostatistics/adjusted_curves1.png" position="center">}}
{{<image src="/biostatistics/adjusted_curves2.png" position="center">}}

# Accelerated Failure Time (AFT) model
In order to get a confirmation of what was observed, we fitted an accelerated failure time model with the same covariates found to be relevant in the Cox model. The fitted model resulted in the survival time being distributed as a Weibull with coefficients (0.89, 1.12), with the AFT interpretation of the covariates being in line with the proportional hazard intepretation: variables found to increase the hazard were found to reduce the time of survival and vice versa.

{{<image src="/biostatistics/aft.png" position="center">}}

| Coefficient       | AFT interpretation | PH interpretation |
|:-----------------:|:------------------:|:-----------------:|
| hx=1              | 0.64               | 1.65              |
| DES 1.0 mg        | 1.36               | 0.71              |
| Ekg: heart strain | 0.67               | 1.55              |
| Wt                | 1.01               | 0.98              |
| Bm=1              | 0.65               | 1.61              |
| Agegrp=76-89      | 0.61               | 1.74              |

# Logistic regression model
We then used logistic regression to model the survival outcome based on the variables related to the patients. We managed to get a model with a sufficiently high AUC=0.72 and by using a sperimental threshold of 0.35 we obtained specificity and sensitivity at around 65%.
Also in this case, the results were concordant with what was observed before: treatment with 5.0 mg DES and a history of CVD were found to significantly increase the odds of the event.

{{<image src="/biostatistics/logistic_regression.png" position="center">}}

| Threshold           | MIS rate | Sensitivity | Specificity |
|:-------------------:|:--------:|:-----------:|:-----------:|
|  p = 0.5            | 26.6%    |  23.2%      |  94%        |
| Sperimental = 0.348 | 34.8%    | 68.1%       | 63.9%       |

| Variable   | Odds                  |
|:----------:|:---------------------:|
| DES 5.0 mg | 1.839 [1.036 , 3.266] |
| Hx = 1     | 4.081 [2.645 , 6.298] |

# Competing risk analysis
To complete our analysis, we decided to use a more advanced model to correctly estimate marginal probability of an event in the presence of competing events for the same outcome. By considering cardiovascular events and cancer complications events as two competing risks for death, we evaluated how each covariate affected the two events for the patient.
As can be seen in the table below, some interesting results can be evaluated from the competing risk analysis (note that hazard ratios significantly different from 1 are highlighted).
First of all, it is evident that 1.0 mg of DES is a protective factor for death by prostate cancer, while 5.0 mg of DES is a risk factor for death by cardiovascular events. Moreover, as expected, we can see that cardiovascular-related covariates (e.g. problematic ECG) have an influence on cardiovascular events, while tumor-related covariates (e.g. tumore size, bone methastasis) have an influence on death by prostate cancer.

{{<image src="/biostatistics/car_coefficients.png" position="center">}}

With this model we could then evaluate the cumulative incidence for each death cause stratified by treatment. With this plots, it's possible to see that the probability of a cardiovascular event is highly increased with 5.0 mg of DES, while 1.0 mg of DES decreases the probability of death due to prostate cancer with respect to both other doses and placebo.

{{<image src="/biostatistics/car_plots.png" position="center">}}

# Conclusions
We evaluated the presence of correlation between DES exposure and cardiovascular events, observing that 5.0 mg of DES are positively related to death caused by cardiovascular events.
Through the Cox PH and AFT model we evaluated the significant variables that inluence the survival of patients. Through these variables, we found that the optimal dose of DES is 1.0 mg as it significantly reduces the hazard.
With logistic regression, we could build a model to classify a patient at risk for cardiovascular events, confirming that the treatment with 5.0 mg of DES are positively related to cardiovascular events.
Lastly, through competing risk analysis we could observe the differences of the impact of variables on the different events of interest: death due to cardiovascular events and due to prostate cancer. Thanks to this analysis, we found that some variables not relevant in the standard Cox PH model becomes relevant for the event-specific model, like ECG read, cancer size and treatment with 5.0 mg DES. With this model, we could more percisely assess the parameters to consider and the best possible treatment.