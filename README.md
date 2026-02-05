# Evolutionary Markov Chain for Adverse Drug Reaction

R package implementing an Evolutionary Monte-Carlo Markov Chain algorithm (an adaptation of Metropolis-Hastings). This package is designed to be used with medical data, especially with patients using medications.

*Supervisor : Mr. Birmele Etienne*

## Install the package via CRAN repository

You can simply install the package by using the following command:

```
install.packages('emcAr')
```

## Get the modified ATC tree (containing every medications)

The algorithm uses a modified medication tree which include an upper bound that locates the last drug in the drug family represented by the current node (if the current node is not a leaf). You can find the original drugs tree in the ```emcAdr/data``` folder. You can also use your own tree but an upper bound for each node is **mandatory** (the upper bound of a leaf is the index of this leaf in your 2D array).

## Build or use a dataset of patients

The algorithm requires a data frame of patient. Every line of this Data frame represents a single patient, the medications they are taking and a boolean representing whether they have the adverse drug effect under consideration.

They are 2 mandatory columns : *patientATC* and *patient ADR*. Respectively the index of the drugs taken by the patient in the tree of drugs (indexes start at 0) and the boolean representing whether he has the ADR.

There is an example of a row for a patient who takes 3 drugs (having 12, 56 and 798 as indexes) and doesn't have the adverse drug reaction under consideration : 

| patientATC | patientADR	|
|---	       |---	        |
|12, 56, 798 |      0    	|

## Use the EMC function 

You are now ready to use the EMC function contained in the package. Here is an example 
```
res <- EMC(n = 100,nbIndividuals = 5,ATCtree = ATC_Tree_UpperBound_2014, observations = simulPatient_df, startingIndividuals = c(), startingTemperatures = c())
```
