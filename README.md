# Phylogenetic identification of influenza virus candidates for seasonal vaccines
This repository contains the codes and materials for the paper: "Phylogenetic identification of influenza virus candidates for
seasonal vaccines"

The seasonal influenza (flu) vaccine is designed to protect against those influenza viruses predicted to circulate during the upcoming flu season. Due to ongoing genetic drift, in particular in the hemagglutinin (HA) and neuraminidase (NA) genes, these vaccines must be updated annually. We use phylogenetic trees reconstructed from HA and NA sequences from human influenza viruses (using isolates from 1980-2020), together with counts of epitope site polymorphisms in hemagglutinin, to predict which influenza virus strains are likely to persist in the next season. We define a set of features, computed for each taxon, using the epitope polymorphisms and structures of the HA and NA phylogenies. We train a support vector machine to classify taxa as likely to be “successful” (circulate in the coming seasons) or not. We obtain accuracies of 0.75-0.89 and a classifier ‘area under the curve’ (AUC) of 0.83-0.91, testing for years 2016-2020. Among the set of successful strains in each year, we explore several ways to select potential candidates for inclusion in the seasonal vaccine. We compare strains’ proximity to relevant ”future” populations, and find that the machine learning model has a moderate ability to select strains that are close to future populations. However, consensus sequences among the most recent three years also do well at this task. We identify similar candidate strains to those proposed by the World Health Organization, suggesting that this approach can help inform vaccine strain selection

# code

This folder contains all the codes to compute the features, train and test the machine learing models, predict the next year's flu vaccine candidates, and generate the plots for the paper. 
 
"AUC_subtree.R" code to make the SVM models for each year.
"AllClassifiers.R" code to compare different binary classifiers on each year datasets.
"Data_5_1_3.R" code to compute the features.
"LBI_cor.R" code to calculate the correlation between LBI and diversificationRate statistic (tree-based).
"LBI_cor_tip.R" code to calculate the correlation between LBI and diversificationRate statistic. 
"Tree_Statistics.R" code to compute tree shape statistics.
"VaccineDist.R" code to compute the distance between our suggested candidates and the next year's circulating sequences.



# 2016-2020

These folders contain all the influenza trees and other data for the years 2016 to 2020 respectively. 

# General_data

This folder contains general data that is not specific to any particular year.  

# Figures

This folder contains all the figures of the paper.






