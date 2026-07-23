# Machine learning of Ising criticality with spin-shuffling

## Authors: Pallab Basu, Jyotirmoy Bhattacharya , Dileep Pavan Surya Jakka, Chuene Mosomane, Vishwanath Shukla. 
### email: jyotirmoy.phys@gmail.com

This repository contains codes which were used in the paper `Machine learning of Ising criticality with spin-shuffling' [https://arxiv.org/abs/2203.04012] 

In this project various elementary ANNs (elementary classifiers) were created to identify the phases and in particular the critical region of the 2D Ising model with NN as well as NNN interactions. 

# Data

High quality data were generated with the help of codes kept within the fodler titled `Data_generate_Code' in this repository. 

The Dataset which was directly used can be found on this public Kaggle repository - https://doi.org/10.34740/kaggle/dsv/2957834

# Codes

All the codes which were used to make the plots in the paper can be found within the folder "Codes". These codes were created to run on Kaggle. The initialization and Data file inputs reflect this. 

**final-v0-mag-chi-xi.ipynb** : Generates the figures 2,3 and 4 in the paper. It takes the generated data as input and computes the magnetization, susceptibility and Correlation length. 

**Critical_prediction**: This folder contains the codes which tests the prediction of the critical region. The plots in the Figures 6, 7 and 8 of the paper were constructed using the codes in this folder. 

**Phase Prediction**: The Phase prediction accuracies reported in tables II and III were obtained with the codes within this folder. 
