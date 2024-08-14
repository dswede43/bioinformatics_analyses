#!/usr/bin/python3

#Machine learning model testing
#---
#Script to test the various machine learning classifier models.


#Table of contents:
#---
#1. Data preparation
#2. Model testing
#3. Visualize ROC curves
#4. Visualize confusion matrices


#load libraries
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score


#define global variables
DIR = "/mnt/c/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
MODEL_NAMES = ['lasso_log','random_forest','AdaBoost']
TARGET = 'cell_group' #target variable name
SCALED_MODELS = ['lasso_log'] #list of model names that require scaled data inputs

#load data
test_data = pd.read_csv(f"{DIR}/Data/ML_test_data.csv") #testing dataset

#load models
models = {}
for model_name in MODEL_NAMES:
    with open(f"{DIR}/Models/{model_name}_model.pkl", 'rb') as file:
        models[model_name] = pickle.load(file)


#1. Data preparation
#---
#convert the target variable into 0s and 1s
test_data[TARGET] = test_data[TARGET].apply(lambda x: 1 if x == 'T_cell' else 0)

#define target and feature variables
features = [col for col in test_data.columns if TARGET not in col]
X_test = test_data[features].values
y_test = test_data[TARGET]

#scale the feature data
scaler = StandardScaler()
X_test_scaled = scaler.fit_transform(X_test)


#2. Model testing
#---
#calculate the ROC curve for a model with random performance
r_probs = [0 for _ in range(len(y_test))]
r_fpr, r_tpr, _ = roc_curve(y_test, r_probs)
r_auc = roc_auc_score(y_test, r_probs)

#define empty dictionaries
cm_data = {}
roc_data = {}
for model_name in models.keys():
    #define the current model
    model = models.get(model_name)
    
    #create model predictions
    if model_name in SCALED_MODELS:
        y_pred = model.predict(X_test_scaled)
    else:
        y_pred = model.predict(X_test)
    
    #create the confusion matrix
    cm_data[model_name] = confusion_matrix(y_test, y_pred)
    
    #calculate the predictive probabilities for each model
    if model_name in SCALED_MODELS:
        model_probs = model.predict_proba(X_test_scaled)[:, 1]
    else:
        model_probs = model.predict_proba(X_test)[:, 1]
    
    #caclulate the model true positive rate (TPR) and false positive rate (FPR)
    model_fpr, model_tpr, _ = roc_curve(y_test, model_probs)
    
    #calculate the model ROC AUC score
    model_auc = roc_auc_score(y_test, model_probs)
    
    #store the results
    roc_data[model_name] = [model_fpr, model_tpr, model_auc]


#3. Visualize ROC curves
#---
#create a subplot
fig, axes = plt.subplots(nrows = 1, ncols = len(roc_data), figsize = (12, 4))

#plot the ROC curves
for i, model_name in enumerate(roc_data.keys()):
    roc = roc_data.get(model_name)
    model_fpr = roc[0]
    model_tpr = roc[1]
    model_auc = roc[2]
    axes[i].plot(r_fpr, r_tpr, linestyle = '--', label = f"Random prediction (AUC: {r_auc:.2f})")
    axes[i].plot(model_fpr, model_tpr, marker = '.', label = f"{model_name} (AUC: {model_auc:.2f})")
    axes[i].set_title(f"{model_name} model")
    axes[i].set_xlabel('FPR (1 - specificity)')
    axes[i].set_ylabel('TPR (sensitivity)')
    axes[i].legend()

#save the plots as PDF
plt.tight_layout()
plt.savefig(f"{DIR}/Visualizations/ROC_curves.pdf")
plt.close()


#4. Visualize confusion matrices
#---
fig, axes = plt.subplots(nrows = 1, ncols = len(cm_data), figsize = (12, 4))

#plot the confusion matrix
for i, model_name in enumerate(cm_data.keys()):
    cm = cm_data.get(model_name)
    sns.heatmap(cm, ax = axes[i], annot = True, fmt = '.0f', cmap = 'Blues', cbar = False)
    axes[i].set_title(f"{model_name} model")
    axes[i].set_xticks([0.5, 1.5], ['non T-cell','T-cell'])
    axes[i].set_yticks([0.5, 1.5], ['non T-cell','T-cell'])
    axes[i].set_xlabel('Predicted label')
    axes[i].set_ylabel('True label')

#save plots as PDF
plt.tight_layout()
plt.savefig(f"{DIR}/Visualizations/confusion_matrices.pdf")
plt.close
