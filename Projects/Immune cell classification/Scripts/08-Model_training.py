#!/usr/bin/python3

#Machine learning model training
#---
#Script to train various machine learning classifier models.


#Table of contents:
#---
#1. Data preparation
#2. Model training
    #a) LASSO logistic regression
    #b) random forest
    #c) adaptive boost


#load libraries
import pandas as pd
import pickle
import time

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score


#define global variables
DIR = "/mnt/c/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
TARGET = 'cell_group' #target variable name
CV = 5 #k-folds cross-validation
SCALED_MODELS = ['lasso_log'] #list of model names that require scaled data inputs
ENSEMBLE_MODELS = ['random_forest','AdaBoost'] #list of ensemble models that use feature importance

#load data
train_data = pd.read_csv(f"{DIR}/Data/ML_train_data.csv") #training dataset
test_data = pd.read_csv(f"{DIR}/Data/ML_test_data.csv") #testing dataset
gene_names = pd.read_csv(f"{DIR}/Data/gene_biotypes-names.csv") #gene names and biotypes


#1. Data preparation
#---
#convert the target variable into 0s and 1s
train_data[TARGET] = train_data[TARGET].apply(lambda x: 1 if x == 'T_cell' else 0)
test_data[TARGET] = test_data[TARGET].apply(lambda x: 1 if x == 'T_cell' else 0)

#define target and feature variables
features = [col for col in train_data.columns if TARGET not in col]
X_train = train_data[features].values
X_test = test_data[features].values
y_train = train_data[TARGET]
y_test = test_data[TARGET]

#scale the feature data
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.fit_transform(X_test)


#2. Model training
#---
#define the base estimator for AdaBoost
base_estimator = DecisionTreeClassifier(max_depth = 1)

#define the models
models = {
    'lasso_log': LogisticRegression( #LASSO logistic regression
        penalty = 'l1',
        solver = 'liblinear',
        random_state = 42),
    'random_forest': RandomForestClassifier( #random forest
        random_state = 42),
    'AdaBoost': AdaBoostClassifier( #adaptive boost
        estimator = base_estimator,
        algorithm = 'SAMME',
        random_state = 42)
    }

#define the model parameter grids to fine-tune
param_grids = {
    'lasso_log': {
        'C': [0.01, 0.1, 1, 10, 100] #regularization penalty term
        },
    'random_forest': {
        'n_estimators': [50, 100, 200, 300], #number of trees
        'max_depth': [10, 20, 30, None], #max depth of each tree
        'min_samples_split': [2, 5, 10], #min number of samples required to split a node
        'min_samples_leaf': [1, 2, 4], #min number of samples required at a leaf node
        'max_features': ['sqrt', 'log2', None], #max number of features considered for splitting a node
        'bootstrap': [True, False] #bootstrap samples when building trees
        },
    'AdaBoost': {
        'n_estimators': [50, 100, 200], #number of trees
        'learning_rate': [0.01, 0.1, 1] #learning rate (scaling of weak learner)
        }
    }

#define an empty data frame
used_genes = pd.DataFrame(columns = ['ensembl','model_name'])
for model_name in models.keys():
    #define the current model and parameter grid space
    model = models.get(model_name)
    param_grid = param_grids.get(model_name)

    #set the grid search space for fine-tuning parameters
    grid_search = GridSearchCV(
        estimator = model,
        param_grid = param_grid,
        cv = CV,
        scoring = 'accuracy',
        n_jobs = -1)

    #fit the model using the training data
    start_time = time.time()
    print(f"Fine-tuning parameters for {model_name} model...")
    if model_name in SCALED_MODELS:
        grid_search.fit(X_train_scaled, y_train)
    else:
        grid_search.fit(X_train, y_train)

    #print the training time
    end_time = time.time()
    train_time = end_time - start_time
    print(f"Training time for {model_name} model: {train_time} seconds")

    #print the optimal parameters
    optimal_params = grid_search.best_params_
    print(f"Optimal parameters for {model_name} model: {optimal_params}")

    #define the most optimal model
    best_model = grid_search.best_estimator_

    #determine the model accuracy
    if model_name in SCALED_MODELS:
        y_pred = best_model.predict(X_test_scaled)
    else:
        y_pred = best_model.predict(X_test)

    #print the model accuracy
    model_accuracy = accuracy_score(y_test, y_pred)
    print(f"Accuracy of {model_name} model: {model_accuracy}")

    #obtain the features used in the model
    if model_name in ENSEMBLE_MODELS:
        #obtain model feature importances
        importances = best_model.feature_importances_

        #obtain the used model features
        used_features = [feature for feature, importance in zip(features, importances) if importance != 0]
        print(f"Number of gene features used by {model_name} model: {len(used_features)}")
    else:
        #obtain model feature coefficients
        coefs = best_model.coef_[0]

        #obtain used model features
        used_features = [feature for feature, coef in zip(features, coefs) if coef != 0]
        print(f"Number of gene features used by {model_name} model: {len(used_features)}")

    #store the results
    used_features = pd.DataFrame(used_features)
    used_features['model_name'] = model_name
    used_features.columns = ['ensembl','model_name']
    used_genes = pd.concat([used_genes, used_features], ignore_index = True)

    #save model as PKL
    with open(f"{DIR}/Models/{model_name}_model.pkl", 'wb') as file:
        pickle.dump(best_model, file)

#obtain gene names and biotypes
used_genes = pd.merge(gene_names, used_genes, on = 'ensembl', how = 'inner')

#save results as CSV
used_genes.to_csv(f"{DIR}/Data/used_genes.csv", index = False)
