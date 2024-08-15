#!/usr/bin/python3

#Machine learning model visualizations
#---
#Script to visualize the decision boundaries for various machine learning classifier models.


#Table of contents:
#---
#1. Data preparation
#2. Plotting model decision-boundaries


#load libraries
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#define global variables
DIR = "/mnt/c/Users/dstra/Nextcloud/Documents/Computers/Code/Github repos/bioinformatics analyses/Machine learning/Immune cell classification" #working directory
TARGET = 'cell_group' #target variable name
MODEL_NAMES = ['lasso_log','random_forest','AdaBoost']
SCALED_MODELS = ['lasso_log'] #list of model names that require scaled data inputs
CMAP = plt.cm.RdYlBu #color mapping for plots

#load data
test_data = pd.read_csv(f"{DIR}/Data/ML_test_data.csv") #testing dataset

#load models
models = {}
for model_name in MODEL_NAMES:
    with open(f"{DIR}/Models/{model_name}_model.pkl", 'rb') as file:
        models[model_name] = pickle.load(file)


#1. Data preparation
#---
#define the target variable levels
groups = test_data[TARGET]

#convert the target variable into 0s and 1s
test_data[TARGET] = test_data[TARGET].apply(lambda x: 1 if x == 'T_cell' else 0)

#define target and feature variables
features = [col for col in test_data.columns if TARGET not in col]
X_test = test_data[features].values
y_test = test_data[TARGET]

#scale the feature data
scaler = StandardScaler()
scaler.fit(X_test)
X_test_scaled = scaler.transform(X_test)

#calculate PC scores
pca = PCA(n_components = 2)
X_test_pca = pca.fit_transform(X_test_scaled)

#find the variance explained by each PC
var = pca.explained_variance_ratio_


#2. Plotting model decision-boundaries
#---
#define function to create the decision boundary for a model
def create_decision_boundaries(X, model, pca, scaler, scaled_features = False, plot_step = 0.2):
    #define the grid space for plotting
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                         np.arange(y_min, y_max, plot_step))

    #define the feature data for the grid space
    if scaled_features == True:
        grid_features = pca.inverse_transform(np.c_[xx.ravel(), yy.ravel()])
    else:
        grid_features = pca.inverse_transform(np.c_[xx.ravel(), yy.ravel()])
        grid_features = scaler.inverse_transform(grid_features)

    #make model predictions for the grid space
    Z = model.predict(grid_features)
    Z = Z.reshape(xx.shape)

    #return variables
    return Z, xx, yy

#define an empty dictionary
grid_space = {}
for model_name in models.keys():
    #define the current model
    model = models.get(model_name)

    #create the decision boundary grid space
    if model_name in SCALED_MODELS:
        Z, xx, yy = create_decision_boundaries(X_test_pca, model, pca, scaler, scaled_features = True)
    else:
        Z, xx, yy = create_decision_boundaries(X_test_pca, model, pca, scaler)

    #store the results
    grid_space[model_name] = [Z, xx, yy]

#create legend
color_map = {'non T-cell': 'red', 'T-cell': 'blue'}
handles = [plt.Line2D([0], [0], marker = 'o', color='w', markerfacecolor = color_map[cat], markersize = 10, label = cat) for cat in color_map]

#create a subplot
fig, axes = plt.subplots(nrows = 1, ncols = len(grid_space), figsize = (10, 4))

#plot the model decision bounderies
for i, model_name in enumerate(grid_space.keys()):
    grid_data = grid_space.get(model_name)
    Z = grid_data[0]
    xx = grid_data[1]
    yy = grid_data[2]
    axes[i].contourf(xx, yy, Z, alpha=0.3, cmap = CMAP)
    axes[i].scatter(X_test_pca[:, 0], X_test_pca[:, 1], c = y_test, alpha = 0.8, edgecolor = 'k', s = 25, cmap = CMAP)
    axes[i].set_title(f"{model_name} decision boundary")
    axes[i].set_xlabel(f"PC1 - {var[0] * 100:.2f}%")
    axes[i].set_ylabel(f"PC2 - {var[1] * 100:.2f}%")
    axes[i].legend(handles = handles)

#save the plots as PDF
plt.tight_layout()
plt.savefig(f"{DIR}/Visualizations/model_decision_boundaries.pdf")
plt.close()
