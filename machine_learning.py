# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 19:00:38 2017

@author: Karim
"""

import scipy
import numpy
import matplotlib
import pandas as pd
import sklearn
from sklearn import model_selection
print('sklearn: {}'.format(sklearn.__version__))

#Dataset without X values (without resampling):
dataset = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_machine_learning/ml_data_v2.csv', sep = ',')
dataset = dataset.drop(['Unnamed: 0'], axis=1)

buffer = []
print(dataset.shape)

# descriptions
print(dataset.describe())

# class distribution
print(dataset.groupby('motifs').size())
# There is a bias in the data. Each class should have the same amount of instance!
# Try with 1000 of each motifs levels

###############################################################################
# Validation set
###############################################################################
# Prediction of DAIC knowing motifs and conservation
conserv = dataset['conserv']
motifs = dataset['motifs']
one_hot = pd.get_dummies(motifs, 'motifs')
X = pd.concat([one_hot, conserv], axis=1)
Y = dataset['DAIC']
validation_size = 0.10 # 10% for the validation set
seed = 7


# (Variable motifs as Y)
one_hot = pd.get_dummies(motifs, 'motifs')
X = dataset[['DAIC', 'conserv']]
Y = one_hot
validation_size = 0.10
seed = 7

X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)

###############################################################################
# Model
###############################################################################
# Load libraries
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor

###############################################################################    
# Prediction of DAIC knowing secondary structure patterns and conservation
###############################################################################

#Decision Tree Regressor
cart = DecisionTreeRegressor()
cart.fit(X_train, Y_train)
predictions = cart.predict(X_validation)

results = model_selection.cross_val_score(cart, X, Y, cv=10) #cv = 10-fold ?????????? Demander qu'est-ce que X et Y.
print("Accuracy: %.3f%% (%.3f%%)" % (results.mean()*100.0, results.std()*100.0))

pred_vs_valid = pd.DataFrame(Y_validation)
pred_vs_valid['predictions'] = predictions
pred_vs_valid.to_csv('pred_vs_valid.csv', sep = ',')

###############################################################################    
# Predictions of secondary_structures knowing the DAIC and conservation.
###############################################################################

#Decision tree (best model compared to the others bellow)
cart = DecisionTreeClassifier()
cart.fit(X_train, Y_train)
predictions = cart.predict(X_validation)
print(accuracy_score(Y_validation, predictions)) #For classification task only
print(confusion_matrix(Y_validation, predictions))
print(classification_report(Y_validation, predictions))

clf = DecisionTreeRegressor()
clf = clf.fit(X_train, Y_train)
tree.export_graphviz(clf, out_file='tree_clf.dot')
