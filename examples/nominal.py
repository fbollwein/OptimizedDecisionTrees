#include the lib folder where the DecisionTree module is stored
import sys
sys.path.append('../lib')
#import the DecisionTree module
import DecisionTree
#import for visualization
import graphviz
import numpy as np
import pandas as pd

#Load the credit approval dataset from the UCI Machine Learning Repository
df = pd.read_csv('https://archive.ics.uci.edu/ml/machine-learning-databases/credit-screening/crx.data',header=None,na_values="?",keep_default_na=False).dropna()

#Encode nominal features and the target variable
nominal_encodings={}
nominal=[0,3,4,5,6,8,9,11,12,15]
for i in nominal:
    nominal_encodings[i]={}
    enc=0
    for v in set(df[i]):
        nominal_encodings[i][v]=enc
        enc+=1
df=df.replace(nominal_encodings)

#Extract dataset
y=df[15].to_numpy()
X=df[[i for i in range(15)]].to_numpy()

#Initialize decision tree
nominalTree=DecisionTree.DecisionTree()

#Specify attribute names
attributeNames=["A"+str(i) for i in range(len(X[0]))]
nominalTree.setAttributeNames(attributeNames)

#Specify names for the target labels
categories=[c for c in nominal_encodings[15]]
values=[nominal_encodings[15][c] for c in categories]
nominalTree.setYNames(values,categories)

#Specify names for the categories of the nominal attributes
nominal_attributes=[i for i in nominal if i!=15]
for i in nominal_attributes:
    categories=[c for c in nominal_encodings[i]]
    values=[nominal_encodings[i][c] for c in categories]
    nominalTree.setNominalNames(i,values,categories)

#Specify maximal depth for better visualization
nominalTree.setParam("max_depth",4)
#Fit the decision tree
nominalTree.fit(X,y,nominal_attributes)
#Visualize the decision tree
graph = graphviz.Source(nominalTree.export_graphviz(False))
graph.render("nominal",view=True,cleanup=True)


