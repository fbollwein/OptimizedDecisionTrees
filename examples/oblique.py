#include the lib folder where the DecisionTree module is stored
import sys
sys.path.append('../lib')
#import the DecisionTree module
import DecisionTree
#to load the dataset
from sklearn import datasets
#import for visualization
import graphviz

#Load iris dataset
iris = datasets.load_iris()
X = iris.data
y = iris.target

#Split into training and test data
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)

#Initialize decision tree
bivariateTree=DecisionTree.DecisionTree()
#Enable bivariate oblique splits
bivariateTree.setParam("2d",1)
#Specify the names of the classes
bivariateTree.setYNames([0,1,2],iris.target_names)
#Specify the names of the attributes
bivariateTree.setAttributeNames(iris.feature_names)
#Fit the decision tree
bivariateTree.fit(X_train,y_train)
#Visualize decision tree
graph = graphviz.Source(bivariateTree.export_graphviz(False))
graph.render("bivariate",view=True,cleanup=True)
#Predict labels of test set
y_pred=bivariateTree.predict(X_test)
#Calculate accuracy
accuracy=1-len([i for i in range(len(y_test)) if y_test[i]!=y_pred[i]])/len(y_test)
print("Bivariate oblique decision tree accuracy: ",accuracy*100,"%",sep='')

#Fit oblique decision tree with both the cross-entropy and the simulated annealing algorithm enabled
obliqueTree=DecisionTree.DecisionTree()
#Enable the cross-entropy method
obliqueTree.setParam("ce",1)
#Enable the simulating annealing method
obliqueTree.setParam("sa",1)
#Enables post-processing of oblique splits by solving a generalized svm problem
obliqueTree.setParam("max_margin",1)
#Specify the names of the classes
obliqueTree.setYNames([0,1,2],iris.target_names)
#Specify the names of the attributes
obliqueTree.setAttributeNames(iris.feature_names)
#Fit the decision tree
obliqueTree.fit(X_train,y_train)
#Visualize decision tree
graph = graphviz.Source(obliqueTree.export_graphviz(False))
graph.render("oblique",view=True,cleanup=True)
#Predict labels of test set
y_pred=obliqueTree.predict(X_test)
#Calculate accuracy
accuracy=1-len([i for i in range(len(y_test)) if y_test[i]!=y_pred[i]])/len(y_test)
print("Oblique decision tree accuracy: ",accuracy*100,"%",sep='')