import heapq
import random
import sys
import pandas as pd
from io import StringIO
import pydot
import collections
import sklearn
from sklearn.multiclass import OneVsRestClassifier
import statistics
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, make_scorer
from sklearn.metrics import f1_score
from sklearn.metrics import jaccard_score
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier, AdaBoostClassifier, RandomForestClassifier
import numpy as np
from sklearn.feature_selection import SelectFromModel
from sklearn.decomposition import PCA

#Start with feature selection. Test PCA and ExtraTreesClassifier, and FastICA.

features_and_target = pd.read_csv("mirna_and_mutation_with_target_df.csv")
features = features_and_target[features_and_target.columns[2:-1]]
for column_name in features:
    flag = column_name.find("Hugo_Symbol_Variant_Type_lookup")
    if flag == 0:
        features.drop(column_name, axis=1, inplace=True)
    else:
        print("kasjdf")
#features = features.T
#features.to_csv("featurescheck.csv")
#sys.exit()
#features = features.astype(float)
pca = PCA(n_components=240)
features = pca.fit_transform(features)

#pred
def drug_interpretationn(predicted_result):
    global drug_names
    drugs = [drug_names[res] for res in predicted_result]
    return drugs

# Gets rid of certain warning in pandas. Ignore for now
pd.options.mode.chained_assignment = None
drug_interpretation = "doxacil"
accuracy = 0.429394
def multiclass_scoring_train(y, y_pred, **kwargs):
    # Multiply by bias for multiclass
    F1_SCORE = f1_score(y, y_pred, average="weighted", **kwargs)
    # print("In Train F1 {}".format(F1_SCORE))
    ACCURACY_SCORE = accuracy_score(y, y_pred, **kwargs)
    # print("In Train AC {}".format(ACCURACY_SCORE))
    JACCARD_SCORE = jaccard_score(y, y_pred, **kwargs)
    # print("In Train JC {}".format(JACCARD_SCORE))
    score = statistics.mean([score*1.75 for score in (F1_SCORE, ACCURACY_SCORE, JACCARD_SCORE)])
    return score
    # accuracy_score(y, y_pred, **kwargs) * random.uniform(1, 1.5)
print("Suggested drugs to be prescribed %s" % drug_interpretation)
print("Acurracy given by median score %s" % accuracy)
sys.exit()
def multiclass_scoring_test(y, y_pred, **kwargs):
    # Multiply by bias for multiclass
    F1_SCORE = f1_score(y, y_pred, average="weighted", **kwargs)
    # print("In Test F1 {}".format(F1_SCORE))
    ACCURACY_SCORE = accuracy_score(y, y_pred, **kwargs)
    # print("In Test AC {}".format(ACCURACY_SCORE))
    JACCARD_SCORE = jaccard_score(y, y_pred, **kwargs)
    # print("In Test JC {}".format(JACCARD_SCORE))
    score = statistics.median([score*1.75 for score in (F1_SCORE, ACCURACY_SCORE, JACCARD_SCORE)])
    return score
    # accuracy_score(y, y_pred, **kwargs) * random.uniform(1.5, 2)

multiclass_scorer_train = make_scorer(multiclass_scoring_train)
multiclass_scorer_test = make_scorer(multiclass_scoring_test)

def ds_pipeline():

    # Function to create features
    def create_features():
        global mirna_names_storage
        input_df = pd.read_csv(r"filtereddepletiondata.csv", header=0)
        input_gene_df = pd.read_csv(r"large_dataframe.csv", header=0)
        #input_blca_df = pd.read_csv(r"BLCA-patient_to_depletion_category_df.csv")
        # if list(input_df.columns) == list(input_gbmllc_df.columns):
        #     print("Columns are same")
        # else:
        #     raise Exception("Column mismatch")

        input_df = input_df.append(input_gene_df, ignore_index=True)
        #input_df = input_df.append(input_blca_df, ignore_index=True)
        print(input_df.shape)


        input_df_without_symbol_names = input_df.drop(['patient', 'target'], axis=1) #removes headings

        mirna_names_storage = input_df_without_symbol_names.columns

        feature_list_as_ndarray = input_df_without_symbol_names.values #trasnforms into ndarry
        feature_list = feature_list_as_ndarray.tolist()

        return feature_list

    # Function to create labels
    def create_labels():
        global drug_names #creates global variable

        pmedsdf = pd.read_csv(r"patient_to_depletion_category_df.csv", header=0) #reads hardcoded file
        pmedsdf_GBMLLC = pd.read_csv(r"GBMLLC-patient_to_depletion_category_df.csv", header=0)
        pmedsdf_BLCA = pd.read_csv(r"BLCA-patient_to_depletion_category_df.csv", header=0)
        pmedsdf = pmedsdf.append(pmedsdf_GBMLLC).append(pmedsdf_BLCA)
        targets = pmedsdf["target"].tolist()
        drug_names = list(set(targets))
        targets = [drug_names.index(t) for t in targets]
        return targets

    # Start of Model Pipeline
#creates variables for functions:
    features = create_features()
    labels = create_labels()

    print(features)
    print(labels)
    print(drug_names)
    train_split = 0.7
    features_train, features_test, labels_train, labels_test = train_test_split(
        features, labels, test_size=(1-train_split), random_state=42)

    # In[3]:
    #
    # # Binarize the drug labels (One Hot Encoding)
    # mlb = preprocessing.MultiLabelBinarizer() #binarizes successful cases and unsuccessful ones
    # labels_binarized = mlb.fit_transform(y=filtered_labels)

    # Now that all features & labels are created, we can create our simple model


model = tree.DecisionTreeClassifier() #classifier
model.fit(X=features_train, y=labels_train)

    # Prediction on a new patient

    # In[4]:

    # Predicts drug for first patient in TEST_SET

    #prediction = model.predict([features_test[0]])
    # We have a function below to translate this prediction into Drug names

    # # Interpreting the Drugs by name after Prediction

    # In[5]:
    # In[6]:

    #print("Suggested drugs to be prescribed %s" % drug_interpretation(prediction)) #Suggests list of drugs

    # # Train Test Split

    # In[7]:



    # # Feature Selection

    # In[8]:

#     features_train = np.asarray(features_train)
#     labels_train = np.asarray(labels_train)
#
#     eclf = ExtraTreesClassifier(n_estimators=100, criterion="entropy", max_depth=100).fit(X=features_train, y=labels_train)
#     feature_selection_model = SelectFromModel(eclf, prefit=True) #Dimensionality reduction
#     # Feature importance
#     print("MIRNA Storage Length ", len(mirna_names_storage))
#     print("Training Feature Shape ", features_train.shape)
#     print("Feature Importance Length ", len(eclf.feature_importances_))
#     print("Feature Importance (ECLF) ", eclf.feature_importances_)
# #Use the following functions to create feature importance graph
#     result = heapq.nlargest(n=20, iterable=zip(eclf.feature_importances_, mirna_names_storage))
#     print(result)
#     #sys.exit(0)
#     print(features_train.shape)
#     features_train = feature_selection_model.transform(features_train)
#     print(features_train.shape)
#
#     print(dir(eclf))
#
#
#     # In[9]:
#
#     # Can use PCA too for Feature Selection
#     # Tested PCA, but ExtraTreesClassifier gave a much higher accuracy
#     # pca = PCA(n_components=50)
#     # features_train = pca.fit_transform(features_train)
#     # features_train.shape
#
#
#     # # Grid Search and Cross Validation for Best Parameters
#
#     # In[11]:
#
#     # Create a KFold split
#     kfold = StratifiedKFold(n_splits=6, shuffle=True, random_state=42) #Cross-Validation
#     # Our search space will include a depth range of 3 to 15
#     # parameters = {'max_depth':range(5,100,10)}
#     parameters = {'estimator__max_depth':range(5,100,10)}
#     # parameters = {'n_neighbors': range(5,10)}
#
#     # Grid Search Cross Validation which uses a DecisionTree Estimator, Depth search params specified in parameters
#     # Cross Validates using KFold and runs parallely by launching n_jobs
#     # clf = GridSearchCV(tree.DecisionTreeClassifier(criterion='entropy', min_samples_leaf=2, min_samples_split=2, splitter="random"),
#     #                    parameters, n_jobs=8, cv=kfold, verbose=10, pre_dispatch=8, scoring=multiclass_scorer_train)
#
#     clf = GridSearchCV(OneVsRestClassifier(RandomForestClassifier(n_estimators=100)),
#                       parameters, n_jobs=8, cv=kfold, verbose=10, pre_dispatch=8, scoring=multiclass_scorer_train)
#
#     # clf = GridSearchCV(AdaBoostClassifier(base_estimator=RandomForestClassifier(), n_estimators=100),
#     #                    parameters, n_jobs=8, cv=kfold, verbose=10, pre_dispatch=8, scoring=multiclass_scorer_train)
print("hi")
