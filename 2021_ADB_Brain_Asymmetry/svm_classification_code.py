import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn import svm
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from scipy import interp
import matplotlib.pyplot as plt
import os

# Define classifier
def run_svm(X,y,savename):    

    rs = np.random.RandomState(111)
    cv = StratifiedKFold(n_splits=10,random_state=rs)

    Cs = np.logspace(-3,3,15)
    gammas = np.logspace(-3,3,15)
    param_grid = {'C': Cs, 'gamma' : gammas}
    clf = GridSearchCV(svm.SVC(kernel='rbf',probability=True),  param_grid, scoring='roc_auc', cv=cv,n_jobs=4)

    ###############################################################################
    # Classification and ROC analysis
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []

    plt.figure(figsize=(15,15))
    for i, (train, test) in enumerate(cv.split(X,y)): # Nested CV
        probas_ = clf.fit(X[train], y[train]).predict_proba(X[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')

    mean_tpr /= i+1
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, 'k--',
            label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic Dependence prediction')
    plt.legend(loc="lower right")
    plt.savefig(savename)
    #plt.show()


os.chdir(r'F:\Google Drive\post-doc\Laterality\Manuscript\Brain_Asymmetry_upload_ADB\code2upload')
# For Model 1
df=pd.read_csv('residulized_aidata.csv')
X=np.array(df[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df['Dependentanydrug'])
run_svm(X,y,'ML_model_case_control.tiff')
del X, y

# For Model 2
df_alc=df.loc[(df['PrimaryDrug']==0)|(df['PrimaryDrug']==1)]
X=np.array(df_alc[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df_alc['PrimaryDrug'])
run_svm(X,y,'ML_model_alc_control.tiff')
del X, y

df_nic=df.loc[(df['PrimaryDrug']==0)|(df['PrimaryDrug']==2)]
X=np.array(df_nic[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df_nic['PrimaryDrug'])
y[y==2]=1
run_svm(X,y,'ML_model_nic_control.tiff')

df_coc=df.loc[(df['PrimaryDrug']==0)|(df['PrimaryDrug']==3)]
X=np.array(df_coc[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df_coc['PrimaryDrug'])
y[y==3]=1
run_svm(X,y,'ML_model_coc_control.tiff')

df_met=df.loc[(df['PrimaryDrug']==0)|(df['PrimaryDrug']==4)]
X=np.array(df_met[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df_met['PrimaryDrug'])
y[y==4]=1
run_svm(X,y,'ML_model_met_control.tiff')

df_can=df.loc[(df['PrimaryDrug']==0)|(df['PrimaryDrug']==5)]
X=np.array(df_can[[n for n in df.columns if n.startswith('AI')]])
y=np.array(df_can['PrimaryDrug'])
y[y==5]=1
run_svm(X,y,'ML_model_can_control.tiff')


