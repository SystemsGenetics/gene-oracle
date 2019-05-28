import numpy as np
import tensorflow as tf
import sys, argparse
import os
import matplotlib.pyplot as plt



from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn import preprocessing


def confusion_heatmap(conf_arr, labels=None):
    norm_conf = preprocessing.normalize(conf_arr, axis=1, norm='l1')

    fig = plt.figure(figsize=(14,8))
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)

    res = ax.imshow(norm_conf, cmap=plt.cm.jet,
            interpolation='nearest')

    width, height = conf_arr.shape

    # for x in range(width):
    #     for y in range(height):
    #         ax.annotate(str(conf_arr[x][y]), xy=(y, x),
    #                     horizontalalignment='center',
    #                     verticalalignment='center')

    cb = fig.colorbar(res)
    if labels == None:
        plt.xticks(range(width), np.arange(0,conf_arr.shape[0]), rotation='vertical')
        plt.yticks(range(height), np.arange(0,conf_arr.shape[0]))
    else:
        plt.xticks(range(width), labels, rotation='vertical')
        plt.yticks(range(height), labels)
    #plt.show()
    plt.savefig('confusion_matrix.png', format='png')

#Tutorial Here: http://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html#sphx-glr-auto-examples-model-selection-plot-roc-py
#y_score = pred
#y_test = actual
def roc_plt(n_classes,y_test,y_score,labels):
    # Compute ROC curve and ROC area for each class
    fig = plt.figure(figsize=(24,16))
    fpr = dict()
    tpr = dict()
    lw =2
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    #colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
    #for i, color in zip(range(n_classes), colors):
    for i in range(n_classes):
        #plt.legend(fontsize=50) # using a size in points
        #plt.legend(fontsize="x-large") # using a named size
        plt.plot(fpr[i], tpr[i], lw=lw,
                 label='{0} (area = {1:0.2f})'
                 ''.format(labels[i], roc_auc[i]))
        #plt.legend(loc=2, prop={'size': 20})

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('panTCGA ROC Curves')
    plt.legend(loc="lower right")
    plt.savefig('roc.png', format='png')
