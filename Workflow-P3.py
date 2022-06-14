import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import tensorflow as tf
import math
from numpy import asarray
from numpy import save
from numpy import load
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, MaxPool2D, Dropout, Flatten, Input
from sklearn.preprocessing import StandardScaler
import os
import tempfile
from tensorflow import keras
from tensorflow.keras.optimizers import RMSprop
import matplotlib as mpl
import sklearn
from collections import Counter
from google.colab import drive

drive.mount('/content/drive')

dict_train = load('/content/drive/MyDrive/M2-int/TestCCN/frames_train.npz')

X_train = dict_train['arr_0']
print(len(X_train))
y_flt_train = load('/content/drive/MyDrive/M2-int/TestCCN/labels_train.npy')
y_train = y_flt_train.astype(int)
#test
dict_test = load('/content/drive/MyDrive/M2-int/TestCCN/frames_test.npz')
X_test = dict_test['arr_0']


y_flt_test = load('/content/drive/MyDrive/M2-int/TestCCN/labels_test.npy')
y_test = y_flt_test.astype(int)

#Group specific type of NtCs

listeAA = [13,14,18,23,24,26,30,38,40,42,46,53,62,71]
#make a binary problem 
for e in range(len(y_train)):
  if y_train[e] in listeAA:
    y_train[e] = 1
  else:
    y_train[e] = 0

for e in range(len(y_test)):
  if y_test[e] in listeAA:
    y_test[e] = 1
  else:
    y_test[e] = 0 

set(y_train)
set(y_test)

print(y_train)

print(len(y_train))

"""
dict_train = load('/content/drive/MyDrive/M2-int/TestCCN/frames_train.npz')

X_train = dict_train['arr_0']
print(len(X_train))
y_flt_train = load('/content/drive/MyDrive/M2-int/TestCCN/labels_train.npy')
y_train = y_flt_train.astype(int)
#test
dict_test = load('/content/drive/MyDrive/M2-int/TestCCN/frames_test.npz')
X_test = dict_test['arr_0']


y_flt_test = load('/content/drive/MyDrive/M2-int/TestCCN/labels_test.npy')
y_test = y_flt_test.astype(int)

#Group specific type of NtCs

listeOP = [5,6,7,8,9,10,11,12,17,19,21,25,35,37,39,41,44,55,56,57,58,60,61,63,64,65,72,73]
#make a binary problem 
for e in range(len(y_train)):
  if y_train[e] in listeOP:
    y_train[e] = 1
  else:
    y_train[e] = 0

for e in range(len(y_test)):
  if y_test[e] in listeOP:
    y_test[e] = 1
  else:
    y_test[e] = 0 
    
set(y_train)
set(y_test)

print(y_train)

print(len(y_train))
"""

#positive represent the desired NtCs
neg, pos = np.bincount(y_train)+np.bincount(y_test)
total = neg + pos
print('Examples:\n    Total: {}\n    Positive: {} ({:.2f}% of total)\n'.format(
    total, pos, 100 * pos / total))

print('Training labels shape:', y_train.shape)
print('Test labels shape:', y_test.shape)

print('Training features shape:', X_train.shape)
print('Test features shape:', X_test.shape)


METRICS = [
      keras.metrics.TruePositives(name='tp'),
      keras.metrics.FalsePositives(name='fp'),
      keras.metrics.TrueNegatives(name='tn'),
      keras.metrics.FalseNegatives(name='fn'),
      keras.metrics.BinaryAccuracy(name='accuracy'),
      keras.metrics.Precision(name='precision'),
      keras.metrics.Recall(name='recall'),
      keras.metrics.AUC(name='auc'),
      keras.metrics.AUC(name='prc', curve='PR'), # precision-recall curve
]

def make_model(metrics=METRICS, output_bias=None):
  if output_bias is not None:
    output_bias = tf.keras.initializers.Constant(output_bias)
  model = tf.keras.models.Sequential([
    # Note the input shape is the desired size of the image 11x11 with 25 bytes color
    # This is the first convolution
    tf.keras.layers.Conv2D(16, (3,3), activation='relu', input_shape=(11, 11, 25)),
    tf.keras.layers.MaxPooling2D(2, 2),
    tf.keras.layers.Dropout(0.25),
    # The second convolution
    tf.keras.layers.Conv2D(32, (3,3), activation='relu'),
    tf.keras.layers.MaxPooling2D(2,2),
    tf.keras.layers.Dropout(0.25),
    # Flatten the results to feed into a DNN
    tf.keras.layers.Flatten(),
    # 512 neuron hidden layer
    tf.keras.layers.Dense(512, activation='relu', bias_initializer=output_bias),
    tf.keras.layers.Dropout(0.5),
    # Only 1 output neuron. It will contain a value from 0-1 where 0 for class ('not contact') and 1 for the other ('contact')
    tf.keras.layers.Dense(1, activation='sigmoid')
  ])

  model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-3), #optimizer=RMSprop(lr=0.001),
    loss=keras.losses.BinaryCrossentropy(),
    metrics=metrics)

  return model

EPOCHS = 30
BATCH_SIZE = 2048 #large batch size to ensure each batch contains a few contacts




mpl.rcParams['figure.figsize'] = (15, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def plot_loss(history, label, n):
  # Use a log scale on y-axis to show the wide range of values.
  plt.semilogy(history.epoch, history.history['loss'],
               color=colors[n], label='Train ' + label)
  plt.semilogy(history.epoch, history.history['loss'],
               color=colors[n], label='Test ' + label,
               linestyle="--")
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  
  
def plot_cm(labels, predictions, p=0.5):
  cm = confusion_matrix(labels, predictions > p)
  plt.figure(figsize=(5,5))
  sns.heatmap(cm, annot=True, fmt="d")
  plt.title('Confusion matrix')
  plt.ylabel('Actual label')
  plt.xlabel('Predicted label')

  print('Not NtCs Detected (True Negatives): ', cm[0][0])
  print('Not NtCs Incorrectly Detected (False Positives): ', cm[0][1])
  print('NtCs Missed (False Negatives): ', cm[1][0])
  print('NtCs Truely Detected (True Positives): ', cm[1][1])
  print('Total Contacts: ', np.sum(cm[1]))

  def plot_metrics(history):
  metrics = ['loss', 'prc', 'precision', 'recall']
  for n, metric in enumerate(metrics):
    name = metric.replace("_"," ").capitalize()
    plt.subplot(2,2,n+1)
    plt.plot(history.epoch, history.history[metric], color=colors[0], label='Train')
    #plt.plot(history.epoch, history.history[metric],color=colors[1], linestyle="--", label='Test')
    plt.xlabel('Epoch')
    plt.ylabel(name)
    if metric == 'loss':
      plt.ylim([0, plt.ylim()[1]])
    elif metric == 'auc':
      plt.ylim([0.8,1])
    else:
      plt.ylim([0,1])

    plt.legend()

def plot_roc(name, labels, predictions, **kwargs):
  fp, tp, _ = sklearn.metrics.roc_curve(labels, predictions)

  plt.plot(100*fp, 100*tp, label=name, linewidth=2, **kwargs)
  plt.xlabel('False positives [%]')
  plt.ylabel('True positives [%]')
  plt.xlim([-5,100])
  plt.ylim([0,100])
  plt.grid(True)
  ax = plt.gca()
  ax.set_aspect('equal')
  
def plot_prc(name, labels, predictions, **kwargs):
    precision, recall, _ = sklearn.metrics.precision_recall_curve(labels, predictions)

    plt.plot(precision, recall, label=name, linewidth=2, **kwargs)
    plt.xlabel('Recall  - TP / (TP + FP)')
    plt.ylabel('Precision - TP / (TP + FN)')
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    
#retrain with class weights
# Scaling by total/2 helps keep the loss to a similar magnitude.
# The sum of the weights of all examples stays the same.
weight_for_0 = (1 / neg) * (total / 2.0)
weight_for_1 = (1 / pos) * (total / 2.0)

class_weight = {0: weight_for_0, 1: weight_for_1}

print('Weight for class 0: {:.2f}'.format(weight_for_0))
print('Weight for class 1: {:.2f}'.format(weight_for_1))


weighted_model = make_model()
#weighted_model.load_weights(initial_weights)

weighted_history = weighted_model.fit(
    X_train,
    y_train,
    batch_size=BATCH_SIZE,
    epochs=EPOCHS,

    # The class weights go here
    class_weight={0: 0.50, 1: 144.48})
weighted_results = weighted_model.evaluate(X_test, y_test,
                                           batch_size=BATCH_SIZE, verbose=0)
for name, value in zip(weighted_model.metrics_names, weighted_results):
  print(name, ': ', value)
print()

plot_cm(y_test, test_predictions_weighted)


plot_roc("Train Weighted", y_train, train_predictions_weighted, color=colors[4])
plot_roc("Test Weighted", y_test, test_predictions_weighted, color=colors[4], linestyle='--')
plt.legend(loc='lower right')
#plt.savefig('ROC.png')

plot_prc("Train Weighted", y_train, train_predictions_weighted, color=colors[4])
plot_prc("Test Weighted", y_test, test_predictions_weighted, color=colors[4], linestyle='--')


plt.legend(loc='lower right')
#plt.savefig('PRC.png')
