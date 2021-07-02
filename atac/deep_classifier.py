import tensorflow as tf
from tensorflow.keras import optimizers, regularizers
from tensorflow.keras.models import Model, load_model, Sequential
from tensorflow.keras.layers import Dense, Activation, BatchNormalization, Dropout
import numpy as np


class CLASSIFIER:
    def __init__(self, input_size, class_num=2, path=''):
        self.input_size = input_size
        self.classifier = None
        self.initializers = "glorot_uniform"
        self.optimizer = optimizers.Adam(lr=0.01)
        self.validation_split = 0.1
        self.class_num = class_num
        self.dropout_rate = 0.05

    def build(self):
        model = Sequential()
        model.add(Dense(64, activation='relu', input_shape=(self.input_size,)))
        model.add(Dropout(rate=self.dropout_rate))
        model.add(Dense(32, activation='relu'))
        model.add(Dropout(rate=self.dropout_rate))
        model.add(Dense(self.class_num, activation='softmax'))
        self.classifier = model

    def compile(self):
        self.classifier.compile(optimizer=self.optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
        self.classifier.summary()

    def train(self, x, label, batch_size=100, epochs=300):
        history = self.classifier.fit(x, label,
                                    epochs=epochs, batch_size=batch_size,
                                    validation_split=self.validation_split, shuffle=True)
        return history

    def prediction(self, x):
        label = self.classifier.predict(x)
        label = np.argmax(label, axis=1)
        return label


# ----------test----------------
from keras.utils import to_categorical
import scanpy as sc
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from collections import Counter

base_path = '/Users/zhongyuanke/data/'

file1 = base_path+'seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix_8000.csv'
file2 = base_path+'seurat_data/sc_atac/pbmc_10k_v3_8000.csv'
davae_path = base_path+'dann_vae/atac/davae_save03.h5ad'

celltype_x_path = base_path + 'seurat_data/sc_atac/pbmc_10k_v3_celltype.csv'
celltype_x = pd.read_csv(celltype_x_path, index_col=0)
celltype_x = celltype_x.values
print(celltype_x)

seurat_celltype_path = base_path + 'dann_vae/atac/seurat_pred_type.csv'
celltype_seurat = pd.read_csv(seurat_celltype_path, index_col=0)
celltype_seurat = list(celltype_seurat.values.flatten())
print(celltype_seurat)

encoder = LabelEncoder()
orig_label = encoder.fit_transform(celltype_x)
orig_label.dtype='int64'

batch_size = 100
epochs = 25

adata1 = sc.read_csv(file1)
adata2 = sc.read_csv(file2)
adata_davae = sc.read_h5ad(davae_path)
data = adata_davae.X
# data = adata_davae.obsm['davae']

len1 = adata1.shape[0]
len2 = adata2.shape[0]

test_set = data[0:len1, ]
train_set = data[len1:len1+len2, ]

label = to_categorical(orig_label)
class_num = label.shape[1]

net_x = CLASSIFIER(input_size=train_set.shape[1], class_num=class_num)
net_x.build()
net_x.compile()
his = net_x.train(x=train_set, label=label, epochs=epochs, batch_size=batch_size)
pred_label = net_x.prediction(test_set)
pred_label.dtype='int64'
pred_type = encoder.inverse_transform(pred_label)

df = pd.DataFrame(pred_type)
df.to_csv('/Users/zhongyuanke/data/dann_vae/atac/pred_type_save03.csv')
np.savetxt('/Users/zhongyuanke/data/dann_vae/atac/pred_label_save03.csv', pred_label, delimiter=',')

all_label = np.concatenate([pred_label, orig_label])
all_type = encoder.inverse_transform(all_label)

print(pred_type)
type_list = list(pred_type)

print(Counter(type_list))
print(Counter(celltype_seurat))

adata_davae.obs['label'] = all_label
adata_davae.obs['cell type'] = all_type
# adata_davae.write_h5ad(base_path+'dann_vae/atac/davae_save03_para.h5ad')
