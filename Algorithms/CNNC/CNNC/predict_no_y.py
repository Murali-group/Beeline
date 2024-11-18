# Usage: python predict_no_y.py  number_of_separation NEPDF_pathway model_pathway
# command line in developer's linux machine :
# python predict_no_y.py  9 /home/yey3/cnn_project/code3/NEPDF_data   /home/yey3/cnn_project/code3/trained_model/models/KEGG_keras_cnn_trained_model_shallow2.h5
from __future__ import print_function
import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping,ModelCheckpoint
import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy import interp
################################
#num_classes = 3 ################################# the number of classes might vary for your special application
#length_TF =2 # number of data parts divided
def load_data_TF2(indel_list,data_path): # cell type specific  ## random samples for reactome is not enough, need borrow some from keggp
    import random
    import numpy as np
    xxdata_list = []
    yydata = []
    count_set = [0]
    count_setx = 0
    for i in indel_list:#len(h_tf_sc)):
        xdata = np.load(data_path+'/Nxdata_tf' + str(i) + '.npy')
        for k in range(xdata.shape[0]):
            xxdata_list.append(xdata[k,:,:,:])
        count_setx = count_setx + xdata.shape[0]
        count_set.append(count_setx)
    return((np.array(xxdata_list),count_set))


length_TF =int(sys.argv[1]) # number of data parts divided
data_path = sys.argv[2]
num_classes = int(sys.argv[3])
model_path = sys.argv[4] ## KEGG or Reactome or TF
print ('select', type)
whole_data_TF = [i for i in range(length_TF)]
test_TF = [i for i in range (length_TF)]
(x_test, count_set) = load_data_TF2(test_TF,data_path)
print(x_test.shape, 'x_test samples')
############

#########################  model structure
model = Sequential()
model.add(Conv2D(32, (3, 3), padding='same',input_shape=x_test.shape[1:]))
model.add(Activation('relu'))
model.add(Conv2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(64, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(64, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(128, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(128, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(512))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes))
model.add(Activation('softmax'))
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(optimizer=sgd,loss='categorical_crossentropy',metrics=['accuracy'])
###########################################################3
save_dir = os.path.join(os.getcwd(),'predict_results_no_y_1')
model.load_weights(model_path)
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
print ('load model and predict')
y_predict = model.predict(x_test)
#np.save(save_dir+'/y_test.npy',y_test)
np.save(save_dir+'/y_predict.npy',y_predict)
s = open (save_dir+'/gene_index.txt','w')
for i in count_set:
    s.write(str(i)+'\n')
s.close()
######################################