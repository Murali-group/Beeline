from __future__ import print_function
# Usage  python train_with_labels_three_fold.py number_of_data_parts_divided NEPDF_pathway number_of_category
# command line in developer's linux machine :
# module load cuda-8.0 using GPU
#srun -p gpu --gres=gpu:1 -c 2 --mem=20Gb python train_with_labels_three_foldx.py 9 /home/yey3/cnn_project/code3/NEPDF_data 3 > results.txt
#######################OUTPUT
# it will generate three-fold cross validation results
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
####################################### parameter settings
data_augmentation = False
# num_predictions = 20
batch_size = 1024 # mini batch for training
#num_classes = 3   #### categories of labels
epochs = 20     #### iterations of trainning, with GPU 1080, 600 for KEGG and Reactome, 200 for tasks for GTRD
#length_TF =3057  # number of divide data parts
# num_predictions = 20
model_name = 'keras_cnn_trained_model_shallow.h5'
###################################################


def load_data_TF2(indel_list,data_path): # cell type specific  ## random samples for reactome is not enough, need borrow some from keggp
    import random
    import numpy as np
    xxdata_list = []
    yydata = []
    count_set = [0]
    count_setx = 0
    for i in indel_list:#len(h_tf_sc)):
        xdata = np.load(data_path+'/Nxdata_tf' + str(i) + '.npy')
        ydata = np.load(data_path+'/ydata_tf' + str(i) + '.npy')
        for k in range(len(ydata)):
            xxdata_list.append(xdata[k,:,:,:])
            yydata.append(ydata[k])
        count_setx = count_setx + len(ydata)
        count_set.append(count_setx)
        print (i,len(ydata))
    yydata_array = np.array(yydata)
    yydata_x = yydata_array.astype('int')
    print (np.array(xxdata_list).shape)
    return((np.array(xxdata_list),yydata_x,count_set))

if len(sys.argv) < 4:
    print ('No enough input files')
    sys.exit()
length_TF =int(sys.argv[1]) # number of data parts divided
data_path = sys.argv[2]
num_classes = int(sys.argv[3])
whole_data_TF = [i for i in range(length_TF)]
###################################################################################################################################
for test_indel in range(1,4): ################## three fold cross validation                                                     ## for  3 fold CV              
    test_TF = [i for i in range (int(np.ceil((test_indel-1)*0.333*length_TF)),int(np.ceil(test_indel*0.333*length_TF)))]         #
    train_TF = [i for i in whole_data_TF if i not in test_TF]                                                                    #
###################################################################################################################################
#####################################################################
    (x_train, y_train,count_set_train) = load_data_TF2(train_TF,data_path)
    (x_test, y_test,count_set) = load_data_TF2(test_TF,data_path)
    print(x_train.shape, 'x_train samples')
    print(x_test.shape, 'x_test samples')
    save_dir = os.path.join(os.getcwd(),str(test_indel)+'YYYY_saved_models_T_32-32-64-64-128-128-512_e'+str(epochs)) ## the result folder 
    if num_classes >2:
        y_train = keras.utils.to_categorical(y_train, num_classes)
        y_test = keras.utils.to_categorical(y_test, num_classes)
    print(y_train.shape, 'y_train samples')
    print(y_test.shape, 'y_test samples')
    ############
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    ############
    model = Sequential()
    model.add(Conv2D(32, (3, 3), padding='same',
                     input_shape=x_train.shape[1:]))
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
    if num_classes <2:
        print ('no enough categories')
        sys.exit()
    elif num_classes ==2:
        model.add(Dense(1, activation='sigmoid'))
        sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        model.compile(optimizer=sgd,loss='binary_crossentropy',metrics=['accuracy'])
    else:
        model.add(Dense(num_classes))
        model.add(Activation('softmax'))
        sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        model.compile(optimizer=sgd,loss='categorical_crossentropy',metrics=['accuracy'])

    early_stopping = keras.callbacks.EarlyStopping(monitor='val_acc', patience=50, verbose=0, mode='auto')
    checkpoint1 = ModelCheckpoint(filepath=save_dir + '/weights.{epoch:02d}-{val_loss:.2f}.hdf5', monitor='val_loss',
                                  verbose=1, save_best_only=False, save_weights_only=False, mode='auto', period=1)
    checkpoint2 = ModelCheckpoint(filepath=save_dir + '/weights.hdf5', monitor='val_acc', verbose=1,
                                  save_best_only=True, mode='auto', period=1)
    callbacks_list = [checkpoint2, early_stopping]
    if not data_augmentation:
        print('Not using data augmentation.')
        history = model.fit(x_train, y_train,
                  batch_size=batch_size,
                  epochs=epochs,validation_split=0.2,
                  shuffle=True, callbacks=callbacks_list)

    # Save model and weights

    model_path = os.path.join(save_dir, model_name)
    model.save(model_path)
    print('Saved trained model at %s ' % model_path)
    # Score trained model.
    scores = model.evaluate(x_test, y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])
    y_predict = model.predict(x_test)
    np.save(save_dir+'/end_y_test.npy',y_test)
    np.save(save_dir+'/end_y_predict.npy',y_predict)
############################################################################## plot training process
    plt.figure(figsize=(10, 6))
    plt.subplot(1,2,1)
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.grid()
    plt.legend(['train', 'val'], loc='upper left')
    plt.subplot(1,2,2)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.grid()
    plt.savefig(save_dir+'/end_result.pdf')
    ###############################################################  evaluation without consideration of data separation
    if num_classes == 3:  ## here we only focus on three category tasks
        plt.figure(figsize=(10, 6))
        for i in range(3):
            y_test_x = [j[i] for j in y_test]
            y_predict_x = [j[i] for j in y_predict]
            fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
            plt.subplot(1, 3, i + 1)
            plt.plot(fpr, tpr)
            plt.grid()
            plt.plot([0, 1], [0, 1])
            plt.xlabel('FP')
            plt.ylabel('TP')
            plt.ylim([0, 1])
            plt.xlim([0, 1])
            auc = np.trapz(tpr, fpr)
            print('AUC:', auc)
            plt.title('label' + str(i) + ', AUC:' + str(auc))
        plt.savefig(save_dir + '/end_3labels.pdf')
        plt.figure(figsize=(10, 6))
        y_predict_x = [j[1] + j[2] for j in y_predict]
        y_test_x = [1 - j[0] for j in y_test]
        fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
        # Print ROC curve
        plt.plot(fpr, tpr)
        plt.plot([0, 1], [0, 1])
        # Print AUC
        auc = np.trapz(tpr, fpr)
        print('AUC:', auc)
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        plt.grid()
        plt.title('label 1+2, AUC:' + str(auc))
        plt.xlabel('FP')
        plt.ylabel('TP')
        plt.savefig(save_dir + '/end_1+2.pdf')
        #######################################
        plt.figure(figsize=(10, 6))
        y_predict1 = []
        y_test1 = []
        x = 2
        for i in range(int(len(y_predict) / 3)):
            y_predict1.append(y_predict[3 * i][x] - y_predict[3 * i + 1][x]) #### here we prepared the data as (GeneA,GeneB),(GeneB,GeneA) and (GeneA,GeneX) as label 1, 2, 0, That is why we can predict direaction using this code
            y_predict1.append(-y_predict[3 * i][x] + y_predict[3 * i + 1][x])
            y_test1.append(y_test[3 * i][x])
            y_test1.append(y_test[3 * i + 1][x])
        fpr, tpr, thresholds = metrics.roc_curve(y_test1, y_predict1, pos_label=1)
        # Print ROC curve
        plt.plot(fpr, tpr)
        plt.plot([0, 1], [0, 1])
        # Print AUC
        auc = np.trapz(tpr, fpr)
        print('AUC:', auc)
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        plt.grid()
        plt.title('label 1 vs 2,direction diff, AUC:' + str(auc))
        plt.xlabel('FP')
        plt.ylabel('TP')
        plt.savefig(save_dir + '/end_1vs2.pdf')
        #############################################################
        ################################################ evaluation with data separation

        fig = plt.figure(figsize=(5, 5))
        plt.plot([0, 1], [0, 1])
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        plt.xlabel('FP')
        plt.ylabel('TP')
        # plt.grid()
        AUC_set = []
        y_testy = y_test
        y_predicty = y_predict
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        s = open(save_dir + '/divided_AUCs1vs2.txt', 'w')
        for jj in range(len(count_set) - 1):  # len(count_set)-1):
            if count_set[jj] < count_set[jj + 1]:
                print(jj, count_set[jj], count_set[jj + 1])
                y_test = y_testy[count_set[jj]:count_set[jj + 1]]
                y_predict = y_predicty[count_set[jj]:count_set[jj + 1]]
                y_predict1 = []
                y_test1 = []
                x = 2
                for i in range(int(len(y_predict) / 3)):
                    y_predict1.append(y_predict[3 * i][x] - y_predict[3 * i + 1][x])
                    y_predict1.append(-y_predict[3 * i][x] + y_predict[3 * i + 1][x])
                    y_test1.append(y_test[3 * i][x])
                    y_test1.append(y_test[3 * i + 1][x])
                fpr, tpr, thresholds = metrics.roc_curve(y_test1, y_predict1, pos_label=1)
                tprs.append(interp(mean_fpr, fpr, tpr))
                tprs[-1][0] = 0.0
                plt.plot(fpr, tpr, color='0.5', lw=0.1)
                auc = np.trapz(tpr, fpr)
                s.write(str(jj) + '\t' + str(count_set[jj]) + '\t' + str(count_set[jj + 1]) + '\t' + str(auc) + '\n')
                print('AUC:', auc)
                AUC_set.append(auc)
        mean_tpr = np.median(tprs, axis=0)
        mean_tpr[-1] = 1.0
        per_tpr = np.percentile(tprs, [25, 50, 75], axis=0)
        mean_auc = np.trapz(mean_tpr, mean_fpr)
        plt.plot(mean_fpr, mean_tpr, 'k', lw=3, label='median ROC')
        plt.title(str(mean_auc))
        plt.fill_between(mean_fpr, per_tpr[0, :], per_tpr[2, :], color='g', alpha=.2, label='Quartile')
        plt.legend(loc='lower right')
        plt.savefig(save_dir + '/divided_ROCs1vs2_percentile.pdf')
        del fig
        fig = plt.figure(figsize=(5, 5))
        plt.hist(AUC_set, bins=50)
        plt.savefig(save_dir + '/divided_ROCs1vs2_hist.pdf')
        del fig
        s.close()


        #########################
        fig = plt.figure(figsize=(5, 5))
        plt.plot([0, 1], [0, 1])
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        plt.xlabel('FP')
        plt.ylabel('TP')
        # plt.grid()
        AUC_set = []
        s = open(save_dir + '/divided_RPKM_AUCs1+2.txt', 'w')
        tprs = []
        mean_fpr = np.linspace(0, 1, 100) # 3068
        for jj in range(len(count_set) - 1):  # len(count_set)-1):
            if count_set[jj] < count_set[jj + 1]:
                print(test_indel, jj, count_set[jj], count_set[jj + 1])
                y_test = y_testy[count_set[jj]:count_set[jj + 1]]
                y_predict = y_predicty[count_set[jj]:count_set[jj + 1]]
                # Score trained model.
                y_predict_x = [j[1] + j[2] for j in y_predict]
                y_test_x = [1 - j[0] for j in y_test]
                fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
                tprs.append(interp(mean_fpr, fpr, tpr))
                tprs[-1][0] = 0.0
                # Print ROC curve
                plt.plot(fpr, tpr, color='0.5', lw=0.001, alpha=.2)
                auc = np.trapz(tpr, fpr)
                s.write(str(jj) + '\t' + str(count_set[jj]) + '\t' + str(count_set[jj + 1]) + '\t' + str(auc) + '\n')
                print('AUC:', auc)
                AUC_set.append(auc)

        mean_tpr = np.median(tprs, axis=0)
        mean_tpr[-1] = 1.0
        per_tpr = np.percentile(tprs, [25, 50, 75], axis=0)
        mean_auc = np.trapz(mean_tpr, mean_fpr)
        plt.plot(mean_fpr, mean_tpr, 'k', lw=3, label='median ROC')
        plt.title(str(mean_auc))
        plt.fill_between(mean_fpr, per_tpr[0, :], per_tpr[2, :], color='g', alpha=.2, label='Quartile')
        plt.plot(mean_fpr, per_tpr[0, :], 'g', lw=3, alpha=.2)
        plt.legend(loc='lower right')
        plt.savefig(save_dir + '/divided_ROCs1+2_percentile.pdf')
        del fig
        fig = plt.figure(figsize=(5, 5))
        plt.hist(AUC_set, bins=50)
        plt.savefig(save_dir + '/divided_kegg_ROCs1+2_hist.pdf')
        del fig
        s.close()
    ###########################################################3 if we select the min point as the trained model
    # y_test = y_testy
    # y_predict = y_predicty
    # model.load_weights(save_dir + '/weights.hdf5')
    # scores = model.evaluate(x_test, y_test, verbose=1)
    # print('Test loss:', scores[0])
    # print('Test accuracy:', scores[1])
    # y_predict = model.predict(x_test)
    # np.save(save_dir+'/min_y_test.npy',y_test)
    # np.save(save_dir+'/min_y_predict.npy',y_predict)
    # plt.figure(figsize=(10, 6))
    # for i in range(3):
    #     y_test_x = [j[i] for j in y_test]
    #     y_predict_x = [j[i] for j in y_predict]
    #     fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
    #     plt.subplot(1, 3, i + 1)
    #     plt.plot(fpr, tpr)
    #     plt.grid()
    #     plt.plot([0, 1], [0, 1])
    #     plt.xlabel('FP')
    #     plt.ylabel('TP')
    #     plt.ylim([0, 1])
    #     plt.xlim([0, 1])
    #     auc = np.trapz(tpr, fpr)
    #     print('AUC:', auc)
    #     plt.title('label' + str(i) + ', AUC:' + str(auc))
    # plt.savefig(save_dir + '/min_3labels.pdf')
    # plt.figure(figsize=(10, 6))
    # y_predict_x = [j[1] + j[2] for j in y_predict]
    # y_test_x = [1 - j[0] for j in y_test]
    # fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
    # # Print ROC curve
    # plt.plot(fpr, tpr)
    # plt.plot([0, 1], [0, 1])
    # # Print AUC
    # auc = np.trapz(tpr, fpr)
    # print('AUC:', auc)
    # plt.ylim([0, 1])
    # plt.xlim([0, 1])
    # plt.grid()
    # plt.title('label 1+2, AUC:' + str(auc))
    # plt.xlabel('FP')
    # plt.ylabel('TP')
    # plt.savefig(save_dir + '/min_1+2.pdf')
    # #################################################### data speration
    # plt.figure(figsize=(10, 6))
    # y_predict1 = []
    # y_test1 = []
    # x = 2
    # for i in range(int(len(y_predict) / 3)):
    #     y_predict1.append(y_predict[3 * i][x] - y_predict[3 * i + 1][x])
    #     y_predict1.append(-y_predict[3 * i][x] + y_predict[3 * i + 1][x])
    #     y_test1.append(y_test[3 * i][x])
    #     y_test1.append(y_test[3 * i + 1][x])
    # fpr, tpr, thresholds = metrics.roc_curve(y_test1, y_predict1, pos_label=1)
    # # Print ROC curve
    # plt.plot(fpr, tpr)
    # plt.plot([0, 1], [0, 1])
    # # Print AUC
    # auc = np.trapz(tpr, fpr)
    # print('AUC:', auc)
    # plt.ylim([0, 1])
    # plt.xlim([0, 1])
    # plt.grid()
    # plt.title('label 1 vs 2,direction diff, AUC:' + str(auc))
    # plt.xlabel('FP')
    # plt.ylabel('TP')
    # plt.savefig(save_dir + '/min_1vs2.pdf')
###########################################################################

