# Model Training Imports
import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping, ModelCheckpoint

# System Access Imports
import sys, os 

# Matplotlib and Numpy Imports
import numpy as np 
import matplotlib
matplotlib.use('Agg') # Matplotlib settings
import matplotlib.pyplot as plt 

# Other
from sklearn import metrics
from scipy import interp

### Parameter settings
data_augmentation = False                           # Look into hard-coding, this is only used once
batch_size = 1024
epochs = 20
model_name = "keras_css_trained_model_shallow.h5"   # May change to take user cmdline argument for model name

def load_data_TF2(indel_list: list[int], data_path: str):
    '''
    Loads the data using the indices passed in from indel_list
    to access the files stored under the directory `data_path`.

    Arguments:
        indel_list: List[int]
            - A list of indices
        data_path: str
            - A string describing the path to the .npy files
    
    Returns:
        A 1-element tuple which is a 3-tuple:
            1. A numpy array of the data - np.array(xxdata_list)
            2. A numpy array of the data - yydata_x = np.array(yydata).astype('int')
            3. The set of counts for each index - count_set
    '''
    import random
    xxdata_list = []
    yydata = []
    count_set = [0]
    count_setx = 0
    for i in indel_list:
        # Loads data from the respective .npy files stored in the data_path folder
        xdata = np.load(f"{data_path}/Nxdata_tf{i}.npy")
        ydata = np.load(f"{data_path}/ydata_tf{i}.npy")
        for k in range(len(ydata)):
            xxdata_list.append(xdata[k,:,:,:])
            yydata.append(ydata[k])
        count_setx = count_setx + len(ydata)
        count_set.append(count_setx)
        # Info-printing purpose
        print(f"Index {i} has length {len(ydata)}.")
    yydadta_array = np.array(yydata)
    yydata_x = yydata_array.astype('int')
    # Info-printing purpose
    print(np.array(xxdata_list).shape)
    return ((np.array(xxdata_list), yydata_x, count_set))

def get_model(x_train):
    '''
    Creates the base model and returns it.

    Arguments:
        x_train
            - The np.array of training data returned from load_data_TF2
    
    Returns:
        The base model used for training
    '''
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

    return model

def run(length_TF: int, data_path: str, num_classes: int):
    '''
    Gets the input data and creates/trains/evaluates the model.
    This function also saves the model to the predefined location
    given by `model_name`.

    Arguments:
        length_TF: int
            The number of data parts divided
        data_path: str
            The string containing the path to the directory of the .npy files
        num_classes: int
            The number of categories

    Returns:
        None
    '''
    whole_data_TF = [i for i in range(length_TF)]

    # Data retrieval:
    for test_indel in range(3):
        # Questionable splitting: This literally does three iterations:
        # one iteration per third of the indices with no randomness
        # whole_data_TF = [ Segment-1 | Segment-2 | Segment-3 ] becomes
        # Test: Segment-1, Train: Segment-2+3
        # Test: Segment-2, Train: Segment-1+3
        # Test: Segment-3, Train: Segment-1+2
        test_TF = [i for i in range(int(np.ceil(0.333 * test_indel * length_TF)),
                                    int(np.ceil(0.333 * (test_indel + 1) * length_TF)))]
        train_TF = [i for i in whole_data_TF if i not in test_TF]
        (x_train, y_train, count_set_train) = load_data_TF2(train_TF, data_path)
        (x_test, y_test, count_set) = load_data_TF2(test_TF, data_path)

        # Info-printing purpose
        print(f"{x_train.shape} x_train samples")
        print(f"{x_test.shape} x_test samples")

        # Save Directory
        save_dir = os.path.join(os.getcwd(), 
                                f"{test_indel}YYYY_saved_models_T_32-32-64-64-128-128-512_e{epochs}")

        if num_classes > 2:
            y_train = keras.utils.to_categorical(y_train, num_classes)
            y_test = keras.utils.to_categorical(y_test, num_classes)
        
        # Info-printing purpose
        print(f"{y_train.shape} y_train samples")
        print(f"{y_test.shape} y_test samples")

        # Make the save directory:
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        # Exit and fail if num_classes < 2:
        if num_classes < 2:
            print("Not enough categories!")
            sys.exit()

        # Get the base model:
        model = get_model(x_train)
        sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        if num_classes == 2:
            model.add(Dense(1, activation="sigmoid"))
            model.compile(optimizer=sgd, loss='binary_crossentropy', metrics=['accuracy'])
        else:
            model.add(Dense(num_classes))
            model.add(Activation('softmax'))
            model.compile(optimizer=sgd, loss='categorical_crossentropy', metrics=['accuracy'])
        # Early Stop Points
        early_stopping = EarlyStopping(monitor='val_acc', patience=50, verbose=0, mode='auto')
        # Checkpoints
        checkpoint1 = ModelCheckpoint(filepath=f"{save_dir}/weights.{epoch:02d}-{val_loss:.2f}.hdf5",
                                      monitor='val_loss', verbose=1, save_best_only=False, 
                                      save_weights_only=False, mode='auto', period=1)
        checkpoint2 = ModelCheckpoint(filepath=f"{savedir}/weights.hdf5", monitor='val_acc', verbose=1,
                                      save_best_only=True, mode='auto', period=1)
        callbacks_list = [checkpoint2, early_stopping]      # Not sure why checkpoint1 exists
        
        # Train (unconditionally, apparently, since data_augmentation was hard-coded as False)
        history = None
        if not data_augmentation:
            # Info-printing purposes
            print("Not using data augmentation.")
            # Complete the training (fitting)
            history = model.fit(x_train, y_train, 
                                batch_size=batch_size, epochs=epochs,
                                validation_split=0.2, shuffle=True,
                                callbacks=callbacks_list)
        
        # Save the model and weights
        model_path = os.path.join(save_dir, model_name)
        model.save(model_path)
        
        # Info-printing purpose
        print(f"Saved the trained model at {model_path}!")

        # Evaluate the trained model
        scores = model.evaluate(x_test, y_test, verbose=1)

        # Info-printing purpose
        print(f"Test loss: {scores[0]}")
        print(f"Test accuracy: {scores[1]}")

        # Validation?
        y_predict = model.predict(x_test)

        # Save test data and prediction data
        np.save(f"{save_dir}/end_y_test.npy", y_test)
        np.save(f"{save_dir}/end_y_predict.npy", y_predict)

        # Optionally (unconditionally for now) included for results:
        plt.figure(figsize=(10, 6))
        plt.subplot(1, 2, 1)
        plt.plot(history.history['acc'])
        plt.plot(history.history['val_acc'])
        plt.title("Model Accuracy")
        plt.ylabel("Accuracy")
        plt.xlabel("Epoch")
        plt.grid()
        plt.legend(["train", "val"], loc='upper left')
        plt.subplot(1, 2, 2)
        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title("Model Loss")
        plt.ylabel("Loss")
        plt.ylabel("Epoch")
        plt.legend(["train", "val"], loc='upper left')
        plt.grid()
        plt.savefig(f"{save_dir}/end_result.pdf")