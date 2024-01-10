####################################################################################################################################################################################################
##### IMPORTS
####################################################################################################################################################################################################


from scipy.stats import entropy
import scipy
from scipy.linalg import sqrtm
from sklearn.metrics.pairwise import euclidean_distances
from tensorflow.keras.preprocessing.sequence import pad_sequences
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
# Tensorflow / Keras
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from math import sqrt
import os
import glob
# for building Neural Networks
from tensorflow import keras
# for model evaluation
import sklearn
from sklearn.manifold import TSNE
from sklearn import metrics
# For performing clustering
from sklearn.cluster import DBSCAN
# For rescaling metrics to fit into 0 to 1 range
from sklearn.preprocessing import MinMaxScaler, LabelEncoder

# from hdbscan import HDBSCAN
# from hdbscan.flat import (HDBSCAN_flat,
#                           approximate_predict_flat,
#                           membership_vector_flat,
#                           all_points_membership_vectors_flat)

import seaborn as sns

import argparse
import os
import glob
import random

import warnings
warnings.filterwarnings('ignore')
print('sklearn: %s' % sklearn.__version__)  # print version
print('Tensorflow/Keras: %s' % keras.__version__)  # print version






####################################################################################################################################################################################################
## VAE FUNCTION
#####################################################################################################################################################################################################

tfk = tf.keras
tfkl = tfk.layers

# Defining input and output
# you can modify the input and output path here as well.
# config = {
#     'input_file_location': '/home/jasmine/Data/patient_1733321.csv',
#     'responsible_variants_column_name': 'responsible_variants',
#     'output_file_location': '/home/justine/'}
# parser = argparse.ArgumentParser()
# parser.add_argument('--ds', help='Input file path to be entered', type=str)

def seed_everything(seed=779):
    """"
    Seed everything.
    """
    tf.random.set_seed(seed)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    
def make_features(data):
    """
    This function performs feature selection and normalization by keeping only the important features for modelling.
    Note: it might throw an error if the name of features differs from the function and the dataset.
    :param data: input file
    :return: returns a normalized features in a dataframe format
    """
    
    #Feature list for hg37 genome
    features_data = data[
        ['ConsScore', 'SIFTval', 'PolyPhenVal', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
         'verPhyloP', 'EncH3K27Ac', 'EncH3K4Me1', 'EncH3K4Me3', 'EncExp', 'EncNucleo', 'EncOCC', 'EncOCCombPVal',
         'EncOCDNasePVal', 'EncOCFairePVal', 'EncOCpolIIPVal', 'EncOCctcfPVal', 'EncOCmycPVal', 'EncOCDNaseSig',
         'EncOCFaireSig', 'EncOCpolIISig', 'EncOCctcfSig', 'EncOCmycSig', 'SpliceAI.acc.gain', 'SpliceAI.acc.loss',
         'SpliceAI.don.gain'
            , 'SpliceAI.don.loss', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
         'MMSp_donorIntron',
         'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp',
         'Sngl10000bp', 'RawScore'
            , 'PHRED']]

    

    # #Feature list for hg38 genome
    # features_data = data[
    #      ['priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
    #         'verPhyloP', 'EncodeH3K27ac.max', 'EncodeH3K27ac.sum','EncodeH3K4me1.max', 'EncodeH3K4me1.sum','EncodeH3K4me3.sum',
    #         'EncodeH3K4me3.max', 'EncodeDNase.sum','EncodeDNase.max', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
    #         'MMSp_donorIntron', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp',
    #           'Rare10000bp', 'Sngl10000bp', 'RawScore', 'PHRED']]
    
    
    lb_make = LabelEncoder()
    #features_data['Type'] = lb_make.fit_transform(features_data['Type'])

    features_data_treated = features_data.fillna(features_data.median())

    scaler = MinMaxScaler()
    print(scaler.fit(features_data_treated))

    scaled_data = scaler.transform(features_data_treated)
    scaled_dataframe = pd.DataFrame(scaled_data, columns=features_data_treated.columns)

    return scaled_dataframe

#attention, change epoch pour test mais original value : 500

def Vae_function(file, read_count=True,
                 batch_size=16, epochs=50, kernel="lecun_normal",
                 kl_loss_weight=0.5, edl1=25,
                 latent_size=16, ddl1=25):
    """
     This function performs variational autoencoder architecture
    :param file: the output file from the make_features function
    :param batch_size: Batch size defined according to dataset size
    :param epochs: You can define more /less number of epochs depending on the data
    :param edl1: encoder layer 1 , defined based on the number of features
    :param edl2: encoder layer 2
    :param edl3: encoder layer 3
    :param edl4: encoder layer 4
    :param latent_size: defined for 16 dimensions, can be changed later incase needed
    :param ddl1: decoder layer 1
    :param ddl2: decoder layer 2
    :param ddl3: decoder layer 3
    :param ddl4: decoder layer 4
    :return: returns prediction , latent space and history of the models
    """
    
    seed_everything()
    def sampling(args):
        z_mean, z_log_var = args
        batch = tfk.backend.shape(z_mean)[0]
        dim = tfk.backend.int_shape(z_mean)[1]
        epsilon = tfk.backend.random_normal(shape=(batch, dim), seed=779)
        out = z_mean + tf.keras.backend.exp(0.5 * z_log_var) * epsilon
        return out

    data = file

    n_cols = data.shape[1]
    start_time = time.time()

    #####ENCODER######

    # Input
    inputs = tfk.Input(shape=(n_cols,), name='encoder_input')

    # First hidden layer
    hidden_encoder_1 = tfkl.Dense(edl1, kernel_initializer=kernel, name="hidden_encoder_1")(inputs)
    encoder_norm_1 = tfkl.BatchNormalization(name="encoder_norm_1")(hidden_encoder_1)
    hidden_encoder_1_activation = tfkl.ELU()(encoder_norm_1)

    # Mean for the sampling
    z_mean = tfkl.Dense(latent_size, name='z_mean')(hidden_encoder_1_activation)

    # Var for the sampling
    z_log_var = tfkl.Dense(latent_size, name='z_log_var')(hidden_encoder_1_activation)

    # Sample from the values belowfrom sklearn.model_selection import train_test_split
    z = tfkl.Lambda(sampling, output_shape=(latent_size,), name='z')([z_mean, z_log_var])

    # Encoder model
    encoder = tfk.Model(inputs, [z_mean, z_log_var, z], name='encoder')
    # encoder.save('/home/jasmine/Data/VAEncoder_rare.h5')

    #####DECODER#####

    # Input
    latent_inputs = tfk.Input(shape=(latent_size,), name='z_sampling')

    # Fifth hidden layer
    hidden_decoder_5 = tfkl.Dense(ddl1, kernel_initializer=kernel, name="hidden_decoder_4")(latent_inputs)
    decoder_norm_5 = tfkl.BatchNormalization(name="decoder_norm_4")(hidden_decoder_5)
    hidden_encoder_5_activation = tfkl.ELU()(decoder_norm_5)

    # Output
    outputs = tfkl.Dense(n_cols)(hidden_encoder_5_activation)

    # Decoder model
    decoder = tfk.Model(latent_inputs, outputs, name="decoder")

    #####VAE####
    outputs = decoder(encoder(inputs)[2])
    vae = tfk.Model(inputs, outputs=[outputs, z], name="vae")

    #####Personalized loss
    msle = tfk.losses.MSLE(inputs, outputs[0])
    msle *= n_cols
    kl_loss = 1 + z_log_var - tfk.backend.square(z_mean) - tfk.backend.exp(z_log_var)
    kl_loss = tfk.backend.sum(kl_loss, axis=-1)
    kl_loss *= -kl_loss_weight
    vae.add_loss(msle + kl_loss)
    vae.add_metric(msle, name='msle', aggregation='mean')
    vae.add_loss(kl_loss)
    vae.add_metric(kl_loss, name='kl_loss', aggregation='mean')
    vae.compile(optimizer="adam")#, loss='mean_squared_error', metrics=['MeanSquaredError'] ,experimental_run_tf_function=False)
    vae.summary()

    #####Train
    history = vae.fit(data, epochs=epochs, batch_size=batch_size)
    end_time = time.time() - start_time
    print("-----", end_time // 3600, "hours", (end_time % 3600) // 60, "minutes", round((end_time % 3600) % 60, 2),
          "secondes", "-----")
    # To save the latent space
    pred, latent_space = vae.predict(data)
    latent_space = pd.DataFrame(latent_space)

    pred = pd.DataFrame(pred)
    pred.index = data.index
    pred.columns = data.columns

    
    # generate the samples using the VAE model
    latent_samples = np.random.normal(size=(len(data), 16))
    generated_samples = decoder.predict(latent_samples)

    return latent_space


####################################################################################################################################################################################################
##DBSCAN
####################################################################################################################################################################################################

### aff function takes only the 1st 16 columns of the latent space from the data
def aff (data):
    data_db = data.iloc[:, :16]
    return (data_db)

### finds the right epsilon value following the patient data that is used   
def find_eps(d, n=2500):
    # Calculate the average distance to k-nearest neighbors for each point
    k = int(sqrt(d.shape[0])) # You can adjust this value
    neigh = NearestNeighbors(n_neighbors=k)
    neigh.fit(d)
    distances, indices = neigh.kneighbors(d)
    avg_distances = np.mean(distances, axis=1)

    # Sort the array of average distances and save the y values
    avg_distances_sorted = np.sort(avg_distances)
    y_values = avg_distances_sorted
    # Sort the array of average distances and save the y values
    avg_distances_sorted = np.sort(avg_distances)
    y_values = avg_distances_sorted

   
    # Extract the y values where x is between 0 and 1
    x_values = np.arange(len(y_values))
    y_values_subset = y_values[(x_values >= 0) & (x_values <= n)]
    eps =  sum(y_values_subset) / len(y_values_subset)
    return eps  



### dbs function applies DBSCAN clustering for the patients with identified responsible variants.
### m is the value of min_samples and it should be 2 x latent_size
def dbs (d, e = 10, m = 32):
    scaler = StandardScaler()
    db = DBSCAN(eps=e, min_samples=m)
    d_db = aff(d)
    d_std = scaler.fit_transform(d_db)
    model= db.fit(d_std)
    clusters_d = pd.DataFrame(model.fit_predict(d_std))
    d["Cluster"] = clusters_d
    print(d["Cluster"].value_counts())
    print(e)
    return (d)




####################################################################################################################################################################################################
##MAIN
####################################################################################################################################################################################################

if __name__ == '__main__':
    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file_path', type=str, required=True,
                        help='Path to the file to process')
    parser.add_argument('-o', '--output_folder_path', type=str, required=True,
                    help='Path to the output folder containing the results')
    args = parser.parse_args()

    # Use the specified parent folder path
    filepath = args.file_path
    output = args.output_folder_path
    
    dim = pd.DataFrame()
    dim['patient'] = ''
    dim['n_variant'] = ''
    dim['n_variant_filtered'] = ''
    dim['n_variant_ls'] = ''
    #dim['epsilon'] = ''
    #dim['n_variant_dbs'] = ''

    # Loop through the files and do whatever you need to do with them
    
    #Feature list for hg37 genome
    list_features = ['ConsScore', 'SIFTval', 'PolyPhenVal', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
         'verPhyloP', 'EncH3K27Ac', 'EncH3K4Me1', 'EncH3K4Me3', 'EncExp', 'EncNucleo', 'EncOCC', 'EncOCCombPVal',
         'EncOCDNasePVal', 'EncOCFairePVal', 'EncOCpolIIPVal', 'EncOCctcfPVal', 'EncOCmycPVal', 'EncOCDNaseSig',
         'EncOCFaireSig', 'EncOCpolIISig', 'EncOCctcfSig', 'EncOCmycSig', 'SpliceAI.acc.gain', 'SpliceAI.acc.loss',
         'SpliceAI.don.gain'
            , 'SpliceAI.don.loss', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
         'MMSp_donorIntron',
         'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp',
         'Sngl10000bp', 'RawScore'
            , 'PHRED']

    # #Feature list for hg38 genome
    # list_features = ['priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
    #              'verPhyloP', 'EncodeH3K27ac.max', 'EncodeH3K27ac.sum','EncodeH3K4me1.max', 'EncodeH3K4me1.sum','EncodeH3K4me3.sum',
    #              'EncodeH3K4me3.max', 'EncodeDNase.sum','EncodeDNase.max', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
    #              'MMSp_donorIntron', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp',
    #              'Rare10000bp', 'Sngl10000bp', 'RawScore', 'PHRED']


    if filepath.endswith('.csv'):
        
        #extract the name of the patient
        filename = filepath.split("/")[-1]
        name= filename.split('_')[0]
        print(name)
        
        data = pd.read_csv(filepath)
        
        #drop the duplicate variants
        d = data[~data[list_features].duplicated()]
        d.reset_index(drop=True, inplace=True)
    
        #drop the variants with no associated Gene
        d = d.dropna(subset=['GeneName'])
        d.reset_index(drop=True, inplace=True)
        
        
        #Running the VAE
        ls = Vae_function(make_features(d))
        
        #Latent space result
        res = pd.concat([ls, d], axis=1, join='inner')
        
        res.to_csv(output+name + "_rare_variant_latent_space_filtered.csv", index=False)
        
        #DBSCAN result
        e = find_eps(aff(res), 1500)
        res_db = dbs(res, e)
        res_db = res_db[res_db['Cluster'] == -1]

        res_db.to_csv(output+name +'_res_dbscan.csv', index=False)
        
        #saving the dimensions of each result
        dim.loc[len(dim)]=[name, data.shape[0], d.shape[0],res.shape[0]]#, e, res_db[res_db['Cluster'] == -1].shape[0]]
        
        dim.to_csv(output+name +'_dim.csv', index=False)
    else:
        print('File needs to be in CSV format')
