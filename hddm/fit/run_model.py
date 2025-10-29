#####
#####  RUN NOTES
#####

# HDDM analysis
# Written by Jae-Young Son (2019)
# Adapted by Harrison Ritz (2025)



#####
#####  SET UP ENVIRONMENT
#####

# Import requisite packages
import matplotlib
matplotlib.use('Agg') # Run before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import hddm
import pandas as pd
import pickle
# from patsy import dmatrix
import kabuki
# from kabuki.analyze import gelman_rubin
# from kabuki.utils import concat_model
# from kabuki.analyze import post_pred_gen
import pathlib
import numpy as np
import sys
import os

# get model version
modelVersion = '24-07-04_NN_'
modelName = sys.argv[1]

# os.chdir('/mnt/')

# print info
print(' ')
print(' ----- ')
print('CURRENT WORKING DIRECTORY')
print(os.getcwd())
print(' ----- ')
print(' ')

print(' ')
print(' ----- ')
print('PYTHON VERSION')
print(sys.version)
print(' ----- ')
print(' ')

print(' ')
print(' ----- ')
print('HDDM VERSION')
print(hddm.__version__)
print(' ----- ')
print(' ')

print(' ')
print(' ----- ')
print('MODEL VERSION')
print(modelVersion)
print(' ----- ')
print(' ')

print(' ')
print(' ----- ')
print('MODEL NAME')
print(modelName)
print(' ----- ')
print(' ')


# Get around a problem with saving regression outputs in Python 3
def savePatch(self, fname):
    import pickle
    with open(fname, 'wb') as f:
        pickle.dump(self, f)
hddm.HDDM.savePatch = savePatch


# Load data from csv file into a NumPy structured array
if modelName == 'VD_orig':
    data = hddm.load_csv('./data/vcr_orig.csv')
else:
    data = hddm.load_csv('./data/vcr_respCoded.csv')

print(data.head())
# for col in data.columns:
#     print(col)
print(data.shape)



# Check whether save directories exist; if not, create them
pathlib.Path('./models/' + modelVersion + 'model/').mkdir(parents=True, exist_ok=True)
pathlib.Path('./results/' + modelVersion + 'results/').mkdir(parents=True, exist_ok=True)
pathlib.Path('./ppc/'  + modelVersion + 'ppc/').mkdir(parents=True, exist_ok=True)
pathlib.Path('./rhat/'  + modelVersion + 'rhat/').mkdir(parents=True, exist_ok=True)


#####
#####  RUN MODEL
#####


# parameter setting tutorial: https://hddm.readthedocs.io/en/latest/tutorial_set_parameter_defaults.html
# singularity exec --bind $PWD/:/mnt/ ~/data/VCR/hddm.sif python3 /mnt/run_model.py 'VD_collapse'
print("")
print("NOTE: modified to use HDDM 1.0.1")
print("BASIC CONFIGURATION:")
print(hddm.model_config.model_config['ddm_hddm_base'])
print("ANGLE CONFIGURATION:")
print(hddm.model_config.model_config['angle'])
print("")

base_config = hddm.model_config.model_config['angle']

# fix theta to 0
base_config['params_default'][4] = 0.0

# increase parameter range
# v range
base_config['param_bounds'][0][0] = -4.0
base_config['param_bounds'][1][0] =  4.0
# a range
base_config['param_bounds'][0][1] = .250
base_config['param_bounds'][1][1] =  4.0
# z range
base_config['param_bounds'][0][2] = .05
base_config['param_bounds'][1][2] = .95

print("modified base_config to set theta=0")
print(base_config)

print("Running model: " + modelName)


models = []

for i in range(5):
   

    # SUBJECT-LEVEL MODELS 


    if modelName == 'VDev_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + signVD', 'a ~ 1 + absVD'], 
            include = {'v','a','t','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )
        

    elif modelName == 'VDev_sAll_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + signVD', 'a ~ 1 + absVD'], 
            include = {'v','a','t','z', 'sv','st','sz'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )
        
    elif modelName == 'VDev_collapse_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + signVD', 'a ~ 1 + absVD', 'theta ~ 1 + absVD'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )
        
    elif modelName == 'VDOVev_collapse_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + signVD', 'a ~ 1 + absVD + OV', 'theta ~ 1 + absVD + OV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )








    elif modelName == 'VD_static_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD'], 
            include = {'v','a','t','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )


    elif modelName == 'VD_staticVar_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD'], 
            include = {'v','a','t','z','sv','st','sz'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )


    elif modelName == 'VD_collapse_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD', 'theta ~ 1 + absMaxVD + absMinVD'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'VD_collapseA_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'VD_collapseT_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'theta ~ 1 + absMaxVD + absMinVD'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )








    elif modelName == 'OV_static_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + maxOV + minOV'], 
            include = {'v','a','t','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )


    elif modelName == 'OV_collapse_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + maxOV + minOV', 'theta ~ 1 + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'OV_collapseA_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'OV_collapseT_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'theta ~ 1 + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )







    elif modelName == 'VDOV_static_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
            include = {'v','a','t','z'},
            informative = False,
            group_only_regressors=False,
            model='angle',
            model_config=base_config,
            )


    elif modelName == 'VDOV_collapse_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD + maxOV + minOV', 'theta ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'VDOV_collapseA_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )


    elif modelName == 'VDOV_collapseT_subj':

        m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'theta ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )
        





    # sample
    print('finding starting values...')
    m.find_starting_values(cycles=2)

    print('fiting model...')
    m.sample(6000, burn=2000,
        dbname='./models/' + modelVersion + 'model/' + modelName + '_%s.db'%i,
        db='pickle')
    # m.sample(10, burn=5,
    #     dbname='./models/' + modelVersion + '/' + modelName + '_%s.db'%i,
    #     db='pickle')
    print('model fit')
    models.append(m)



# print stats
catM = kabuki.utils.concat_models(models)
catM.print_stats()


# rhat
pd.DataFrame.from_dict(kabuki.analyze.gelman_rubin(models), orient='index').to_csv('./rhat/' + modelVersion + 'rhat/' + modelName + '_rhat.csv')

# Save traces
traces = catM.get_traces()
traces.to_csv('./results/' + modelVersion + 'results/' + modelName + '_traces.csv')

# save summary
summary_df = traces.mean(axis=0)
summary_df['DIC'] = catM.dic
print(summary_df)
summary_df.to_csv('./results/' + modelVersion + 'results/' + modelName + '_summary.csv')


# PPC
print('generating PPC')
kabuki.analyze.post_pred_gen(catM, samples=1000, append_data=True).to_csv('./ppc/' + modelVersion + 'ppc/' + modelName + '_ppc.csv')