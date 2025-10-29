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




# print info
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





modelVersion = '24-07-04_NN_'








# get model version ============================================================
model_idx =  int(sys.argv[1])

print(' ')
print(' ----- ')
print('MODEL IDX')
print(model_idx)
print(' ----- ')
print(' ')



# set generative model
if (model_idx > 0) & (model_idx <= 200):
    generating_model = 'VDev_subj'

elif (model_idx > 200) & (model_idx <= 400):
    generating_model = 'VD_static_subj'

elif (model_idx > 400) & (model_idx <= 600):
    generating_model = 'VD_collapse_subj'

elif (model_idx > 600) & (model_idx <= 800):
    generating_model = 'VDOV_collapseA_subj'

else:
    print('ERROR: model_idx out of range for generating_model')
    sys.exit()


print(' ')
print(' ----- ')
print('GENERATING MODEL')
print(generating_model)
print(' ----- ')
print(' ')



# fitted model name
if (np.mod(model_idx, 200) > 0) & (np.mod(model_idx, 200) <= 50):
    fitted_model = 'VDev_subj'

elif (np.mod(model_idx, 200) > 50) & (np.mod(model_idx, 200) <= 100):
    fitted_model = 'VD_static_subj'

elif (np.mod(model_idx, 200) > 100) & (np.mod(model_idx, 200) <= 150):
    fitted_model = 'VD_collapse_subj'

elif (np.mod(model_idx, 200) > 150) | (np.mod(model_idx, 200) == 0):
    fitted_model = 'VDOV_collapseA_subj'

else:
    print('ERROR: model_idx out of range for fitting model')
    sys.exit()


print(' ')
print(' ----- ')
print('FITTED MODEL')
print(fitted_model)
print(' ----- ')
print(' ')


# sample number
sample_number = np.mod(model_idx, 50)
if sample_number == 0:
    sample_number = 50


print(' ')
print(' ----- ')
print('SAMPLE NUMBER')
print(sample_number)
print(' ----- ')
print(' ')

# ==============================================================================





model_name = modelVersion + 'gen-' + generating_model + '__fit-' + fitted_model
print(model_name)


# Get around a problem with saving regression outputs in Python 3
def savePatch(self, fname):
    import pickle
    with open(fname, 'wb') as f:
        pickle.dump(self, f)
hddm.HDDM.savePatch = savePatch


# Load data from csv file into a NumPy structured array
data = hddm.load_csv(modelVersion + 'ppc/' + generating_model +'_ppc.csv')
print('data loaded from: ' + generating_model +'_ppc.csv')

data = data[data['sample'] == int(sample_number)]
data["rt"] = data["rt_sampled"]
data["response"] = data["response_sampled"]
data = data[data['rt'] > .100]
data = data[data['rt'] < .750] # remove responses past deadline

data.drop(['node', 'Unnamed: 2'], axis=1, inplace=True)

print(data.head())
# for col in data.columns:
#     print(col)
print(data.shape)




# Check whether save directories exist; if not, create them
pathlib.Path(modelVersion + '/models/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/results/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/ppc/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/rhat/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/recovery/stats/' + model_name + '/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/recovery/traces/' + model_name + '/').mkdir(parents=True, exist_ok=True)
pathlib.Path(modelVersion + '/recovery/summary/' + model_name + '/').mkdir(parents=True, exist_ok=True)


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
base_config['param_bounds'][0][0] = -10.0
base_config['param_bounds'][1][0] =  10.0
# a range
base_config['param_bounds'][0][1] = .250
base_config['param_bounds'][1][1] =  10.0
# z range
base_config['param_bounds'][0][2] = .05
base_config['param_bounds'][1][2] = .95

print("modified base_config to set theta=0")
print(base_config)



models = []

   



# SUBJECT-LEVEL MODELS


if fitted_model == 'VDev_subj':

    m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + signVD', 'a ~ 1 + absVD'], 
            include = {'v','a','t','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle',
            model_config=base_config,
            )


elif fitted_model == 'VD_static_subj':

    m = hddm.HDDMnnRegressor(data, 
        ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD'], 
        include = {'v','a','t','z'},
        informative = False,
        group_only_regressors=False,
        model = 'angle',
        model_config=base_config,
        )


elif fitted_model == 'VD_collapse_subj':

    m = hddm.HDDMnnRegressor(data, 
            ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
            include = {'v','a','t','theta','z'},
            informative = False,
            group_only_regressors=False,
            model = 'angle'
            )



elif fitted_model == 'VDOV_collapseA_subj':

    m = hddm.HDDMnnRegressor(data, 
        ['v ~ -1 + maxVD + minVD', 'a ~ 1 + absMaxVD + absMinVD + maxOV + minOV'], 
        include = {'theta', 'z', 'p_outlier'},
        informative = False,
        group_only_regressors=False,
        model = 'angle'
        )






# sample
print('finding starting values...')
m.find_starting_values()

print('starting sampling...')
m.sample(12000, burn=2000,
    dbname= modelVersion + '/models/' + fitted_model + '_%s.db'%int(sys.argv[1]),
    db='pickle')

print('finished sampling...')




# save stats
print("STATS ================================================================")
m.print_stats()
m.print_stats(fname= modelVersion + '/recovery/stats/' + model_name + '/' + 'gen-' + generating_model + '__fit-' + fitted_model + '__sample-' + str(sample_number) + '_stats.csv')


# Save traces
traces = m.get_traces()
traces.to_csv( modelVersion + '/recovery/traces/' + model_name + '/' + 'gen-' + generating_model + '__fit-' + fitted_model + '__sample-' + str(sample_number) + '_traces.csv')


# save summary
summary_df = traces.mean(axis=0)
summary_df['DIC'] = m.dic
print(summary_df.head())
summary_df.to_csv( modelVersion + '/recovery/summary/' + model_name + '/' + 'gen-' + generating_model + '__fit-' + fitted_model + '__sample-' + str(sample_number) + '_summary.csv')