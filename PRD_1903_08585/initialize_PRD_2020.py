
import os
import run as r
from dirs import read_dirs as rd
import numpy as np

def correct_ini1_3():

    print('\n')
    print('CORRECTING NORMALIZATION OF EGW IN RUNS ini1, ini2, and ini3',
          'BY A FACTOR 6 pi \n')

    rrs = ['ini1', 'ini2', 'ini3']
    for i in rrs:
        run = runs.get(i)
        EEGW = run.ts.get('EEGW')
        EEGW *= 16*np.pi/3
        run.ts.update({'EEGW':EEGW})

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'
os.chdir(HOME)

# dictionary with the name identifying
# the runs and pointing to the corresponding directory
dirs = rd('PRD_2020_ini')
dirs = rd('PRD_2020_hel', dirs)
dirs = rd('PRD_2020_noh', dirs)
dirs = rd('PRD_2020_ac', dirs)
R = [s for s in dirs]

# set quiet to False to see the spectra available, the runs read,
# and some characteristic info of the run
runs = r.initialize_runs(R, dir0, dirs, quiet=False)
r.characterize_runs(runs, quiet=False)

# runs ini1, ini2, ini3 had an error in the normalization of the time series
# of the GW energy density EEGW, so we correct it here before
# saving the new variables, such that the pickle variables already contain
# the corrected values

#### ONLY RUN THIS ONCE AFTER READING THE DATA FILES!!
correct_ini1_3()

# we now compute the average and maximum spectra over time and
# we study which is the index of the time array after which the spectra
# are no longer growing but just oscillating
# (this can be done separately and one by one for all the runs)

print('Computing the maximum value of the spectrum at each wave number',
      ' over times and the averaged spectrum in the oscillatory regime.',
      ' The resulting indices for the current runs are:\n')
print('   runs: ini1, ini2, ini3, hel1, hel2, hel3, hel4, noh1, noh2,'
      ' ac1, ac2, ac3 \n')
# once we now the indices corresponding to each run we can run all of them
indsst = np.array([3, 50, 300, 10, 20, 80, 15, 10, 30, 11, 10, 20])
print('   indt:', [s for s in indsst], ' \n')

j = 0
for i in runs:
    run = runs.get(i)
    t = run.spectra.get('t_GWs')
    indt = indsst[j]
    j += 1
    run.min_max_stat(indt=indt)
    run.OmGWsat = np.trapz

# we can save the variables runs now to avoid repeating the
# computation of the averaged and maximum spectra
for i in runs:
    run = runs.get(i)
    run.save(dir0=dir0)
