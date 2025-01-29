Contains all the testing for computing the mineco among frame pairs from ANM 
This version is made to run on PHX supercomputing cluster at ASU

### Updates
## 1
''' 
sbatch md.sh
'''

Contains the setup.sh script call to avoid manual loading

## 2
eco_anm.yml has the appropriate packages - with versions - and can be used to create the required python environment:
'''
mamba create -f eco_anm.yml
'''
