# ======================================================
#       Author :      Jacob Miratsky
#                     Chemistry Ph.D.
#       Affiliation : ASU, School of Molecular Sciences
#       Contact :     jacoba1@asu.edu
#       Usage : Change 'filename'
# ======================================================
import os
import re
import math
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.cluster import DBSCAN
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import add_distance_to_edges, add_peptide_bonds, add_distance_threshold
from graphein.protein.graphs import construct_graph
import sys

class PerFrameANMECO:

    def __init__(self, starting_dir):
        self.starting_dir = starting_dir

    def run_eco(self, pdb1, pdb2, override=0):

        # This script will run minECO. Provide your two PDBs as
        # command line arguments. It will auto-detect which PDB
        # has fewer contacts and make that the starting one. If
        # you want to override this behavior, a third command line
        # argument is provided to "flip" the ordering. If the
        # override is set to 0, then it will go from the structure
        # with fewer contacts to the one with more. If the override
        # is 1, it will go in the reverse direction.

        # Run with python eco.py PDB1.pdb PDB2.pdb 0
        # or with python eco.py PDB1.pdb PDB2.pdb 1

        # Output is written to eco.log.

        f = open(f"eco_{pdb1}_{pdb2}.log", "w")

        # protein1 = sys.argv[1]
        # protein2 = sys.argv[2]
        # override = int(sys.argv[3])
        protein1 = pdb1
        protein2 = pdb2

        # obtain contact graphs
        config = ProteinGraphConfig(**{"granularity": "CA", "edge_construction_functions": [add_peptide_bonds]})

        P1 = construct_graph(config=config, path=protein1)
        #add_distance_threshold(P1, long_interaction_threshold=6, threshold=7)
        add_distance_threshold(P1, long_interaction_threshold=1, threshold=5)
        add_distance_to_edges(P1)

        P2 = construct_graph(config=config, path=protein2)
        add_distance_threshold(P2, long_interaction_threshold=1, threshold=5)
        add_distance_to_edges(P2)

        P1_edges = P1.number_of_edges()
        P2_edges = P2.number_of_edges()

        if override == 0:
            if P1_edges > P2_edges:
                C = P1
                O = P2
            elif P2_edges > P1_edges:
                C = P2
                O = P1
        elif override != 0:
            if P1_edges > P2_edges:
                C = P2
                O = P1
            elif P2_edges > P1_edges:
                C = P1
                O = P2

        # identify unique contacts
        U = C.copy()
        U.remove_edges_from(n for n in C.edges if n in O.edges)

        # filter unique contacts by distance change
        pairs = []
        O2 = O.copy()
        for a, b in U.edges:
            O2.add_edge(a,b)

        add_distance_to_edges(O2)

        for a, b in U.edges:
            if np.abs(O2.edges[a,b]["distance"]-C.edges[a,b]["distance"]) < 1:
                U.remove_edge(a,b)

        # cluster unique contacts
        for a, b, d in U.edges(data=True):
            ri1 = int(re.findall(r'\d+', a)[0])
            ri2 = int(re.findall(r'\d+', b)[0])
            pairs.append([ri1, ri2])

        #clustering = DBSCAN(eps = 6, min_samples = 1).fit(np.array(pairs))
        clustering = DBSCAN(eps = 1, min_samples = 1).fit(np.array(pairs))
        labels = np.array(clustering.labels_)

        clusters = np.unique(labels)

        totalECO = 0

        min_eco_clusters = []

        min_eco_contacts = []

        while clusters.any():

            allContacts = {}
            allECOs = {}
            minECOs = []
            allENEs = []
            for cluster in clusters:
                clusterContacts = []
                clusterECOs = []
                cluster_edges = np.array(U.edges)[labels==cluster]
                for a, b in cluster_edges:
                    clusterContacts.append([a, b])
                    clusterECOs.append(nx.shortest_path_length(O, a, b))
                allContacts[cluster] = clusterContacts
                allECOs[cluster] = clusterECOs
                minECOs.append(min(clusterECOs))

            minCluster = clusters[np.argmin(minECOs)]
            min_eco_clusters.append(minCluster)
            min_eco_contacts.append(allContacts[minCluster])
            minCluster_edges = np.array(U.edges)[labels==minCluster]
            for mce in minCluster_edges:
                O.add_edge(mce[0], mce[1])
                pair = mce[0].split(":")[1]+","+mce[1].split(":")[1]
            f.write("formed cluster" + " " + str(minCluster) + " " + "with contacts \n" + 
                    str(np.column_stack([allContacts[minCluster], allECOs[minCluster]])) + "\n")

            totalECO += min(minECOs)
            clusters = np.delete(clusters, np.where(clusters==minCluster), 0)
                
        f.write("total ECO for this pathway:" + str(totalECO) + "\n")

        f.write("the minECO pathway (by cluster):" + str(min_eco_clusters) + "\n")

        f.write("the minECO pathway (by contact):" + str(min_eco_contacts) + "\n")

        f.close()

        return totalECO, min_eco_contacts
        

    def iterative_eco(self):
        ''' 
        This function will take pathway.pdb and compute the eco score between all the frames.

        Returns:
        --------
        .csv file with the following columns:
        name, total_eco, min_eco_clusters, min_eco_contacts

        Row 1 is the end to end eco from eco.log
        The subsequent rows are the eco scores between each frame 
        that are compatible with eco
        '''

        # Should house all the systems where eco / ANM has been attempted
        os.chdir(self.starting_dir) 

        # Access all the system directories, specifically step5
        for sub_dir in os.listdir(self.starting_dir):
            os.chdir(f'{sub_dir}/step5')
            with open('pathway.pdb') as file:
                lines = file.readlines()
                # Models are saved in pathway.pdb
                model_count = 0
                sub_file = None
                for line in lines:
                    if line.startswith('MODEL'):
                        if sub_file:
                            sub_file.close()
                        model_count += 1
                        sub_file = open(f'pathway_{model_count}.pdb', 'w')
                    # I don't want to write the MODEL and ENDMDL lines to the new file
                    if sub_file and not (line.startswith('MODEL') or line.startswith('ENDMDL')):
                        sub_file.write(line)
                if sub_file:
                    sub_file.close()
            print(f'Finished parsing pathway.pdb for system: {sub_dir}')
            # Start the pandas datframe with the end to end eco from eco.log
            df = pd.DataFrame(columns=['name', 'total_eco','min_eco_contacts']).to_csv(f'anm_eco_{sub_dir}.csv')
            with open('../eco.log') as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith('total ECO for this pathway:'):
                        total_eco = line.split(':')[-1].strip()
                    if line.startswith('the minECO pathway (by contact):'):
                        min_eco_contacts = line.split('):')[-1].strip()

                df = pd.DataFrame(data={'name': [sub_dir], 'total_eco': [total_eco], 'min_eco_contacts': [min_eco_contacts]}).to_csv(f'anm_eco_{sub_dir}.csv', mode='a', header=False)
            print(f'Created anm_eco_{sub_dir}.csv')
            print(f'Beginning eco calculations for system: {sub_dir}')
            # name, total_eco, min_eco_clusters, min_eco_contacts for headers
            # I need to try to run eco, and if it fails, I need to skip to the next pdb
            files = sorted([f for f in os.listdir() if f.startswith('pathway_') and f.endswith('.pdb')],
                        key=lambda x: int(x.split('_')[1].split('.')[0]))  # Sort numerically
            
            i = 0  # Start from the first file
            while i < len(files) - 1:  # Ensure there's a subsequent file to compare
                for j in range(i + 1, len(files)):  # Iterate over subsequent files
                    try:
                        total_eco, min_eco_contacts = self.run_eco(files[i], files[j], 0)

                        if total_eco != 0:
                            # Log the necessary data
                            df = pd.DataFrame(data={'name': [f'{files[i]}_{files[j]}'], 'total_eco': [total_eco], 
                                                    'min_eco_contacts': [min_eco_contacts]})
                            df.to_csv(f'anm_eco_{sub_dir}.csv', mode='a', header=False)

                            # Start again from the new index
                            i = j  
                            break  # Restart from new `i`
                    except Exception as e:
                        print(f'Error processing {files[i]} and {files[j]}: {e}')
                        continue

                else:
                    # If no total_eco != 0 was found, move to the next file
                    i += 1

            # Iterate through the models and run the eco calculation for each subsequent pair
            # Remove the .pdbs after the calculations are done
            for file in os.listdir():
                if file.startswith('pathway_') or file.endswith('.log'):
                    os.remove(file)
            os.chdir('../../')
        

test = PerFrameANMECO('/scratch/jacoba1/ECO-ANM-1.0.0/successful_eco')
test.iterative_eco()
