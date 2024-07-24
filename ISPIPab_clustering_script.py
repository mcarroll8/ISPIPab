#import statements
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math 
from pathlib import Path
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as shc
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans
import os, sys
from numpy import mean
import subprocess
import pandas as pd
from pathlib import Path
import sys


#localize the components of ISPIPab input and output
output_input = "/Users/moshecarroll/Downloads/Manuscript_Results_Finalized/ISPIPab_output_files"
#bin frame is a file outputted by ISPIPab, containing the per-residue epitope probability scores for each method, including XGBoost
#input below the path to the ISPIP bin frame results
df = pd.read_csv(f"{output_input}/bin_frame.csv")
#cutoff CSV (required input to ISPIPab) contains: PDB_ID,Number of Surface Residues,Dynamic Threshold Value,Number of Annotated Residues
#insert below the path to the ISPIP cutoff file
cutoff_path = "/Users/moshecarroll/Downloads/Antigen_Research/unbound_cutoff_2_23.csv"
#insert path to clustering outfile below
clustering_outfile = "/Users/moshecarroll/Downloads/clustering"
#insert path to pdb files below
pdb_folder = Path(f"/Users/moshecarroll/Downloads/Antigen_Research/unbound_pdbs")
predictor = "INSERT NAME OF PREDICTOR OF INTEREST"

predictors =[predictor]
#create directory to host clustering results
os.mkdir(f"{clustering_outfile}/{predictor}")
cutoff_csv = pd.read_csv(cutoff_path)

#function will use bin_frame to return a list of predicted residues for the method of interest for each antigen
def pml_predicted(protein,cutoff_csv,df,item_predictor):
    print(protein)
    frame = df[df["protein"] == protein]
    cutoff_row = cutoff_csv[cutoff_csv["Protein"] == protein]
    threshhold = cutoff_row["cutoff res"].values[0]
    predictedframesort = frame.sort_values(by=[item_predictor], inplace =False, ascending=False)
    thresholdframe = predictedframesort.head(threshhold) 
    predicted_res = thresholdframe.residue.values.tolist()
    predicted_res = [str(i) for i in predicted_res]
    pred_res = [i.split("_")[0] for i in predicted_res]
    pred_res_list = "\n".join(pred_res)
    #print(pred_res_list, protein)
    return pred_res_list
df["protein"] = [x.split('_')[1] for x in df.residue]
proteins = df["protein"].unique()
#below, the function will run, and create an inner folder to host clustering data for each protein as well as create a txt file listing the predicted residues
for item_predictor in predictors:
    for protein in proteins:

        pred_res_list= pml_predicted(protein,cutoff_csv,df,item_predictor)
        os.mkdir(f"{clustering_outfile}/{predictor}/{protein}")
        os.mkdir(f"{clustering_outfile}/{predictor}/{protein}/predicted_residues/")

        with open(f"{clustering_outfile}/{predictor}/{protein}/{protein}_unbound_predicted_residues.txt", 'a') as f:
            f.write(f"{pred_res_list}")

    print("pred_list_complete")


#function below defines spitting of the PDB file to capture each component of lines beginning with ATOM - ex: x,y,z coordinates of each atom in each residue of interest
def split_pdb_line(line):
    pdb_parts = [line[:6], line[6:11],
                line[12:16], line[16:20], line[21],
                line[22:28], line[30:38],
                line[38:46], line[46:54],
                line[54:60], line[60:66]]
    pdb_parts = list(map(lambda x: x.strip(), pdb_parts))
    return pdb_parts


pred_annotated_folder = Path(f"{clustering_outfile}/{predictor}")
for item in pred_annotated_folder.iterdir():
    for item_1 in item.iterdir():
        if "_predicted_residues" in item_1.name:
            print(item_1.name)
            for file in pdb_folder.iterdir():
                #print(file.name[12:16])
                #OPEN PDB FILE AND CONFIRM PDB IDs
                if file.name[12:16] == item.name:
                    with open(file) as pdb_file:
                        #DEPENDING ON PDB NAMING, CHANGE LINES 89,92-93 TO MATCH PROTEIN NAME TO PDB FILE NAME
                        protein_name_complete = file.name[12:25]
                        protein_bound_name = file.name[12:16]
                        #EXTRACT PDB DATA
                        for line in pdb_file:
                            if line.startswith("ATOM"):
                                line_data = split_pdb_line(line)
                                res_num = line_data[5]
                                res_name = line_data[3]
                                atom_name = line_data[2]
                                xcoordinated = line_data[6]
                                ycoordinated = line_data[7]
                                zcoordinated = line_data[8]
                                ann_protein_name = item_1.name[:4]
                                if ann_protein_name == protein_bound_name:
                                    with open(item_1) as annotated_file:
                                        for line1 in annotated_file:
                                            annotated_res_num = line1.strip()
                                            #print(annotated_res_num, res_num, protein)
                                            if (str(annotated_res_num) == str(res_num)): #change between pred and annotated
                                                #CREATES NEW FILE FOR EACH PREDICTED RESIDUE WITH PDB DATA DESCRIBED BELOW, INCLUDING COORDINATES FOR EACH ATOM THAT COMPRISES RESIDUE
                                                with open(f"{clustering_outfile}/{predictor}/{item.name}/predicted_residues/{res_num}.txt", 'a') as f:
                                                    f.write(f"{atom_name},{res_name},{res_num},{xcoordinated},{ycoordinated},{zcoordinated}\n")


#DETERMINES THE AVERAGE X, Y, AND Z COORDINATES FOR RESIDUE OF INTEREST BASED ON ATOM DATA (GEOMETRIC CENTERS)
folder_infile_1 = Path(f"{clustering_outfile}/{predictor}")
for item in folder_infile_1.iterdir():
    for item_1 in item.iterdir():
        if "predicted_residues" == item_1.name:
            for file in item_1.iterdir():
                list_x = []
                list_y = []
                list_z = []
                with open (file) as infile:
                    for line in infile:
                        xcoordinate = line.strip().split(",")[3]
                        ycoordinate = line.strip().split(",")[4]
                        zcoordinate = line.strip().split(",")[5]
                        list_x.append(float(xcoordinate))
                        list_y.append(float(ycoordinate))
                        list_z.append(float(zcoordinate))
            
                x_mean = round(mean(list_x),3)
                y_mean = round(mean(list_y),3)
                z_mean = round(mean(list_z),3)
                residue = file.name.split(".")[0]
                #CREATES NEW FOLDER WITH GEOMETRIC CENTER DATA
                with open (f"{item}/{item.name}_unbound_averaged_results.txt", "a") as outfile:
                    outfile.write(f"{residue}_Residue_Predicted,{x_mean},{y_mean},{z_mean}\n")



#code supports up to five clusters in organizing data, though can be expanded below. our results show max of three clusters in formation of epitopes based on hierarchical clustering
os.mkdir(f"{clustering_outfile}/{predictor}/hierarchical_cluster_1")
os.mkdir(f"{clustering_outfile}/{predictor}/hierarchical_cluster_2")
os.mkdir(f"{clustering_outfile}/{predictor}/hierarchical_cluster_3")
os.mkdir(f"{clustering_outfile}/{predictor}/hierarchical_cluster_4")
os.mkdir(f"{clustering_outfile}/{predictor}/hierarchical_cluster_5")







def remove(string):
    # get rid of the spaces in the list and turn it into a string. returns the list
    ones = ''
    for i in string:
            if i != " ":
                ones += str(i)
    #print(ones)
    return ones

def cluster():
    #loading the dataset
    print(file_2.name)
    dataset = pd.read_csv(file_2, header=None)
    
    data = dataset.iloc[:, 1:4].values
    #print(data)
    Z = linkage(data, method='ward')
    
    # fig = plt.figure(figsize=(5, 5))
    dend = shc.dendrogram(shc.linkage(data, method='ward'))
    print(dend['color_list'])
    unique_colors=set(dend['color_list'])
    
    optimal_num = len(unique_colors) -1    
    # dn = dendrogram(Z,
    #             orientation='top',
    #             distance_sort='descending',
    #             show_leaf_counts=True)

    # cluster!
    cluster = AgglomerativeClustering(n_clusters=optimal_num, metric='euclidean', linkage='ward', compute_full_tree=True, distance_threshold=None)
    
    cluster.fit_predict(data)    
    oneli = cluster.labels_
    if len(oneli) != len(data):
        print(" the len(oneli) is not equal to len(data)")
    print(oneli)
    print(data)
    # this is a scatter plot of the data
    # plt.figure(figsize=(10, 7))
    # plt.scatter(data[:,0],data[:,1], c=cluster.labels_, cmap='rainbow')
    # plt.show()
    
    #displaying of dendrogram
    plt.figure(figsize=(10, 7))
    plt.title(file_2.name)
    dend = shc.dendrogram(shc.linkage(data, method='ward'))
    plt.show()

    
    cluster1, cluster2, cluster3, cluster4, cluster5, data = getClusterList(oneli)
    averageXYZ(cluster1, cluster2, cluster3, cluster4, cluster5, data)
    

def getClusterList(ones):

    # ones is a list that corolates position to cluster 
    # ones can == [1 1 0 0 0 0 0 0 0 0 1 1 1 1 1].... IT GOES IN ORDER OF THE TXT FILE
    dataset = pd.read_csv(file_2, header=None)
    data = dataset.iloc[:, :].values
    clu1 = []
    clu2 = []
    clu3 = []
    clu4 = []
    clu5 = []
    ones = remove(ones)

    for row in range(len(data)):
        resNum = data[row][0]
        cluster = ones[row]
        # if the cluster is equal to zero
        if cluster == '0':
            # add up the sum of the points 
            clu1 += [resNum]
            
        elif cluster == '1':
            clu2 += [resNum]
        elif cluster == '2':
            clu3 += [resNum]
        elif cluster == '3':
            clu4 += [resNum]
        elif cluster == '4':
            clu5 += [resNum]
        else:
            pass
    
    print("cluster 1 is", clu1)
    print("cluster 2 is",clu2)
    print("cluster 3 is", clu3)
    return clu1, clu2, clu3, clu4, clu5, data

def averageXYZ(cluster1, cluster2, cluster3, cluster4, cluster5, data):
    #organizing cluster data to write to file

    string_clu1 = "\n".join(cluster1)
    string_clu3 = "\n".join(cluster3)
    string_clu4 = "\n".join(cluster4)
    string_clu5 = "\n".join(cluster5)

    clu3_length = len(cluster3)
    clu4_length = len(cluster4)
    clu5_length = len(cluster5)
    clu1_length = len(cluster1)
    clu2_length = len(cluster2)
    string_clu2 = "\n".join(cluster2)

    #writing files with cluster data
    protein_name = file.name[:4]
    if protein_name == protein:
        with open (cutoff_path) as infile_4:
            for cluster in infile_4:
                if cluster.strip().split(",")[0] == protein_name:
                    dynamic_cutoff = cluster.strip().split(",")[2]
                    with open (f"{clustering_outfile}/{predictor}/hierarchical_cluster_1/{protein_name}_{predictor}.txt", "a") as outfile_1:
                        outfile_1.write(f"{string_clu1}")
                    with open (f"{clustering_outfile}/{predictor}/hierarchical_cluster_2/{protein_name}_{predictor}.txt", "a") as outfile_2:
                        outfile_2.write(f"{string_clu2}")
                    with open (f"{clustering_outfile}/{predictor}/hierarchical_cluster_3/{protein_name}_{predictor}.txt", "a") as outfile_3:
                        outfile_3.write(f"{string_clu3}")
                    with open (f"{clustering_outfile}/{predictor}/hierarchical_cluster_4/{protein_name}_{predictor}.txt", "a") as outfile_4:
                        outfile_4.write(f"{string_clu4}")
                    with open (f"{clustering_outfile}/{predictor}/hierarchical_cluster_5/{protein_name}_{predictor}.txt", "a") as outfile_5:
                        outfile_5.write(f"{string_clu5}")
                    #writing grand file, containing summary data of clustering for each protein within method of interest
                    with open(f"{clustering_outfile}/{predictor}/{predictor}_hierarchical_clustering_results.csv", "a") as outfile:
                        outfile.write(f"{protein},{dynamic_cutoff},{clu1_length},{clu2_length},{clu3_length},{clu4_length},{clu5_length}\n")
#run clustering
def main():
    cluster()






folder = Path(f"{clustering_outfile}/{predictor}")    
for file in folder.iterdir():
    for file_2 in file.iterdir():
        if "_unbound_averaged_results" in file_2.name:
            protein = file.name[:4]
            main()   





#CLUSTER DATA WILL NOW BE AVAILABLE WITHIN PREDICTOR FOLDER FOR EACH ANTIGEN OF INTEREST
#F1-SCORES/MCC OF ANTIGEN CLUSTERS MAY NOW BE DETERMINED
