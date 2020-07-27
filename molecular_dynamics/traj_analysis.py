import pandas as pd 
import mdtraj as md
import numpy as np
import glob
import matplotlib.pylab as plt 
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pytraj as pt

__author__="Houcemeddine Othman"
__email__="houcemoo@gmail.com"

def RMSD(trajectory, adjust_time = 0.01):
    """
    Calculate RMSD
    """
    trajectory.center_coordinates()
    backbone = trajectory.top.select('backbone')
    rmsd = md.rmsd(trajectory, trajectory, 0, atom_indices=backbone)
    time = np.arange(len(rmsd))
    return pd.DataFrame({"Time":time*adjust_time+ adjust_time, "RMSD":rmsd })

def RmsdAvgMinMax(path_to_traj_folder, topology):
    """
    Generate min max vector for traj 
    Suppose that traj has files are NetCDF and 
    use extension ".nc"
    """
    RMSD_list = []
    for traj in glob.glob(path_to_traj_folder+"/*.nc"):
        rmsd = RMSD(traj,topology )
        RMSD_list.append(rmsd.RMSD)
    RMSD_list = np.array(RMSD_list)
    average = np.average(RMSD_list, axis=0)
    min_rmsd = np.min(RMSD_list, axis=0)
    max_rmsd = np.max(RMSD_list, axis=0)
    return pd.DataFrame({"time":rmsd.Time*0.01, "avg_rmsd": average, "min_rmsd": min_rmsd, "max_rmsd": max_rmsd})
    
def PlotRmsd(rmsdDf): 
    """
    Plot RMSD of trajectories simulated from different 
        seeds
    """
    plt.plot(rmsdDf.time, rmsdDf.avg_rmsd)
    plt.fill_between(rmsdDf.time, rmsdDf.max_rmsd, rmsdDf.min_rmsd, alpha= 0.5)
    
def Rmsf(trajectory, topology, frame=0, offset = 0):
    """
    skip: how many frames to skip 
    offset: use this option to adjust the resid according to the 
        reference sequence
    """
    t = trajectory[frame:]
    t.center_coordinates()
    top = t.topology
    ca_selection = t.top.select('name CA')
    resid = range(1+offset, len(ca_selection)+1+offset) 
    return pd.DataFrame( {"resid":resid, "rmsf" :md.rmsf(t, t, atom_indices=ca_selection)} )

def plotRmsf(Dataframe_rmsf):
    plt.plot(Dataframe_rmsf.resid, Dataframe_rmsf.rmsf)
    
def CombineTrajs(path_to_traj_folder, topology, skipframes=1000):
    traj_paths = glob.glob(path_to_traj_folder+"/*.nc")
    # add reference snapshot to traj at the first position 
    concat_traj = md.load(traj_paths[0], top=topology)[0]
    for traj in traj_paths:
        t = md.load(traj, top=topology)
        concat_traj=md.join([concat_traj, t[skipframes:]]  )
    return concat_traj


def PcaTraj(traj_object, n_components=2):
    """
    calculates the principle components of trajectory 
    based on CA atoms coordinates. 
        traj_object: MDTraj trajectory object
        n_components: number of PCs to calculate
    """
    ca_selection = traj_object.top.select('name CA')   # indexes for the CA atoms
    traj_superposed = traj_object.superpose(traj_object, frame=0, atom_indices=ca_selection, ref_atom_indices=ca_selection, parallel=True)
    number_of_coordinates = len(ca_selection)*3
    print("Dimension:       {0},{1}".format(number_of_coordinates, len(traj_superposed) ))
    # reduce the high order array
    coord_array = []
    for frame in traj_superposed: 
        ca_coordinates = frame.xyz[:, ca_selection]
        frame_array = list(np.concatenate(ca_coordinates, axis=None)) 
        coord_array.append(frame_array)
    coord_array = np.array(coord_array)  
    # PCA calculation
    pca_traj = PCA(n_components= n_components)
    principalComponents = pca_traj.fit_transform(coord_array)
    scree = pca_traj.explained_variance_ratio_
    return pd.DataFrame(principalComponents), scree    
    

def parameterPlt(pc_list, col_number=3, size_rows=10):
    """
    Calculates dimensions to plot PCs as a sqaure plot 
    """
    cell_number = len(pc_list)-1
    row_number = cell_number//col_number + cell_number%col_number
    size_cols = (size_rows * col_number)/row_number
    return ((size_cols,size_rows ) , (col_number,row_number  ))
    

def plotPcs(list_pcs, row_number, col_number, index_of_wt = 0, pc1 = 1, pc2=2 ): 
    pc_list=list_pcs.copy()
    projection_wt = pc_list[index_of_wt][0]
    del pc_list[index_of_wt]
    for index, pcs  in enumerate(pc_list): 
        plt.subplot(row_number, col_number, index+1)
        projection = pcs[0]
        #plot PC1 vs PC2
        plt.scatter(projection_wt[pc1-1], projection_wt[pc2-1], color= color_game[0])
        plt.scatter(projection[pc1-1], projection[pc2-1], alpha=0.5, color= color_game[1] )
        plt.tight_layout()
        xlabel="PC"+str(pc1)
        ylabel = "PC"+str(pc2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
