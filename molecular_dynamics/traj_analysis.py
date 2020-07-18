import pandas as pd 
import mdtraj as md
import numpy as np
import glob
import matplotlib.pylab as plt 
from sklearn.decomposition import PCA

__author__="Houcemeddine Othman"
__email__="houcemoo@gmail.com"

def RMSD(trajectory, topology):
    """
    Calculate RMSD
    """
    t = md.load(trajectory, top=topology)
    t.center_coordinates()
    backbone = t.top.select('backbone')
    rmsd = md.rmsd(t, t, 0, atom_indices=backbone)
    time = np.arange(len(rmsd))
    return pd.DataFrame({"Time":time , "RMSD":rmsd })

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
    trajectory.center_coordinates()
    top = trajectory.topology
    ca_selection = trajectory.top.select('name CA')
    resid = range(1+offset, len(ca_selection)+1+offset) 
    return pd.DataFrame( {"resid":resid, "rmsf" :md.rmsf(trajectory, trajectory, frame=frame, atom_indices=ca_selection)} )

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
    # reduce the high order array
    coord_array = []
    for frame in t: 
        ca_coordinates = frame.xyz[:, ca_selection]
        frame_array = list(np.concatenate(ca_coordinates, axis=None)) 
        coord_array.append(frame_array)
    coord_array = np.array(coord_array)  
    # PCA calculation
    pca_traj = PCA(n_components= n_components)
    principalComponents = pca_traj.fit_transform(coord_array)
    return pd.DataFrame(principalComponents)    

    

