import pandas as pd 
import mdtraj as md
import numpy as np
import glob
import matplotlib.pylab as plt 
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pytraj as pt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.cm as cm
from scipy.spatial import distance


__author__="Houcemeddine Othman"
__email__="houcemoo@gmail.com"

def RMSD(trajectory, reference, adjust_time = 0.01):
    """
    Calculate RMSD
    """
    trajectory.center_coordinates()
    backbone = trajectory.top.select('not (resid 234 to 243 or resid 254 to 262 )  and backbone ') 
    backbone_ref = reference.top.select('not (resid 234 to 243 or resid 254 to 262 )  and backbone ') 
    #not (resid 234 to 243 or resid 254 to 262 )  and name CA
    rmsd = md.rmsd(trajectory, reference, 0, atom_indices=backbone, ref_atom_indices=backbone_ref )
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
    
def Rmsf(trajectory, topology, reference,  offset = 0, frame=0 ):
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
    return pd.DataFrame( {"resid":resid, "rmsf" :md.rmsf(t, reference, atom_indices=ca_selection)} )

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
    

def parameterPlt(pc_list, col_number=3, size_rows=10, mode="superpose"):
    """
    Calculates dimensions to plot PCs as a sqaure plot 
    """
    if mode == "superpose" : 
        cell_number = len(pc_list)-1
    elif  mode == "separate":
        cell_number = len(pc_list)    
    row_number = cell_number//col_number + cell_number%col_number
    size_cols = (size_rows * col_number)/row_number
    return ((size_cols,size_rows ) , (col_number,row_number  ))
    
        
def screePlot(list_pcs, labels, output):
    fig, ax = plt.subplots()
    for trajectory, label in zip(list_pcs, labels): 
        eigenvalues = trajectory[1][0]/np.sum( trajectory[1][0])
        eigenvalue_indexes = np.arange(1, len(eigenvalues)+1)
        ax.plot(eigenvalue_indexes, eigenvalues, "--o", label=label)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='y', which='minor', bottom=False)
        #ax.grid(color='gray', linestyle='--', linewidth=0.5,alpha=0.8, which="both")
        plt.ylabel("% of explained variances")
        plt.xlabel("Component index")
        plt.legend()
        plt.savefig(output)

def PcaPytraj(traj_list, references, sector = [0,-1] ): 
    PC_list = []
    for traj, reference in zip(traj_list, references): 
        PC_list.append( pt.pca(traj[sector[0]:sector[1]], ref=reference, mask='@CA,N,C', n_vecs=10)) 
    return PC_list
   
def plotPcs(list_pcs, row_number, col_number, index_of_wt = 0, pc1 = 1, pc2=2 ): 
    pc_list=list_pcs.copy()
    color_palette = ["#4bb1b4", "#B44E4B"]
    projection_wt = pc_list[index_of_wt][0]
    del pc_list[index_of_wt]
    for index, pcs  in enumerate(pc_list): 
        plt.subplot(row_number, col_number, index+1)
        projection = pcs[0]
        #plot PC1 vs PC2
        plt.scatter(projection_wt[pc1-1], projection_wt[pc2-1],alpha=0.5, color= color_palette[0])
        ymin, ymax = plt.gca().get_ylim()
        xmin, xmax = plt.gca().get_xlim()
        plt.scatter(projection[pc1-1], projection[pc2-1], alpha=0.5, color= color_palette[-1] )
        #plt.ylim([-50, 40])
        #plt.vlines(0,-70, 60, linestyles="dashed", color="k", alpha=0.5 )
        #plt.hlines(0,-70, 70, linestyles="dashed", color="k", alpha=0.5 )
        #plt.xlim([-40, 49])
        plt.tight_layout()
        xlabel="PC"+str(pc1)
        ylabel = "PC"+str(pc2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    
def Rmsip(eigs1, eigs2):
    """
    Calculates RMSIP of dimension_D
    """
    sum = 0
    dimension_D=len(eigs1)
    for eigenvector1 in eigs1 : 
        second_sum =  np.sum( [ np.dot(eigenvector1, eigenvector2 )**2 for eigenvector2 in eigs2 ] )
        sum = sum+second_sum
    return  np.round(np.sqrt(sum/dimension_D), 3)


def FelRange(pcs_for_var, pc_index):
    pcs = pcs_for_var[0] 
    pc = pcs[pc_index -1 ]
    return np.min(pc), np.max(pc)

def MinMaxPCsXY(list_pc, pcx_idx , pcy_idx): 
    """
    Returns max and mi values for PC values of X and 
        Y axis
    """
    min_pcx = 0
    min_pcy = 0
    max_pcx= 0 
    max_pcy = 0
    for trajectory_pc in list_pc: 
        pcs = trajectory_pc[0]
        pcx = pcs[pcx_idx-1]
        pcy = pcs[pcy_idx-1]
        if pcx.min() < min_pcx: 
            min_pcx = pcx.min()
        if pcy.min() < min_pcy:
            min_pcy = pcy.min()
        if pcx.max() > max_pcx:
            max_pcx = pcx.max()
        if pcy.max() > max_pcy:
            max_pcy = pcy.max()
    adjustment_factor = 10
    return np.array([[min_pcx-adjustment_factor, max_pcx+adjustment_factor], [min_pcy-adjustment_factor,  max_pcy+adjustment_factor]])

def ReturnConfAtMin(data_pc, pc_x, pc_y, minimum_at_bin ): 
    """
    returns the coordinates at the subspace corresponding
    to the closest conformation to the global minimum, calculated from binning 
    the free energy landscape
    """
    pcs = data_pc[0]
    pc1 = pcs[pc_x-1]
    pc2 = pcs[pc_y-1]
    minimum = np.array(minimum_at_bin)
    pc1_pc2 = np.column_stack((pc1, pc2))
    #distance.euclidean(pc1_pc2, minimum)
    distances = [ distance.euclidean(vec, minimum) for vec in pc1_pc2]
    min_distance = np.min(distances)
    index = int(np.where(distances == min_distance)[0])
    print("Conformation number {} is the closes to the global minimum".format(index+1))
    return pc1_pc2[index] # coordinates in the subspaces pc_x and pc_y

def FindMin(data_pc, pc_x, pc_y, bins, data_range ): 
    """
    returns the projection of the coordinates of the global minimum 
    from PCA data. 'data_pc', is the PCA from pytraj of a trajectory
    returns: 
        0: projection of the global minimum -> x axis 
        1: projection of the global minimum -> y axis
        2: range of pc_x values (min, max)
        3: range of pc_y values (min, max)
        4: Energy calculated from density data (type array ), with pc_x(raws) and pc_y(columns) as reaction coordinates
    """
    pcs = data_pc[0]
    pc1=pcs[pc_x -1 ]
    pc2=pcs[pc_y -1 ]   
    density=np.histogram2d(pc1, pc2, bins=bins, range=data_range) 
    matrix = density[0]
    flat_array = density[0].flatten()
    unique_values = np.unique(  np.sort(flat_array) )    # will convert all 0 values to the lowest density value 
    min_density = unique_values[1]
    matrix[ matrix == 0 ] = min_density
    energy = -0.001985875*298.15*np.log( matrix/matrix.sum() )
    min_energy =  energy.min()
    indexes_min = np.where(energy == min_energy)
    pc_projection_1 = density[1]   # get the pc values array from the second element in density 
    pc_projection_2 = density[2] 
    pc1_range = FelRange(data_pc, pc_x )
    pc2_range = FelRange(data_pc, pc_y )
    return  float(pc_projection_1[indexes_min[0]]), float(pc_projection_2[indexes_min[1]] ), pc1_range, pc2_range,  np.rot90(energy)
    
def getEnergyBoundaries(list_pc, pc_x, pc_y, bins):
    """ Reports the lowest and highest energies in the pc table """
    vmin = None
    vmax=None
    for data_pc in list_pc:
        pcs = data_pc[0]
        pc1=pcs[pc_x -1 ]
        pc2=pcs[pc_y -1 ]
        density=np.histogram2d(pc1, pc2, bins=bins)
        matrix = density[0]
        flat_array = density[0].flatten()
        unique_values = np.unique(  np.sort(flat_array) )    # will convert all 0 values to the lowest density value 
        min_density = unique_values[1]
        matrix[ matrix == 0 ] = min_density
        energy = -0.001985875*298.15*np.log( matrix/matrix.sum() )
        if vmin == None: 
            vmin = energy.min()
        elif vmax == None : 
            vmax = energy.max()
        elif energy.min() < vmin : 
            vmin = energy.min()
        elif energy.max() > vmax:
            vmax = energy.max()
    return vmin, vmax
            

def Fel(list_pc, col_number, row_number, bins, pc_x, pc_y): 
    """
    Plotting the FEL 
    """
    min_max_pcs= MinMaxPCsXY(list_pc, pc_x, pc_y )
    vmin, vmax = getEnergyBoundaries(list_pc,  pc_x, pc_y, bins )
    for index,data_pc in enumerate(list_pc): 
        plt.subplot(row_number, col_number, index+1)
        pcs = data_pc[0]     
        pc1_range = FelRange(data_pc, pc_x )
        pc2_range = FelRange(data_pc, pc_y )
        fel_data = FindMin(data_pc, pc_x, pc_y, bins, min_max_pcs )
        global_min_x = fel_data[0]
        global_min_y = fel_data[1] 
        minimum_at_data =  ReturnConfAtMin( data_pc , pc_x, pc_y, [global_min_x , global_min_y])
        max_energy = fel_data[4].max()
        fel_data[4][ fel_data[4] == max_energy ] = vmax   # make all the highst energies equals to vmax
        plt.imshow(fel_data[4], interpolation='bilinear', cmap=cm.jet, 
                   extent=[min_max_pcs[0][0] , min_max_pcs[0][1], min_max_pcs[1][0], min_max_pcs[1][1]] ,
                   aspect="auto", vmin=vmin, vmax=vmax)
        #plt.colorbar()
        plt.tight_layout()
        plt.scatter(pcs[0], pcs[1], color = "black", alpha =0.2, s=1)
        plt.scatter(minimum_at_data[0], minimum_at_data[1], s=100, color="red", marker="x")
        plt.scatter(pcs[0][0], pcs[1][1], color = "yellow", alpha =1, s=100, marker="x")
        #plt.contour( np.flip(fel_data[4], axis=0), colors='white', extent=[min_max_pcs[0][0] , min_max_pcs[0][1], min_max_pcs[1][0], min_max_pcs[1][1]] , alpha=0.5 )            
        xlabel="PC"+str(pc_x)
        ylabel = "PC"+str(pc_y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

class RmsfCalculation:
    def __init__(self, trajs, reference , offset, labels=[] , start=0, size_rows = 10, index_ref = 0): 
        self.traj_table = trajs
        self.reference = reference
        self.start= start
        self.offset = offset
        self.labels = labels
    
    def CalculateRmsf(self): 
        """
        Claculates the RMSF for all trajectrories in trajs 
        relative to a unique reference sructure. 
        Calculations are made for the CA atoms
        """
        self.rmsfs = []
        indexes_ref = self.reference.top.select("name CA")
        for sliced in self.traj_table :
            traj =sliced[100:200]
            indexes_traj = traj.top.select("name CA")
            if len(indexes_ref) != len(indexes_traj):
                raise Exception("Atoms in traj and reference do not match") 
            rmsf_traj = md.rmsf(traj, self.reference, atom_indices= indexes_traj, ref_atom_indices=indexes_ref)
            self.rmsfs.append(rmsf_traj)
    
    def dumpToCsv(self, output):
        """
        Dump calculated RMSFs go a csv file
        """
        resid = np.arange(len(self.rmsfs[0]))+1+self.offset
        rmsf_dic = { "resid":resid }
        if len(self.labels) != len(self.rmsfs): 
            raise Exception
        if self.labels != []: 
            for label, rmsf in zip(self.labels, self.rmsfs) : 
                rmsf_dic[label] = rmsf
        else: 
            raise Exception("You must provide the labels table")
        df = pd.DataFrame(rmsf_dic)
        df.to_csv(output, index=False)

class RmsfPlots: 
    def __init__(self, rmsf_csv_file, size_rows ): 
        self.rmsf_dataframe = pd.read_csv(rmsf_csv_file)
        self.size_rows= size_rows
                  
    def PlotRmsfVars(self, range_of_data, output ): 
        columns = self.rmsf_dataframe.columns[2:]
        working_df = self.rmsf_dataframe.loc[:, columns ]
        self.graph_dim = parameterPlt(columns, col_number = 1, size_rows=6, mode="separate")
        row_number = self.graph_dim[1][1]
        col_number = self.graph_dim[1][0]
        reference_rmsf = self.rmsf_dataframe["Ref"]        
        for index, col in enumerate(working_df) : 
            plt.subplot(row_number, col_number, index+1 )
            plt.plot(self.rmsf_dataframe.resid, reference_rmsf, color="#4bb1b4", marker="o", markersize=3, linestyle='dashed')
            plt.plot(self.rmsf_dataframe.resid, working_df[col], color="#B44E4B", marker="s" , markersize=3)
            plt.bar(self.rmsf_dataframe.resid, np.sqrt((working_df[col] - reference_rmsf)**2), color="gray"  )
            plt.xlim(range_of_data[0], range_of_data[1])
            plt.ylim(0,0.3)
            plt.tight_layout()
            plt.minorticks_on()
            #plt.vlines([160,170, 180,  190, 200, 210], 0, 0.5)
            plt.xlabel("#Residue")
            plt.ylabel("RMSF (nm)")
            plt.savefig(output)
            

