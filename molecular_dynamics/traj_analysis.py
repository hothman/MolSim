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
import matplotlib as mpl


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
    markers = ["o", "v", "s", "D", "^", "H", ">"]
    for trajectory, label, marker in zip(list_pcs, labels, markers): 
        eigenvalues = trajectory/np.sum( trajectory)
        eigenvalue_indexes = np.arange(1, len(eigenvalues)+1)
        ax.plot(eigenvalue_indexes, eigenvalues, "--", marker=marker, label=label)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='y', which='minor', bottom=False)
        #ax.grid(color='gray', linestyle='--', linewidth=0.5,alpha=0.8, which="both")
        plt.ylabel("% of explained variances")
        plt.xlabel("Component index")
        plt.legend()
        plt.savefig(output)

def cumVar(eigenvalues, to_which_eigenvalue=3): 
    eigenvalues_sum_subset = np.sum(eigenvalues[:to_which_eigenvalue])
    eigenvalues_sum = np.sum(eigenvalues)
    print(np.round( (eigenvalues_sum_subset/eigenvalues_sum)*100 ), "%" )

def PcaPytraj(traj_list, references, sector = [0,-1] ): 
    PC_list = []
    for traj, reference in zip(traj_list, references): 
        PC_list.append( pt.pca(traj[sector[0]:sector[1]], ref=reference, mask='@CA,N,C', n_vecs=10)) 
    return PC_list
   
def plotPcs(list_pcs, row_number, col_number, xmin , xmax, ymin,ymax, index_of_wt = 0, pc1 = 1, pc2=2, point_size=5  ): 
    pc_list=list_pcs.copy()
    color_palette = ["#4bb1b4", "#B44E4B"]
    # generate color map to ditinguish the time evolution in PCA plots
    colormap_var = [(248, 239, 239), (180, 78, 75)]  # dark red to light red 
    colormap_wt = [(219, 239, 240),(45, 106, 108)]  # dark cyan to light cyan
    cmap_wt = create_colormap(colormap_wt, bit=True)
    cmap_var  = create_colormap(colormap_var, bit=True)
    projection_wt = pc_list[index_of_wt]
    del pc_list[index_of_wt]
    for index, pcs  in enumerate(pc_list): 
        plt.subplot(row_number, col_number, index+1)
        projection = pcs
        #plot PC1 vs PC2     
        plt.axvline(x=0, linestyle='--', color='gray', alpha = 0.8, lw=1)
        plt.axhline(y=0, linestyle='--', color='gray', alpha = 0.8, lw=1)
        steps = np.arange(0,len(projection[pc1-1]))
        plt.scatter(projection[pc1-1], projection[pc2-1], c=steps, cmap=cmap_var, s=point_size )
        steps = np.arange(0,len(projection_wt[pc1-1]))  # generate array of time steps 0,1,2, ..... 
        plt.scatter(projection_wt[pc1-1], projection_wt[pc2-1], c=steps, cmap=cmap_wt, s=point_size)
        #plt.ylim([-50, 40])
        #plt.vlines(0,-70, 60, linestyles="dashed", color="k", alpha=0.5 )
        #plt.hlines(0,-70, 70, linestyles="dashed", color="k", alpha=0.5 )
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
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
    adjustment_factor = 32
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
    pcs = data_pc
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
        pcs = data_pc
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
    print(vmin, vmax)
    for index,data_pc in enumerate(list_pc): 
        plt.subplot(row_number, col_number, index+1)
        pcs = data_pc  
        pc1_range = FelRange(data_pc, pc_x )
        pc2_range = FelRange(data_pc, pc_y )
        fel_data = FindMin(data_pc, pc_x, pc_y, bins, min_max_pcs )
        global_min_x = fel_data[0]
        global_min_y = fel_data[1] 
        #minimum_at_data =  ReturnConfAtMin( data_pc , pc_x, pc_y, [global_min_x , global_min_y])
        max_energy = fel_data[4].max()
        fel_data[4][ fel_data[4] == max_energy ] = vmax   # make all the highst energies equals to vmax
        plt.imshow(fel_data[4], interpolation='bilinear', cmap=cm.Spectral_r, 
                   extent=[min_max_pcs[0][0] , min_max_pcs[0][1], min_max_pcs[1][0], min_max_pcs[1][1]] ,
                   aspect="auto", vmin=vmin, vmax=vmax)
        #plt.colorbar()
        plt.tight_layout()
        #plt.scatter(pcs[0], pcs[1], color = "black", alpha =0.05, s=0.001)
        #plt.scatter(minimum_at_data[0], minimum_at_data[1], s=100, color="red", marker="x")
        #plt.scatter(pcs[0][0], pcs[1][1], color = "yellow", alpha =1, s=100, marker="x")
        plt.contour( np.flip(fel_data[4], axis=0), colors='white', extent=[min_max_pcs[0][0] , min_max_pcs[0][1], min_max_pcs[1][0], min_max_pcs[1][1]] , alpha=0.5 )            
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
            plt.plot(self.rmsf_dataframe.resid, reference_rmsf, color="#4bb1b4", marker="o", lw=1, markersize=2, linestyle='dashed')
            plt.plot(self.rmsf_dataframe.resid, working_df[col], color="#B44E4B", marker="s", lw=1 , markersize=2)
            plt.bar(self.rmsf_dataframe.resid, np.sqrt((working_df[col] - reference_rmsf)**2), color="gray"  )
            plt.xlim(range_of_data[0], range_of_data[1])
            plt.ylim(0,0.3)
            plt.tight_layout()
            plt.minorticks_on()
            #plt.vlines([160,170, 180,  190, 200, 210], 0, 0.5)
            plt.xlabel("#Residue")
            plt.ylabel("RMSF (nm)")
            plt.savefig(output)
            
def plotlowerTriangleDccm(matrix, output, cmap, norm): 
    """
    Plots only the lower triangple of the dccm (or any other numpy square matrix)
    """
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    upper_matrix = np.tril(matrix)
    plt.imshow(upper_matrix, cmap=cmap, norm=norm, vmin=-1, vmax=1)
    #plt.colorbar()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.contour(upper_matrix , linewidths=0.2, colors="k", norm=norm)
    ax.set_xlabel("#Residue")
    ax.set_ylabel("#Residue")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='minor', width=1, length=3)
    fig.savefig(output, ddpi=400)
    plt.close(fig)

def cmapDccm(output):
    """
    Generates a colorbar and corresponding norm with 
    specific colors and range
    """
    cmap = mpl.colors.ListedColormap(['#2d6a6c', '#4bb1b4', "#93d0d2", '#dbeff0', 'white',
                                  '#f8efef', "#ddafad", "#c26e6c", "#B44E4B"])
    cmap.set_over('#B44E4B')
    cmap.set_under('#2d6a6c')
    bounds = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(figsize=(5, 0.5))
    fig.subplots_adjust(bottom=0.5)
    fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=ax,
        ticks=bounds,
        spacing='uniform',
        orientation='horizontal',
        label='',
    )
    fig.savefig(output)
    plt.close(fig)
    return cmap, norm

def plotDencityDccm(dccm_var, dccm_wt, output):
    """
    Plot the density of cross-correlation values (only for one triangle)
    """
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    dccm_n, _ = dccm_var.shape
    iu=np.tril_indices(469)
    lower_triangle_elements_var =  dccm_var[iu]
    lower_triangle_elements_wt =  dccm_wt[iu]
    y_var, x_var =np.histogram(lower_triangle_elements_var, bins=50)
    y_wt, x_wt =np.histogram(lower_triangle_elements_wt, bins=50)
    y_var = np.append([0], y_var)
    y_wt = np.append([0], y_wt)
    ax = plt.subplot(111)
    ax.plot( x_var, y_var/y_var.sum(), '--', color="#B44E4B", linewidth=2)
    ax.plot( x_wt, y_wt/y_wt.sum(), color="#4bb1b4", linewidth=2)
    ax.set_xlabel('Cross correlation')
    ax.set_ylabel('Density')
    ax.legend(['Variant', 'WT'])
    fig.savefig(output, ddpi=400)
    plt.close(fig)

def scalingData(value, max_radius): 
    """
    returns a scaled value to draw cylinder 
    radius in VMD for cross correlation analysis
    """
    value = np.abs(value)
    radius = max_radius * value
    return radius

def generateTclDccm(dccm_matrix, output, cutoff=0.5): 
    """
    Gnerates a tcl scrpt for vmd 
    correlation network
    """
    triangle = np.triu(dccm_matrix)
    indexes = np.where((triangle >cutoff) | (triangle <-cutoff))
    with open(output, "w") as tclfile:
        for x, y in zip(indexes[0], indexes[1]) : 
            if x != y:
                if triangle[x][y] > 0: 
                    color='cyan2'
                else:
                    color='red3'
                radius = scalingData(triangle[x][y], 0.2)
                commands =""" 
draw color  {3}
set atom1 [atomselect top "resid {0} and name CA"]
set atom2 [atomselect top "resid {1} and name CA"]
set coord1 [ $atom1 get {2} ]
set coord2 [ $atom2 get {2} ]
lassign $coord1 vector1
lassign $coord2 vector2
draw cylinder $vector1 $vector2 radius {4}
                        """.format(x, y, "{x y z}", color, radius)
                tclfile.write(commands)

#####################################################
# create_colormap() has been created by Chris Slocum
# https://github.com/CSlocumWX/custom_colormap
#####################################################

def create_colormap(colors, position=None, bit=False, reverse=False, name='custom_colormap'):
    """
    returns a linear custom colormap

    Parameters
    ----------
    colors : array-like
        contain RGB values. The RGB values may either be in 8-bit [0 to 255]
        or arithmetic [0 to 1] (default).
        Arrange your tuples so that the first color is the lowest value for the
        colorbar and the last is the highest.
    position : array like
        contains values from 0 to 1 to dictate the location of each color.
    bit : Boolean
        8-bit [0 to 255] (in which bit must be set to
        True when called) or arithmetic [0 to 1] (default)
    reverse : Boolean
        If you want to flip the scheme
    name : string
        name of the scheme if you plan to save it

    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
        cmap with equally spaced colors
    """
    from matplotlib.colors import LinearSegmentedColormap
    if not isinstance(colors, np.ndarray):
        colors = np.array(colors, dtype='f')
    if reverse:
        colors = colors[::-1]
    if position is not None and not isinstance(position, np.ndarray):
        position = np.array(position)
    elif position is None:
        position = np.linspace(0, 1, colors.shape[0])
    else:
        if position.size != colors.shape[0]:
            raise ValueError("position length must be the same as colors")
        elif not np.isclose(position[0], 0) and not np.isclose(position[-1], 1):
            raise ValueError("position must start with 0 and end with 1")
    if bit:
        colors[:] = [tuple(map(lambda x: x / 255., color)) for color in colors]
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    return LinearSegmentedColormap(name, cdict, 256)

class Projection: 
    def __init__(self, trajectory_path, topology, eigenvalues, eigenvectors): 
        """
        A class that calculates the projection of coordinates over a set of 
        protein modes.
        """
        self.pytraj_like_frames = pt.load(trajectory_path, topology) 
        self.trajectory_path = trajectory_path
        self.topology = topology
        self.eigenvalues = eigenvalues 
        self.eigenvectors = eigenvectors 
    
    def reformatTrajObject(self):
        """ Will generate a traj object of 3 frames if 
        the input contains less than 3 snapshots. This is 
        necessary in order for the projection to give 
        reliable results """
        if len(self.pytraj_like_frames ) < 3 : 
            self.processed_traj =  pt.iterload([self.trajectory_path, self.trajectory_path, self.trajectory_path], top=self.topology) 
        else: 
            self.processed_traj = self.pytraj_like_frames
    
    def getProjection(self, wild_card=':1-232,244-253,263-470@CA,N,C'):
        """
        Project the coodinates of snapshots in self.processed_traj in all the 
        subspaces described by the given eigenvectors/eigenvalues 
        """
        self.reformatTrajObject()
        data = pt.projection(self.processed_traj, wild_card, eigenvalues=self.eigenvalues, eigenvectors=self.eigenvectors, scalar_type='covar')
        return data
