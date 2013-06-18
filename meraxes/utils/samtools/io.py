import numpy as np
import h5py as h5

def read_gals(fname, firstfile=None, lastfile=None, snapshot=None, 
                      props=None, **kwargs):
    
    """ Read in a Meraxes hdf5 output file.

    Reads in the default type of HDF5 file generated by the code.

    Args:
        fname: Full path to input hdf5 master file.
        snapshot: (optional) The snapshot to read in.
                  default = The last snapshot present in the output.
                            (usually z=0)
        props: (optional) A list of galaxy properties requested.
               default = All properties.
        verbose: (optional) Print extra info and status during read.
                 default = True
        sim_props: (optional) Output some simulation properties as well.
                   default = False

    Returns:
        A ndarray with the requested galaxies and properties.

        If sim_props==True then output a tuple of form 
        (galaxies, sim_props) where sim_props hold the following information
        as a dictionary:
        ( BoxSize, 
          MaxTreeFiles, 
          ObsHubble_h,
          Volume,
          Redshift )
    """

    verbose = kwargs.get('verbose', True)
    output_sim_props_flag = kwargs.get('output_sim_props', False)

    # Open the file for reading
    fin = h5.File(fname, 'r')

    # Set the snapshot correctly
    if snapshot==None: snapshot=-1
    MaxSnaps = fin['InputParams'].attrs['LastSnapShotNr'][0]+1
    if snapshot<0: snapshot+=MaxSnaps
    print "Reading snapshot %d" % snapshot

    # Select the group for the requested snapshot.
    snap_group = fin['Snap%03d'%(snapshot)]

    # Create a dataset large enough to hold all of the requested galaxies
    ngals = snap_group['Galaxies'].size
    if props!=None:
        gal_dtype = snap_group['Galaxies'].value[list(props)[:]][0].dtype
    else:
        gal_dtype = snap_group['Galaxies'].dtype

    G = np.empty(ngals, dtype=gal_dtype)
    if verbose:
        print "Allocated %.1f MB" % (G.itemsize*ngals/1024./1024.)

    # Loop through each of the requested groups and read in the galaxies
    if ngals>0:
        snap_group['Galaxies'].read_direct(G, dest_sel=np.s_[:ngals])

    # Print some checking statistics
    if verbose:
        print 'Read in %d galaxies.' % len(G)

    # Set some run properties
    if output_sim_props_flag:
        SimHubble_h  = fin['InputParams'].attrs['SimHubble_h'][0]
        ObsHubble_h  = fin['InputParams'].attrs['ObsHubble_h'][0]
        BoxSize      = fin['InputParams'].attrs['BoxSize'][0] / ObsHubble_h
        MaxTreeFiles = fin['InputParams'].attrs['FilesPerSnapshot'][0]
        VolumeFactor = fin['InputParams'].attrs['VolumeFactor'][0]
        Volume       = BoxSize**3.0 * (lastfile - firstfile + 1.) /\
                        MaxTreeFiles * VolumeFactor
        Redshift     = snap_group.attrs['Redshift']
        fin.close()
        return G, {'BoxSize':BoxSize, 
                   'MaxTreeFiles':MaxTreeFiles, 
                   'ObsHubble_h':ObsHubble_h, 
                   'SimHubble_h':SimHubble_h, 
                   'Volume':Volume,
                   'Redshift':Redshift}

    else:
        fin.close()
        return G


def read_input_params(fname, props=None):

    """ Read in the input parameters from a Meraxes hdf5 output file.

    Reads in the default type of HDF5 file generated by the code.

    Args:
        fname: Full path to input hdf5 master file.
        props: (optional) A list of run properties requested.
                default = All properties.
    
    Returns:
        A dict with the requested run properties.
    """

    # Initialise the output dictionary
    props_dict = {}

    # Open the file for reading
    fin = h5.File(fname, 'r')

    group = fin['InputParams']

    if props == None:
        props = group.attrs.keys()

        # Add some extra properties
        props_dict['SimHubble_h']  = group.attrs['SimHubble_h'][0]
        props_dict['ObsHubble_h']  = group.attrs['ObsHubble_h'][0]
        props_dict['BoxSize']      = group.attrs['BoxSize'][0] / props_dict['ObsHubble_h']
        props_dict['MaxTreeFiles'] = group.attrs['FilesPerSnapshot'][0]
        props_dict['VolumeFactor'] = group.attrs['VolumeFactor'][0]
        props_dict['Volume']       = props_dict['BoxSize']**3.0 * (group.attrs['LastFile'] - group.attrs['FirstFile'] + 1.) /\
                                        props_dict['MaxTreeFiles'] * props_dict['VolumeFactor']
    
    for p in props:
        try:
            props_dict[p] = group.attrs[p][0]
        except (KeyError):
            print "Property '%s' doesn't exist in the InputParams group." % p

    fin.close()

    return props_dict


def read_gitref(fname):
    """Read the git ref saved in the master file.
   
    Args:
        fname:  Full path to input hdf5 master file.

    Returns:
        gitref:     git ref of the model
    """

    fin = h5.File(fname, 'r')
    gitref = fin.attrs['GitRef'].copy()[0]
    fin.close()

    return modelname, gitref


def read_snaplist(fname):

    """ Read in the list of available snapshots from the Meraxes hdf5 file. 
    
    Args:
        fname: Full path to input hdf5 master file.
    
    Returns:
        snaps:      ndarray of snapshots
        redshifts:  ndarray of redshifts
        lt_times:   ndarray of light travel times (Gyr)
    """

    fin = h5.File(fname, 'r')

    zlist = []
    snaplist = []
    lt_times = []
    for snap in fin.keys():
        try:
            zlist.append(fin[snap].attrs['Redshift'][0])
            snaplist.append(int(snap[-3:]))
            lt_times.append(fin[snap].attrs['LTTime'][0])
        except KeyError:
            pass

    return np.array(snaplist, dtype=float), np.array(zlist, dtype=float)\
            , np.array(lt_times, dtype=float)

