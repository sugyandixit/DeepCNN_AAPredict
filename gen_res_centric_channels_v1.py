import os
import pickle
import _pickle as cpickle
import bz2
# import h5py
import itertools
import numpy as np
import Bio
import collections
import glob
import torch
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as g_o
import time
from memory_profile import profile

res_name_label_dict = {'HIS': 0,
                       'LYS': 1,
                       'ARG': 2,
                       'ASP': 3,
                       'GLU': 4,
                       'SER': 5,
                       'THR': 6,
                       'ASN': 7,
                       'GLN': 8,
                       'ALA': 9,
                       'VAL': 10,
                       'LEU': 11,
                       'ILE': 12,
                       'MET': 13,
                       'PHE': 14,
                       'TYR': 15,
                       'TRP': 16,
                       'PRO': 17,
                       'GLY': 18,
                       'CYS': 19}

van_der_waal_radius = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}

class PDBAtom:
    """
    class to store atom type, residue, chain ID, x,y,z coordinates, index, and value
    """

    def __init__(self, atom_id, atom_type, atom_index, res_name, res_num, chain_id, model_num,
                 x_coord, y_coord, z_coord, value):
        """
        class to store pdb info for individual atoms
        :param atom_id: element id
        :param atom_type: atom type in proteins 'CA', 'CZ', 'NZ' etc
        :param atom_index: atom index
        :param res_name: residue name
        :param res_num: residue number
        :param chain_id: chain id
        :param model_num: model number
        :param x_coord: x coordinates
        :param y_coord: y coordinates
        :param z_coord: z coordinates
        :param value: value of 1 when filled
        """

        self.atom_index = atom_index
        self.atom_id = atom_id
        self.atom_type = atom_type
        self.res_name = res_name
        self.res_num = res_num
        self.chain_id = chain_id
        self.model_num = model_num
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.value = value

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


def save_pickle_object(obj, fpath):

    with open(fpath, 'wb') as outfile:
        pickle.dump(obj, outfile)


def save_cpickle(obj, fpath):

    with bz2.BZ2File(fpath, 'w') as outfile:
        cpickle.dump(obj, outfile)



def collect_PDB_info(pdb_fpath):
    """
    read the pdb structure and store info using the PDBAtom class
    :param pdb_fpath: pdb filepath
    :return: pdb model dictionary containing PDBAtoms for all atoms (plus some additional data)
    """

    pdb_id = str(os.path.split(pdb_fpath)[1]).split('.')[0]
    ID_dict = collections.OrderedDict()

    pdb_model = dict()

    all_coords = []
    all_x = []
    all_y = []
    all_z = []
    all_atom_type = []
    all_atom_id = []
    all_res_name = []
    all_res_label = []
    PDB_entries = []

    parser = PDBParser()
    struct = parser.get_structure(str(os.path.split(pdb_fpath)[1]), pdb_fpath)

    residues = [res for res in struct[0].get_residues()]
    sequence = ''.join([Bio.PDB.Polypeptide.three_to_one(aa.get_resname()) for aa in residues])

    atom_index = 0
    for model in struct:
        if model._id > 0:
            break
        for chain in model:
            for residue in chain:
                if residue._id[0] == ' ':
                    res_name = residue.resname
                    for atom in residue:
                        if not res_name == 'GLY':
                            if 'CB' not in residue.child_dict.keys():
                                break
                        model_num = model._id
                        chain_id = chain._id
                        res_field = residue._id
                        res_num = res_field[1]
                        res_name = residue.resname
                        res_label = res_name_label_dict[res_name]
                        atom_name = atom.get_name()
                        atom_id = atom.element
                        coords = atom.coord
                        x_coord = float(coords[0])
                        y_coord = float(coords[1])
                        z_coord = float(coords[2])
                        all_x.append(x_coord)
                        all_y.append(y_coord)
                        all_z.append(z_coord)
                        all_coords.append([x_coord, y_coord, z_coord])
                        all_atom_id.append(atom_id)
                        all_atom_type.append(atom_name)
                        all_res_name.append(res_name)
                        all_res_label.append(res_label)

                        if chain_id not in ID_dict.keys():
                            ID_dict[chain_id] = [PDBAtom(atom_id=atom_id,
                                                         atom_type=atom_name,
                                                         atom_index=atom_index,
                                                         res_name=res_name,
                                                         res_num=res_num,
                                                         chain_id=chain_id,
                                                         model_num=model_num,
                                                         x_coord=x_coord,
                                                         y_coord=y_coord,
                                                         z_coord=z_coord,
                                                         value=1)]
                        else:
                            ID_dict[chain_id].append(PDBAtom(atom_id=atom_id,
                                                             atom_type=atom_name,
                                                             atom_index=atom_index,
                                                             res_name=res_name,
                                                             res_num=res_num,
                                                             chain_id=chain_id,
                                                             model_num=model_num,
                                                             x_coord=x_coord,
                                                             y_coord=y_coord,
                                                             z_coord=z_coord,
                                                             value=1))
                        PDB_entries.append(PDBAtom(atom_id=atom_id,
                                                   atom_type=atom_name,
                                                   atom_index=atom_index,
                                                   res_name=res_name,
                                                   res_num=res_num,
                                                   chain_id=chain_id,
                                                   model_num=model_num,
                                                   x_coord=x_coord,
                                                   y_coord=y_coord,
                                                   z_coord=z_coord,
                                                   value=1))
                        atom_index += 1

    pdb_model['pdb_id'] = pdb_id
    pdb_model['sequence'] = sequence
    pdb_model['sequence_length'] = len(sequence)
    pdb_model['ID_dict'] = ID_dict
    pdb_model['PDB_entries'] = PDB_entries
    pdb_model['all_atom_id'] = all_atom_id
    pdb_model['all_atom_type'] = all_atom_type
    pdb_model['all_res_name'] = all_res_name
    pdb_model['all_res_label'] = all_res_label
    pdb_model['all_coords'] = all_coords
    pdb_model['all_x'] = all_x
    pdb_model['all_y'] = all_y
    pdb_model['all_z'] = all_z
    return pdb_model


def get_position_dict(pdb_atoms, res_num):
    """
    get dict info on atom coordinates, chain id, res name, res num, and atom index based atom type
    :param pdb_atoms: list of pdb atoms
    :param res_num: residue number
    :return: position dictionary
    """
    position_dict = {}
    for pdb_atom in pdb_atoms:
        if pdb_atom.res_num == res_num:
            position_dict[pdb_atom.atom_type] = [(pdb_atom.x_coord, pdb_atom.y_coord, pdb_atom.z_coord),
                                                 pdb_atom.chain_id, pdb_atom.res_name, pdb_atom.res_num,
                                                 pdb_atom.atom_index]
    return position_dict


def center_and_transform(res_name, position_dict):
    """
    generate x, y, and z axis based on a reference point and a transformation matrix to transform coordinates so that
    ref point is centered
    :param res_name: reference residue name
    :param position_dict: reference residue position dictionary
    :return: reference point, transformation matrix
    """
    reference_pt = position_dict['CA'][0]
    axis_x = np.array(position_dict['N'][0]) - np.array(position_dict['CA'][0])
    pseudo_axis_y = np.array(position_dict['C'][0]) - np.array(position_dict['CA'][0])
    axis_z = np.cross(axis_x, pseudo_axis_y)
    if not res_name == 'GLY':
        direction = np.array(position_dict['CB'][0]) - np.array(position_dict['CA'][0])
        axis_z *= np.sign(direction.dot(axis_z))
    axis_y = np.cross(axis_z, axis_x)

    axis_x /= np.sqrt(sum(axis_x ** 2))
    axis_y /= np.sqrt(sum(axis_y ** 2))
    axis_z /= np.sqrt(sum(axis_z ** 2))

    transform = np.array([axis_x, axis_y, axis_z]).T
    return reference_pt, transform


def plot_3d_coords_comp_transform(all_coords_1, all_coords_2):
    """
    plot scatter in 3d to compare transformed coordiatnes
    :param all_coords_1: coordinate list 1
    :param all_coords_2: coordinate list 2
    :return:
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(all_coords_1[:, 0], all_coords_1[:, 1], all_coords_1[:, 2], color='red', alpha=0.5)
    ax.scatter(all_coords_2[:, 0], all_coords_2[:, 1], all_coords_2[:, 2], color='blue', alpha=0.5)
    plt.show()


def plot_Xsmooth(density_arr, box_len=21):
    """
    plotting density using plotly volume module
    :param density_arr: density array
    :param box_len: box length
    :return: figure
    """

    X, Y, Z = np.meshgrid(np.arange(box_len), np.arange(box_len), np.arange(box_len))

    values = density_arr

    fig = g_o.Figure(data=g_o.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        opacity=0.1,
        surface_count=25,
        colorscale='Magma'
    ))
    # fig.update_layout(scene_xaxis_showticklabels=False,
    #                   scene_yaxis_showticklabels=False,
    #                   scene_zaxis_showticklabels=False)

    return fig


def gen_grid_points(box_length, pixel_size):
    """
    generate grid points based on number of pixels and pixel size
    :param box_size: box_size
    :param pixel_size: pixel size in Angstroms
    :return: grid points, number of grids in an axis
    """

    grid_axis = []

    num_grids_in_half_axis = box_length / 2
    if not box_length % 2 == 0:
        num_grids_in_half_axis = (box_length - 1) / 2

    for num in range(1, int(num_grids_in_half_axis+1)):
        pt = (num * pixel_size)
        grid_axis.append(pt)
        grid_axis.append(-pt)

    grid_axis = [0] + grid_axis

    grid_axis = sorted(grid_axis)
    num_grids_axis = len(grid_axis)

    grid_points = np.array(list(itertools.product(grid_axis, grid_axis, grid_axis)), dtype=float)

    return grid_points, num_grids_axis


def gen_channel_coords(channel_list, protein_model_dict, center_transformed_pos, box_length, center_res_atoms_exclude_index=[],
                       use_atom_radius_pixel_size=False, pixel_size=1):

    """
    generate grid points and coords based on the channel list (atom id)
    :param channel_list: list of atom ids
    :param protein_model_dict: protein model dict from collect PDB info
    :param center_transformed_pos: center and transformed position
    :param box_length: length of box
    :param center_res_atoms_exclude_index: any residues to exclude from the center residue
    :return: grid_box_valid (boolean to say if a box has any atoms (True) or not), grid_point_channels (grid_points for
    each channel), coords_channels (coordinates of atoms in box for each channel), num_grids_axis (number of grid points
    in one axis of the box)
    """

    grid_box_valid = [True] * len(channel_list)
    grid_points_channels = [[] for _ in range(len(channel_list))]
    coords_channels = [[] for _ in range(len(channel_list))]

    if use_atom_radius_pixel_size:

        for num in range(len(channel_list)):
            channel = channel_list[num]

            grid_points, num_grids_axis = gen_grid_points(box_length=box_length,
                                                              pixel_size=van_der_waal_radius[channel])
            grid_points_channels[num].append(grid_points)

            min_coord, max_coord = np.min(grid_points), np.max(grid_points)

            include_coords_index = np.where((center_transformed_pos[:, 0] > min_coord) &
                                              (center_transformed_pos[:, 0] < max_coord) &
                                              (center_transformed_pos[:, 1] > min_coord) &
                                              (center_transformed_pos[:, 1] < max_coord) &
                                              (center_transformed_pos[:, 2] > min_coord) &
                                              (center_transformed_pos[:, 2] < max_coord))

            final_index = [int(x) for x in include_coords_index[0] if x not in center_res_atoms_exclude_index]
            box_atoms = [protein_model_dict['PDB_entries'][i] for i in final_index]
            box_channel_atoms = [x for x in box_atoms if x.atom_id == channel]
            box_channel_atoms_index = [x.atom_index for x in box_channel_atoms]
            pos_in_box = center_transformed_pos[box_channel_atoms_index]
            coords_channels[num].append(pos_in_box)
            if len(box_channel_atoms) > 0:
                grid_box_valid[num] = True
            else:
                grid_box_valid[num] = False

    else:
        grid_points, num_grids_axis = gen_grid_points(box_length=box_length,
                                                      pixel_size=pixel_size)
        min_coord, max_coord = np.min(grid_points), np.max(grid_points)
        include_coords_index = np.where((center_transformed_pos[:, 0] > min_coord) &
                                        (center_transformed_pos[:, 0] < max_coord) &
                                        (center_transformed_pos[:, 1] > min_coord) &
                                        (center_transformed_pos[:, 1] < max_coord) &
                                        (center_transformed_pos[:, 2] > min_coord) &
                                        (center_transformed_pos[:, 2] < max_coord))

        final_index = [int(x) for x in include_coords_index[0] if x not in center_res_atoms_exclude_index]
        box_atoms = [protein_model_dict['PDB_entries'][i] for i in final_index]

        for num in range(len(channel_list)):
            channel = channel_list[num]
            box_channel_atoms = [x for x in box_atoms if x.atom_id == channel]
            box_channel_atoms_index = [x.atom_index for x in box_channel_atoms]
            pos_in_box = center_transformed_pos[box_channel_atoms_index]
            coords_channels[num].append(pos_in_box)
            grid_points_channels[num].append(grid_points)
            if len(box_channel_atoms) > 0:
                grid_box_valid[num] = True
            else:
                grid_box_valid[num] = False


    return grid_box_valid, grid_points_channels, coords_channels, num_grids_axis


def gen_res_center_density_in_box(protein_model_dict, center_pos_list, channels_list=['N', 'C', 'O', 'S'], box_length=20):

    """
    generate residue center density for all channels. uses torch tensors. if cuda device is available, torch tensors
    are directed to cuda device.
    :param protein_model_dict: protein odel dictionary
    :param center_pos_list: center residue position list
    :param channels_list: channels list
    :param box_length: length of the box
    :return: channel_grid_voxel_density (grid voxel density for each channel)
    """

    center_pos, chain_ID, center_res_name, center_res_num,  = center_pos_list

    center_pos_dict = get_position_dict(protein_model_dict['ID_dict'][chain_ID], center_res_num)
    center_res_side_chain_atoms_index = []
    center_res_backbone_atoms_index = []
    for keys, values in center_pos_dict.items():
        if not keys == 'N' and not keys == 'CA' and not keys == 'C' and not keys == 'O':
            center_res_side_chain_atoms_index.append(values[4])
        else:
            center_res_backbone_atoms_index.append(values[4])

    ref_pos, transform_mat = center_and_transform(center_res_name, center_pos_dict)

    all_pos = np.array(protein_model_dict['all_coords'])
    transformed_pos = ((all_pos - ref_pos).dot(transform_mat))

    grid_box_valid, channel_grid_points, channel_coords, num_grids_axis = gen_channel_coords(
        channel_list=channels_list,
        protein_model_dict=protein_model_dict,
        center_transformed_pos=transformed_pos,
        box_length=box_length,
        center_res_atoms_exclude_index=center_res_side_chain_atoms_index,
        use_atom_radius_pixel_size=False,
        pixel_size=1)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    channel_grid_voxel_density = []

    for index, val_box in enumerate(grid_box_valid):

        if val_box == True:
            grid_points = torch.from_numpy(channel_grid_points[index][0]).to(device=device, dtype=torch.float)
            coords_in_box = torch.from_numpy(channel_coords[index][0]).to(device=device, dtype=torch.float)
            coord_grid_dist = torch.cdist(grid_points, coords_in_box)
            normal_dist = torch.distributions.normal.Normal(torch.tensor([0.0], device=device, dtype=torch.float), torch.tensor([1.0], device=device, dtype=torch.float))
            all_grid_density = torch.exp(normal_dist.log_prob(coord_grid_dist))
            grid_sum_density = torch.sum(all_grid_density, dim=1)
            grid_voxel_density = grid_sum_density.reshape((num_grids_axis, num_grids_axis, num_grids_axis))
        else:
            grid_voxel_density = torch.zeros((num_grids_axis, num_grids_axis, num_grids_axis), dtype=torch.float,
                                             device=device)


        # fig = plot_Xsmooth(grid_voxel_density)
        # fig.show()


        channel_grid_voxel_density.append(grid_voxel_density)

    channel_grid_voxel_density = torch.stack(channel_grid_voxel_density)

    return channel_grid_voxel_density


def gen_res_centric_env_for_list_of_pdbs_out_hd5f(pdb_fpath_list, output_dirpath=None, outfile_label='_', store_hdf_in_chunk=20):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    env_density_list = []
    res_label_list = []
    meta_data_dict = collections.OrderedDict()

    if output_dirpath == None:
        output_dirpath = os.getcwd()

    hdf_fpath = output_dirpath + '/test_hdf_file.h5'


    st_time_0 = time.perf_counter()

    counter = 0
    for index, pdb_fpath in enumerate(pdb_fpath_list):
        st_time_1 = time.perf_counter()
        print(index + 1, ' of ', len(pdb_fpath_list), ' PDBs ...')
        pdb_dir, pdb_name = os.path.split(pdb_fpath)
        pdb_name = str(pdb_name)

        pdb_model = collect_PDB_info(pdb_fpath)

        visited_pos = set()

        res_center_label_list = []
        res_center_density_list = []
        meta_data_list = []

        for ind, pdb_atom in enumerate(pdb_model['PDB_entries']):



            if pdb_atom.res_name in res_name_label_dict.keys():
                res_name = pdb_atom.res_name
                atom_type = pdb_atom.atom_type
                if atom_type == "CA":
                    chain_ID = pdb_atom.chain_id
                    res_num = pdb_atom.res_num
                    ctr_coord = (pdb_atom.x_coord, pdb_atom.y_coord, pdb_atom.z_coord)
                    if ctr_coord not in visited_pos:

                        counter += 1

                        visited_pos.add(ctr_coord)
                        ctr_pos_list = [ctr_coord, chain_ID, res_name, res_num]
                        meta_data = np.array([pdb_name, *ctr_coord, chain_ID, res_name, res_num], dtype=str)
                        # meta_data_str = ','.join([str(x) for x in meta_data])
                        res_center_density = gen_res_center_density_in_box(pdb_model, ctr_pos_list)
                        res_center_label = torch.tensor(res_name_label_dict[res_name], dtype=torch.long, device=device)

                        res_label_shape = res_center_label.shape
                        res_center_density_shape = res_center_density.shape


                        if counter == 1:
                            with h5py.File(hdf_fpath, 'w') as hf:
                                res_label_dataset = hf.create_dataset('res_label', data=res_center_label)
                                res_density_dataset = hf.create_dataset('res_center_density', data=res_center_density)
                                meta_data_dataset = hf.create_dataset('meta_data', data=meta_data)

                        else:
                            with h5py.File(hdf_fpath, 'a') as hf:
                                hf['res_label'].resize((hf['res_label'].shape[0] + res_center_label.shape[0]), axis=0)
                                hf['res_label'][-res_center_label.shape[0]] = res_center_label

                                hf['res_center_density'].resize((hf['res_center_density'].shape[0] + res_center_density.shape[0]), axis=0)
                                hf['res_center_density'][-res_center_density.shape[0]] = res_center_density

                                hf['meta_data'].resize((hf['meta_data'].shape[0] + 1), axis=0)
                                hf['meta_data'][-1] = meta_data






        time_took = time.perf_counter() - st_time_1
        print('Generated density data for ', pdb_name, ' in ', time_took)

    print('Generated data for all pdb in ', time.perf_counter() - st_time_0)

    # if write_output:
    #     print('Saving output')
    #     st_time3 = time.perf_counter()
    #     save_cpickle(meta_data_dict, output_dirpath + '/res_centric_metadata.bz2')
    #     save_cpickle(env_density_list, output_dirpath + '/res_centric_density.bz2')
    #     save_cpickle(res_label_list, output_dirpath + '/res_centric_label.bz2')
    #     print('Saved ... took ', time.perf_counter() - st_time3)
    #
    # print('Total time took ', time.perf_counter() - st_time_0)

    return env_density_list, res_label_list, meta_data_dict


def gen_res_centric_env_for_list_of_pdbs_out_one(pdb_fpath_list, output_dirpath=None, outfile_label='_'):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    env_density_list = []
    res_label_list = []
    meta_data_dict = collections.OrderedDict()

    if output_dirpath == None:
        output_dirpath = os.getcwd()

    env_density_fpath = output_dirpath + '/res_centric_density' + outfile_label + '.pk'
    res_label_fpath = output_dirpath + '/res_centric_label' + outfile_label + '.pk'
    meta_data_fpath = output_dirpath + '/res_centric_metadata' + outfile_label + '.pk'

    # env_density_file = open(env_density_fpath, 'wb')
    # res_label_file = open(res_label_fpath, 'wb')
    # meta_data_file = open(meta_data_fpath, 'wb')


    st_time_0 = time.perf_counter()

    for index, pdb_fpath in enumerate(pdb_fpath_list):
        st_time_1 = time.perf_counter()
        print(index + 1, ' of ', len(pdb_fpath_list), ' PDBs ...')
        pdb_dir, pdb_name = os.path.split(pdb_fpath)
        pdb_name = str(pdb_name)

        pdb_model = collect_PDB_info(pdb_fpath)

        visited_pos = set()

        for ind, pdb_atom in enumerate(pdb_model['PDB_entries']):
            if pdb_atom.res_name in res_name_label_dict.keys():
                res_name = pdb_atom.res_name
                atom_type = pdb_atom.atom_type
                if atom_type == "CA":
                    chain_ID = pdb_atom.chain_id
                    res_num = pdb_atom.res_num
                    ctr_coord = (pdb_atom.x_coord, pdb_atom.y_coord, pdb_atom.z_coord)
                    if ctr_coord not in visited_pos:
                        visited_pos.add(ctr_coord)
                        ctr_pos_list = [ctr_coord, chain_ID, res_name, res_num]
                        meta_data = [pdb_name] + ctr_pos_list
                        res_center_density = gen_res_center_density_in_box(pdb_model, ctr_pos_list)
                        res_center_label = torch.tensor(res_name_label_dict[res_name], dtype=torch.long, device=device)

                        with open(env_density_fpath, 'ab') as ed_file:
                            pickle.dump(res_center_density, ed_file)

                        with open(res_label_fpath, 'ab') as rd_file:
                            pickle.dump(res_center_label, rd_file)

                        with open(meta_data_fpath, 'ab') as md_file:
                            pickle.dump(meta_data, md_file)

                        # env_density_list.append(res_center_density)
                        # res_label_list.append(res_center_label)
                        # meta_data_list.append(meta_data)

                        # if pdb_name not in meta_data_dict.keys():
                        #     meta_data_dict[pdb_name] = ctr_pos_list
                        # else:
                        #     meta_data_dict[pdb_name].append(ctr_pos_list)

        time_took = time.perf_counter() - st_time_1
        print('Generated density data for ', pdb_name, ' in ', time_took)

    print('Generated data for all pdb in ', time.perf_counter() - st_time_0)

    # if write_output:
    #     print('Saving output')
    #     st_time3 = time.perf_counter()
    #     save_cpickle(meta_data_dict, output_dirpath + '/res_centric_metadata.bz2')
    #     save_cpickle(env_density_list, output_dirpath + '/res_centric_density.bz2')
    #     save_cpickle(res_label_list, output_dirpath + '/res_centric_label.bz2')
    #     print('Saved ... took ', time.perf_counter() - st_time3)
    #
    # print('Total time took ', time.perf_counter() - st_time_0)

    return env_density_list, res_label_list, meta_data_dict

@profile
def gen_res_centric_env_for_list_of_pdbs_out_one_dpath(pdb_fpath_list, output_dirpath=None, outfile_label='_'):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    env_density_list = []
    res_label_list = []
    meta_data_dict = collections.OrderedDict()

    if output_dirpath == None:
        output_dirpath = os.getcwd()

    data_dirpath = os.path.join(output_dirpath, 'res_center_data')
    if not os.path.isdir(data_dirpath):
        os.makedirs(data_dirpath)

    data_pickle_fname = 'res_center_data_label_meta_'


    st_time_0 = time.perf_counter()

    counter = 0

    for index, pdb_fpath in enumerate(pdb_fpath_list):
        st_time_1 = time.perf_counter()
        print(index + 1, ' of ', len(pdb_fpath_list), ' PDBs ...')
        pdb_dir, pdb_name = os.path.split(pdb_fpath)
        pdb_name = str(pdb_name)

        pdb_model = collect_PDB_info(pdb_fpath)

        visited_pos = set()

        for ind, pdb_atom in enumerate(pdb_model['PDB_entries']):
            if pdb_atom.res_name in res_name_label_dict.keys():
                res_name = pdb_atom.res_name
                atom_type = pdb_atom.atom_type
                if atom_type == "CA":
                    chain_ID = pdb_atom.chain_id
                    res_num = pdb_atom.res_num
                    ctr_coord = (pdb_atom.x_coord, pdb_atom.y_coord, pdb_atom.z_coord)
                    if ctr_coord not in visited_pos:
                        visited_pos.add(ctr_coord)
                        ctr_pos_list = [ctr_coord, chain_ID, res_name, res_num]
                        meta_data = [pdb_name] + ctr_pos_list
                        res_center_density = gen_res_center_density_in_box(pdb_model, ctr_pos_list)
                        res_center_label = torch.tensor(res_name_label_dict[res_name], dtype=torch.long, device=device)

                        data_dict = dict()
                        data_dict['res_center_density'] = res_center_density
                        data_dict['res_center_label'] = res_center_label
                        data_dict['meta_data'] = meta_data

                        pickle_save_name = data_pickle_fname + str(counter) + '.pk'
                        with open(os.path.join(data_dirpath, pickle_save_name), 'wb') as pk_file:
                            pickle.dump(data_dict, pk_file)

                        counter += 1


        time_took = time.perf_counter() - st_time_1
        print('Generated density data for ', pdb_name, ' in ', time_took)

    print('Generated data for all pdb in ', time.perf_counter() - st_time_0)

    # if write_output:
    #     print('Saving output')
    #     st_time3 = time.perf_counter()
    #     save_cpickle(meta_data_dict, output_dirpath + '/res_centric_metadata.bz2')
    #     save_cpickle(env_density_list, output_dirpath + '/res_centric_density.bz2')
    #     save_cpickle(res_label_list, output_dirpath + '/res_centric_label.bz2')
    #     print('Saved ... took ', time.perf_counter() - st_time3)
    #
    # print('Total time took ', time.perf_counter() - st_time_0)

    return env_density_list, res_label_list, meta_data_dict


def gen_res_centric_env_for_pdb(pdb_fpath):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    pdb_model = collect_PDB_info(pdb_fpath)

    center_pos_list = []
    center_res_env_density = []
    center_res_label = []
    visited_pos = set()

    for ind, pdb_atom in enumerate(pdb_model['PDB_entries']):
        if pdb_atom.res_name in res_name_label_dict.keys():
            res_name = pdb_atom.res_name
            atom_type = pdb_atom.atom_type
            if atom_type == "CA":
                chain_ID = pdb_atom.chain_id
                res_num = pdb_atom.res_num
                ctr_coord = (pdb_atom.x_coord, pdb_atom.y_coord, pdb_atom.z_coord)
                if ctr_coord not in visited_pos:
                    visited_pos.add(ctr_coord)
                    ctr_pos_list = [ctr_coord, chain_ID, res_name, res_num]
                    res_center_density = gen_res_center_density_in_box(pdb_model, ctr_pos_list)
                    center_res_env_density.append(res_center_density)
                    center_res_label.append(res_name_label_dict[res_name])
                    center_pos_list.append(ctr_pos_list)

    center_res_env_density = torch.stack(center_res_env_density)
    center_res_label = torch.tensor(center_res_label, device=device, dtype=torch.long)

    return center_res_env_density, center_res_label, center_pos_list


def gen_res_centric_env_for_list_pdbs(pdb_fpath_list, write_output=False, output_dirpath=None):

    env_density_data = []
    res_label_data = []

    meta_data_dict = collections.OrderedDict()

    st_time1 = time.perf_counter()
    for index, pdb_fpath in enumerate(pdb_fpath_list):
        st_time2 = time.perf_counter()
        print(index + 1, ' of ', len(pdb_fpath_list), ' PDBs ...')
        pdb_dir, pdb_name = os.path.split(pdb_fpath)
        pdb_name = str(pdb_name)
        res_env_density, res_label, center_res_pos_list = gen_res_centric_env_for_pdb(pdb_fpath)
        env_density_data.append(res_env_density)
        res_label_data.append(res_label)
        if pdb_name not in meta_data_dict.keys():
            meta_data_dict[pdb_name] = center_res_pos_list
        else:
            meta_data_dict[pdb_name].append(center_res_pos_list)

        print(pdb_name + ' done in ', time.perf_counter() - st_time2)

    env_density_data = torch.cat(env_density_data, dim=0)
    res_label_data = torch.cat(res_label_data, dim=0)

    print('Generated data in ', time.perf_counter() - st_time2)

    if write_output:
        print('Saving output')
        st_time3 = time.perf_counter()
        save_cpickle(meta_data_dict, output_dirpath + '/res_centric_metadata.bz2')
        save_cpickle(env_density_data, output_dirpath + '/res_centric_density.bz2')
        save_cpickle(res_label_data, output_dirpath + '/res_centric_label.bz2')
        print('Saved ... took ', time.perf_counter() - st_time3)
        # torch.save(env_density_data, output_dirpath + '/res_centric_density.pt')
        # torch.save(res_label_data, output_dirpath + '/res_centric_label.pt')

    return env_density_data, res_label_data, meta_data_dict


def gen_res_centric_env_for_list_pdbs_in_batch(pdb_fpath_list, batch_size=50, write_output=True, output_dirpath=None, outfile_label=''):

    num_pdb = len(pdb_fpath_list)
    batch_index = np.arange(start=0, stop=num_pdb, step=batch_size)

    for num in range(len(batch_index)):
        if num < len(batch_index)-1:
            print(batch_index[num], batch_index[num+1])
            pdb_list = pdb_fpath_list[batch_index[num]:batch_index[num+1]]
        else:
            print(batch_index[num])
            pdb_list = pdb_fpath_list[batch_index[num]:]

        env_density_data = []
        res_label_data = []
        meta_data_dict = collections.OrderedDict()

        for ind, pdb_fpath in enumerate(pdb_list):
            pdb_dir, pdb_name = os.path.split(pdb_fpath)
            pdb_name = str(pdb_name)

            print(ind+1, ' of ', len(pdb_list))
            print(pdb_name)

            res_env_density, res_label, center_res_pos_list = gen_res_centric_env_for_pdb(pdb_fpath)
            env_density_data.append(res_env_density)
            res_label_data.append(res_label)

            if pdb_name not in meta_data_dict.keys():
                meta_data_dict[pdb_name] = center_res_pos_list
            else:
                meta_data_dict[pdb_name].append(center_res_pos_list)

        env_density_data = torch.cat(env_density_data, dim=0)
        res_label_data = torch.cat(res_label_data, dim=0)

        if write_output:
            res_centric_denstiy_file = 'res_centric_density_' + outfile_label + '_' + str(num) + '.bz2'
            res_centric_label_file = 'res_centric_label_' + outfile_label + '_' + str(num) + '.bz2'
            res_centric_metadata_file = 'rec_centric_metadata_'+ outfile_label + '_' + str(num) + '.bz2'
            save_cpickle(env_density_data, output_dirpath + '/'+ res_centric_denstiy_file)
            save_cpickle(res_label_data, output_dirpath + '/' + res_centric_label_file)
            save_cpickle(meta_data_dict, output_dirpath + '/' + res_centric_metadata_file)

    return env_density_data, res_label_data, meta_data_dict


if __name__ == '__main__':
    pdb_dir = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure/trainset'
    output_dir = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure'
    pdb_fpath_list = glob.glob((pdb_dir + '/*.pdb'))
    # env_density_list, res_label_list, meta_data_dict = gen_res_centric_env_for_list_of_pdbs_out_one(pdb_fpath_list[:2], output_dirpath=output_dir)
    env_density_list, res_label_list, meta_data_dict = gen_res_centric_env_for_list_of_pdbs_out_one_dpath(pdb_fpath_list[:20],
                                                                                                    output_dirpath=output_dir)
    # env_density_list, res_label_list, meta_data_dict = gen_res_centric_env_for_list_of_pdbs_out_hd5f(pdb_fpath_list[:3],
    #                                                                                                 output_dirpath=output_dir)

    # env_density_tensor, res_label_tensor, metadata_dict = gen_res_centric_env_for_list_pdbs_in_batch(pdb_fpath_list[:3],
    #                                                                                                  batch_size=2,
    #                                                                                         write_output=True,
    #                                                                                         output_dirpath=output_dir,
    #                                                                                                  outfile_label='test')
