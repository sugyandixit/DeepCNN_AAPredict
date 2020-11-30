def find_pos_in_grid(kd_tree, pos, PDB_entries, num_neighbors=1):
    dist, loc = kd_tree.query(pos, k=num_neighbors)
    pdb_entries = PDB_entries[loc]
    return pdb_entries


def get_min_and_max(data):
    min_data = min(data)
    max_data = max(data)
    return min_data, max_data



def raster_axes(axis_min, num_pts, size):
    grid = []

    k = 0
    pos = axis_min
    while (k < num_pts):
        grid.append(pos)
        pos = pos + size
        k = k + 1

    return grid



def get_box_min_max(box_size):
    box_min = - box_size / 2
    box_max = + box_size / 2
    return box_min, box_max



def gen_grid_points_with_x_y_z_samples(all_x, all_y, all_z, grid_size=10):
    x_min, x_max = get_min_and_max(all_x)
    y_min, y_max = get_min_and_max(all_y)
    z_min, z_max = get_min_and_max(all_z)

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    num_grid_x = x_range / grid_size
    num_grid_y = y_range / grid_size
    num_grid_z = z_range / grid_size

    # get grids

    x_grids = raster_axes(x_min, num_grid_x, grid_size)
    y_grids = raster_axes(y_min, num_grid_y, grid_size)
    z_grids = raster_axes(z_min, num_grid_z, grid_size)

    grid_pos = []

    for i in range(0, len(x_grids)):
        for j in range(0, len(y_grids)):
            for k in range(0, len(z_grids)):
                grid_pos.append([x_grids[i], y_grids[j], z_grids[k]])

    return grid_pos




def points_to_Xsmooth_with_translate(protein_model_dict, center_pos_list, channels_list=['N', 'C', 'O', 'S'], atom_density=0,
                      pixel_size=1, box_size=20):
    bias = [-0.558, -0.73, 1.226]  ## ???? why is there a bias?

    num_3d_pixels = int(box_size / pixel_size)
    center_pos, chain_ID, center_res_name, center_res_num, = center_pos_list
    all_atom_type = protein_model_dict['all_atom_type']
    backbone = protein_model_dict['ID_dict'][chain_ID][0:4]
    deleted_res = protein_model_dict['ID_dict'][chain_ID][4:]
    deleted_res_index = [atom.atom_index for atom in deleted_res]
    backbone_res_index = [atom.atom_index for atom in backbone]

    box = []
    box_ori = []
    reference = []
    new_pos_in_box = []
    atom_count = 0
    valid_box = False
    atom_loc_box = []
    atom_loc_box_translate = []

    box_x_min, box_x_max = get_box_min_max(box_size)
    box_y_min, box_y_max = get_box_min_max(box_size)
    box_z_min, box_z_max = get_box_min_max(box_size)

    center_pos_dict = get_position_dict(protein_model_dict['ID_dict'][chain_ID], center_res_num)

    ref_pos, transform_mat = center_and_transform(center_res_name, center_pos_dict)

    all_pos = np.array(protein_model_dict['all_coords'])
    transformed_pos = ((all_pos - ref_pos).dot(transform_mat))  # - bias ##todo: why do we need a bias???
    # plot_3d_coords_comp_transform(all_pos, transformed_pos) # to investigate if coords are being transformed

    x_index = np.intersect1d(np.where(transformed_pos[:, 0] > box_x_min), np.where(transformed_pos[:, 0] < box_x_max))
    y_index = np.intersect1d(np.where(transformed_pos[:, 1] > box_x_min), np.where(transformed_pos[:, 1] < box_x_max))
    z_index = np.intersect1d(np.where(transformed_pos[:, 2] > box_x_min), np.where(transformed_pos[:, 2] < box_x_max))

    final_index = np.intersect1d(x_index, y_index)
    final_index = np.intersect1d(final_index, z_index)
    final_index = [x for x in final_index]
    # final_index = [ind for ind in final_index if ind not in deleted_res_index]
    # final_index = [ind for ind in final_index if (all_atom_type[ind] == 'N' or all_atom_type[ind] == 'CA' or all_atom_type == 'C' or all_atom_type == 'O')]

    box_ori = [protein_model_dict['PDB_entries'][i] for i in final_index]
    new_pos_in_box = transformed_pos[final_index]
    atom_count = len(box_ori)
    threshold = (box_size ** 3) * atom_density

    # translate grid points
    translate_vector = [(0, 0, 0), (1, 1, 1), (-1, -1, -1),(1, 2, 1),  (-1, -2, -1), (2, 2, 2), (-2, -2, -2),
                        (2, 3, 2), (-2, -3, -2), (3, 3, 3), (-3, -3, -3), (3, 2, 1), (-3, -2, -1)]

    translate_samples = [[] for _ in range(len(translate_vector))]
    translate_xsmooth = [[] for _ in range(len(translate_vector))]
    translate_loc_box = [[] for _ in range(len(translate_vector))]
    translate_atoms_loc_box = [[] for _ in range(len(translate_vector))]

    for tv_num, tranlsate_v in enumerate(translate_vector):

        if atom_count > threshold:
            valid_box = True

            sample = np.zeros((len(channels_list), num_3d_pixels, num_3d_pixels, num_3d_pixels))

            for index in range(0, len(box_ori)):
                box_atoms = box_ori[index]

                x_pos = new_pos_in_box[index][0]
                y_pos = new_pos_in_box[index][1]
                z_pos = new_pos_in_box[index][2]

                x_new = x_pos - box_x_min
                y_new = y_pos - box_y_min
                z_new = z_pos - box_z_min

                bin_x = int(np.floor(x_new / pixel_size))
                bin_y = int(np.floor(y_new / pixel_size))
                bin_z = int(np.floor(z_new / pixel_size))

                if bin_x <= 15 and bin_x >= 5:
                    if bin_y <= 15 and bin_y >= 5:
                        if bin_z <= 15 and bin_z >= 5:
                            bin_x = bin_x + tranlsate_v[0]
                            bin_y = bin_y + tranlsate_v[0]
                            bin_z = bin_z + tranlsate_v[0]
                            translate_atoms_loc_box[tv_num].append([(bin_x, bin_y, bin_z), box_atoms])

                translate_loc_box[tv_num].append([(bin_x, bin_y, bin_z), box_atoms])

                if bin_x == num_3d_pixels:
                    bin_x = num_3d_pixels - 1
                if bin_y == num_3d_pixels:
                    bin_y = num_3d_pixels - 1
                if bin_z == num_3d_pixels:
                    bin_z = num_3d_pixels - 1

                for ind, chan in enumerate(channels_list):
                    if box_atoms.atom_id == chan:
                        sample[ind, bin_x, bin_y, bin_z] = sample[ind, bin_x, bin_y, bin_z] + box_atoms.value
                        break

            X_smooth = np.zeros(sample.shape)

            for num in range(0, len(channels_list)):
                X_smooth[num, :, :, :] = gaussian_filter(sample[num, :, :, :], sigma=0.6, order=0, output=None,
                                                         mode='reflect', cval=0.0, truncate=4.0)
                X_smooth[num, :, :, :] *= 1000

            translate_samples[tv_num].append(sample)
            translate_xsmooth[tv_num].append(X_smooth)


    translate_density_diff_list = [[] for _ in range(len(translate_xsmooth))]
    translate_density_diff_percent_list = [[] for _ in range(len(translate_xsmooth))]

    ref_xsmooth = translate_xsmooth[0][0]

    sum_densities = [[] for _ in range(len(translate_xsmooth))]

    for xs_num in range(len(translate_xsmooth)):
        curr_xsmooth = translate_xsmooth[xs_num][0]

        diff_arr = np.zeros(curr_xsmooth.shape)
        diff_percent_arr = np.zeros(curr_xsmooth.shape)

        sum_density_channel = [[] for _ in range(len(channels_list))]



        for ch_num in range(0, len(channels_list)):

            sum_density_channel = [[] for _ in range(len(channels_list))]
            sum_density = 0

            box_l = len(curr_xsmooth[ch_num])
            for box_num in range(len(curr_xsmooth[ch_num])):



                ref_curr_arr = curr_xsmooth[ch_num][box_num]
                ref_arr = ref_xsmooth[ch_num][box_num]
                difference = ref_curr_arr - ref_arr
                diff_percent = (difference/ref_arr) * 100
                diff_arr[ch_num][box_num] = curr_xsmooth[ch_num][box_num] - ref_xsmooth[ch_num][box_num]
                diff_percent_arr[ch_num][box_num] = diff_percent

                sum_density += np.sum(ref_curr_arr)

            sum_density_channel[ch_num].append(sum_density)


        translate_density_diff_list[xs_num].append(diff_arr)
        translate_density_diff_percent_list[xs_num].append(diff_percent_arr)
        print('heho')


    #plot some figures

    fig1 = plot_Xsmooth(ref_xsmooth[1])
    fig1.update_layout(title_text='Ref Density')
    fig1.write_html('/Users/smd4193/Documents/deep_learning_prot_struct/translate_grid_figs/_ref_density.html')

    for num_k in range(len(translate_vector)):

        translate_xsmooth_diff_percent = translate_density_diff_list[num_k][0]
        translate_xsmooth_diff_percent_ = translate_density_diff_percent_list[num_k][0]

        translate_density = translate_xsmooth[num_k][0]

        fig2 = plot_Xsmooth(translate_density[1])
        fig3 = plot_Xsmooth(translate_xsmooth_diff_percent[1])
        fig4 = plot_Xsmooth(translate_xsmooth_diff_percent_[1])

        fig2.update_layout(title_text='Translate Density | Translate Vector (%d, %d, %d)' % (
        translate_vector[num_k][0], translate_vector[num_k][1], translate_vector[num_k][1]))
        fig3.update_layout(title_text='Difference Density | Translate Vector (%d, %d, %d)' % (
            translate_vector[num_k][0], translate_vector[num_k][1], translate_vector[num_k][1]))
        fig4.update_layout(title_text='Difference Density percent | Translate Vector (%d, %d, %d)' % (
            translate_vector[num_k][0], translate_vector[num_k][1], translate_vector[num_k][1]))

        fig2.write_html('/Users/smd4193/Documents/deep_learning_prot_struct/translate_grid_figs/' + str(num_k) + '_translate_density.html')
        fig3.write_html('/Users/smd4193/Documents/deep_learning_prot_struct/translate_grid_figs/' + str(
            num_k) + '_density_diff.html')
        fig4.write_html('/Users/smd4193/Documents/deep_learning_prot_struct/translate_grid_figs/' + str(
            num_k) + '_density_diff_percent.html')



    print('heho')



def generate_input_coords_channels_labels_e3nn(list_of_pdb_fpaths, broad_channel='atom_id',
                                               sub_channels=['C', 'N', 'O', 'S']):
    """
    create
    :param pdb_fpath:
    :param broad_channel:
    :param sub_channels:
    :return:
    """

    mol_coords = []
    mol_species = []
    mol_atom_type_set = set()
    mol_labels = []

    for pdb_fpath in list_of_pdb_fpaths[:100]:
        pdb_model = collect_PDB_info(pdb_fpath)

        mol_coords.append(torch.tensor(pdb_model['all_coords']))

        atom_species = list(map(str, pdb_model['all_atom_id']))
        mol_atom_type_set = mol_atom_type_set.union(set(atom_species))
        mol_species.append(atom_species)

        mol_labels.append(torch.tensor(pdb_model['all_res_label']))

    mol_atom_type_list = sorted(list(mol_atom_type_set))

    A = len(mol_atom_type_list)
    mol_Rs_in = [(A, 0, 1)]

    mol_features = []
    for coords, species in zip(mol_coords, mol_species):
        N, _ = coords.shape
        feature = torch.zeros(N, A, device='cuda:0')
        atom_type_indicies = [mol_atom_type_list.index(specie) for specie in species]
        feature[range(N), atom_type_indicies] = 1
        mol_features.append(feature)

    mol_label_randn = torch.randn(len(mol_features), 1)

    print('heho')
