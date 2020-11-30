import os
import pickle
import io
import tarfile


output_dir = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure'
pk_fname = 'res_centric_density_.pk'

targz_fpath = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure/res_center_data/resdata.tar.gz'

tar_fpaths = []

with tarfile.open(tar_fpaths, mode='r:gz') as outf:
    print('heho')


# data = []
#
# with open(os.path.join(output_dir, pk_fname), 'rb') as pk_file:
#     obj = pickle.load(pk_file)
#
#     print('heho')
#     # for _ in range(pickle.load(pk_file)):
#     #     data.append(pickle.load(pk_file))
