import os
import bz2
import _pickle as cpickle
import glob
import torch


def load_bz_file(fpath):
    data = bz2.BZ2File(fpath, 'rb')
    data = cpickle.load(data)
    return data

def get_res_centric_data(dirpath):

    density_data_files = glob.glob(dirpath + '/*density*.bz2')
    label_data_files = glob.glob(dirpath + '/*label*.bz2')
    metadata_files = glob.glob(dirpath + '/*metadata*.bz2')

    return density_data_files, label_data_files, metadata_files


def save_cpickle(obj, fpath):

    with bz2.BZ2File(fpath, 'w') as outfile:
        cpickle.dump(obj, outfile)


def combine_torch_files(list_files, output_dir):

    tensor_data = []

    for file in list_files:

        data = load_bz_file(file)
        tensor_data.append(data)

    out_tensor = torch.cat(tensor_data, dim=0)

    save_cpickle(out_tensor, )