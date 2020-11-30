import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import _pickle as cpickle
import bz2


def load_bz_file(fpath):
        data = bz2.BZ2File(fpath, 'rb')
        data = cpickle.load(data)
        return data


def plot_bar_graph(ydata, xlabel, ylabel):

    x = np.arange(len(ydata))

    width = 0.3

    fig, ax = plt.subplots()

    bar = ax.bar(x, ydata, width)

    plt.xticks(x, xlabel)

    plt.ylabel(ylabel)

    plt.show()




def load_pickle_obj(fpath):

    with open(fpath, 'rb') as pk:
        obj = pickle.load(pk)

    return obj


output_dir = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure'
dict_fname = 'res_centric_metadata_50.bz2'

data_dict = load_bz_file(os.path.join(output_dir, dict_fname))

res_list = []

for key in data_dict.keys():

    for list in data_dict[key]:

        res_list.append(list[2])

count_res = Counter(res_list)
print(len(res_list))

res_label = []
res_count = []

for key in count_res.keys():
    res_label.append(key)
    res_count.append(count_res[key])

plot_bar_graph(res_count, res_label, 'count')

print('heho')