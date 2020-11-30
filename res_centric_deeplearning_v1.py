import torch
import plotly.express as px
import os
import glob
import numpy as np
import pickle as pk
# import _pickle as cpickle
# import bz2
import time
import plotly.graph_objects as go
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


# def save_cpickle(obj, fpath):
#     with bz2.BZ2File(fpath, 'w') as outfile:
#         cpickle.dump(obj, outfile)

# class ResCentricDataset(torch.utils.data.Dataset):
#
#     def __init__(self, data_tensor_fpath, label_tensor_fpath, transform=None):
#         """
#         initialize the dataset by giving the data and label tensor filepaths
#         :param data_tensor_fpath: data tensor filepath
#         :param label_tensor_fpath: label tensor filepath
#         :param transform: apply transform. If the data sizes are different sizes, transformation can be applied
#         to make the tensor of the same size
#         """
#         self.transform = transform
#
#         if torch.cuda.is_available():
#             self.device = torch.device("cuda")
#         else:
#             self.device = torch.device("cpu")
#
#         # the below if for previous data saving output
#         # self.data = self.load_bz_file(data_tensor_fpath).to(device=self.device)
#         # self.label = self.load_bz_file(label_tensor_fpath).to(device=self.device)
#
#         self.data = self.load_bz_file(data_tensor_fpath)
#         self.label = self.load_bz_file(label_tensor_fpath)
#
#         self.sample_num = len(self.label)
#
#     def __len__(self):
#         return self.sample_num
#
#     def __getitem__(self, index):
#         out_data = self.data[index].to(device=self.device)
#         out_label = self.label[index].to(device=self.device)
#         return out_data, out_label
#
#     def load_bz_file(self, fpath):
#         data = bz2.BZ2File(fpath, 'rb')
#         data = cpickle.load(data)
#         return data


class ResCentricDatasetFromDirPath(torch.utils.data.Dataset):

    def __init__(self, data_directory_path, endid='.pk', transform=None):
        """
        initialize the dataset by giving the data and label tensor filepaths
        :param data_tensor_fpath: data tensor filepath
        :param label_tensor_fpath: label tensor filepath
        :param transform: apply transform. If the data sizes are different sizes, transformation can be applied
        to make the tensor of the same size
        """
        self.transform = transform

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")

        self.data_fpaths = glob.glob(data_directory_path + '/*'+endid)[:200]

        self.sample_num = len(self.data_fpaths)

    def __len__(self):
        return self.sample_num

    def __getitem__(self, index):
        fpath = self.data_fpaths[index]
        with open(fpath, 'rb') as pk_file:
            data_dict = pk.load(pk_file)
        out_data = data_dict['res_center_density'].to(device=self.device)
        out_label = data_dict['res_center_label'].to(device=self.device)
        return out_data, out_label


class NNet(torch.nn.Module):
    """
    Neural net module
    """
    def __init__(self):
        super(NNet, self).__init__()
        self.relu = torch.nn.ReLU()
        self.pool = torch.nn.MaxPool3d(2)

        self.conv3_1 = torch.nn.Conv3d(4, 100, 3)
        self.conv3_2 = torch.nn.Conv3d(100, 200, 3)
        self.conv3_3 = torch.nn.Conv3d(200, 400, 3)

        self.fc1 = torch.nn.Linear(400 * 3 * 3 * 3, 1000)
        self.fc2 = torch.nn.Linear(1000, 100)
        self.fc3 = torch.nn.Linear(100, 20)

        self.dropout = torch.nn.Dropout(0.3)

    def forward(self, x):
        print('Input data size: ', x.shape)
        x = self.conv3_1(x)
        print('Conv1 outsize: ', x.shape)
        x = self.conv3_2(x)
        print('Conv2 outsize: ', x.shape)

        x = self.pool(x)
        print('MaxPool after conv2 outsize: ', x.shape)

        x = self.conv3_3(x)
        print('Conv3 outsize: ', x.shape)
        x = self.relu(x)
        print('Outsize after RELu: ', x.shape)
        x = self.dropout(x)
        print('Outsize after dropout: ', x.shape)
        x = self.pool(x)
        print('Maxpool after conv3 outsize: ', x.shape)

        x = x.view(x.size()[0], -1)

        x = self.fc1(x)
        print('Outsize FC1: ', x.shape)
        x = self.relu(x)
        x = self.dropout(x)

        x = self.fc2(x)
        print('OUtsize FC2: ', x.shape)
        x = self.relu(x)
        x = self.dropout(x)

        x = self.fc3(x)
        print('Outsize FC3: ', x.shape)
        print('heho')

        return x


# def gen_train_test_dataset(data_fpath, label_fpath, train_size_fraction=0.8):
#     """
#     generate training and test dataset
#     :param data_fpath: x data fpath
#     :param label_fpath: y data fpath
#     :param train_size_fraction: fraction of dataset to be used for trainsize
#     :return: train_dataset, test_dataset
#     """
#     dataset = ResCentricDataset(data_fpath, label_fpath)
#     train_size = int(dataset.sample_num * train_size_fraction)
#     test_size = dataset.sample_num - train_size
#
#     train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])
#
#     return train_dataset, test_dataset


def gen_train_test_dataset_v2(data_directory_path, train_size_fraction=0.8):
    """
    generate training and test dataset
    :param data_directory_path: directory path
    :param train_size_fraction: fraction of dataset to be used for trainsize
    :return: train_dataset, test_dataset
    """
    dataset = ResCentricDatasetFromDirPath(data_directory_path)
    train_size = int(dataset.sample_num * train_size_fraction)
    test_size = dataset.sample_num - train_size

    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])

    return train_dataset, test_dataset


def train_data(nnet, train_dataset_loader, optimizer, criterion):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    for (inputs, labels) in train_dataset_loader:
        inputs, labels = inputs.to(device, dtype=torch.float), labels.to(device, dtype=torch.long)
        optimizer.zero_grad()
        outputs = nnet(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

    return nnet


def test_nnet_with_data(nnet, dataset_loader, batch_size, optimizer, criterion, sum_loss=0, sum_correct=0, sum_total=0):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    label_list = []
    prediction_list = []

    for (inputs, labels) in dataset_loader:
        inputs, labels = inputs.to(device, dtype=torch.float), labels.to(device, dtype=torch.long)
        optimizer.zero_grad()
        outputs = nnet(inputs)
        loss = criterion(outputs, labels)
        sum_loss += loss.item()
        _, predicted = outputs.max(1)
        sum_total += labels.size(0)
        sum_correct += (predicted == labels).sum().item()

        label_list.append(labels)
        prediction_list.append(predicted)

    train_loss_val = sum_loss * batch_size/len(dataset_loader.dataset)
    train_acc_val = float(sum_correct/sum_total)

    return train_loss_val, train_acc_val, label_list, prediction_list


def flatten_label_list(label_list):
    label_list_ = [x.tolist() for x in label_list]
    flat_list = [x for sub_list in label_list_ for x in sub_list]
    return flat_list


def gen_pred_label_hits_sq_matrix(label_list, prediction_list):

    cnt_sq_mat = torch.zeros(20, 20, dtype=torch.int)

    for num in range(len(prediction_list)):

        label_, predict_ = label_list[num], prediction_list[num]
        cnt_sq_mat[label_, predict_] = cnt_sq_mat[label_, predict_] + 1

    return cnt_sq_mat


def plot_label_predict_matrix(label_list, prediction_list):

    flat_label_list = flatten_label_list(label_list)
    flat_prediction_list = flatten_label_list(prediction_list)

    sq_matrix = gen_pred_label_hits_sq_matrix(flat_label_list, flat_prediction_list)

    fig = px.imshow(sq_matrix)
    fig.show()


def plot_train_test_metric_over_epoch(nn_save_dict, output_dir, output_fname):
    x_ = np.arange(nn_save_dict['num_epochs'])

    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=x_, y=nn_save_dict['train_acc_value'], mode='lines+markers', name='train_accuracy'))
    fig1.add_trace(go.Scatter(x=x_, y=nn_save_dict['test_acc_value'], mode='lines+markers', name='test_accuracy'))
    fig1.write_html(output_dir + '/train_test_accuracy_' + output_fname + '.html')

    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=x_, y=nn_save_dict['train_loss_value'], mode='lines+markers', name='train_loss'))
    fig2.add_trace(go.Scatter(x=x_, y=nn_save_dict['test_loss_value'], mode='lines+markers', name='test_loss'))
    fig2.write_html(output_dir + '/train_test_loss_' + output_fname + '.html')


def write_nn_output(output_dir, epoch_num=0, train_acc=0.0, test_acc=0.0, train_loss=0.0, test_loss=0.0,
                    init_file=True):
    """

    :param output_dir:
    :param epoch_num:
    :param train_acc:
    :param test_acc:
    :param train_loss:
    :param test_loss:
    :param init_file:
    :return:
    """

    fpath = os.path.join(output_dir, 'nn_output.csv')
    if init_file:

        header = 'epoch,train_acc,test_acc,train_loss,test_loss\n'

        with open(fpath, 'w') as outfile:
            outfile.write(header)
            outfile.close()
    else:
        data_string = '{},{},{},{},{}\n'.format(epoch_num, train_acc, test_acc, train_loss, test_loss)
        with open(fpath, 'a') as outfile:
            outfile.write(data_string)
            outfile.close()


def make_new_dirpath(dirpath):

    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    return dirpath

@profile
def NN_module(train_dataset, test_dataset, batch_size, loader_num_worker, epoch_size, output_dir):

    # initialize data loader for train and test dataset

    train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True,
                                                    num_workers=loader_num_worker)
    test_data_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False,
                                                   num_workers=loader_num_worker)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    nnet = NNet()
    nnet = nnet.to(device)
    criterion = torch.nn.CrossEntropyLoss()
    optimizer = torch.optim.SGD(nnet.parameters(), lr=0.001, momentum=0.9, weight_decay=0.005)

    # make output directory path
    output_dir = make_new_dirpath(output_dir)

    # initialize the output file
    write_nn_output(output_dir=output_dir, init_file=True)


    ref_time1 = time.perf_counter()

    for epoch in range(epoch_size):

        start_time = time.perf_counter()
        print('\nEpoch ', epoch + 1)

        # dict for saving nn module
        nn_module_dict = dict()

        nnet = train_data(nnet, train_data_loader, optimizer, criterion)

        train_loss_val, train_acc_val, train_label, train_pred = test_nnet_with_data(nnet, train_data_loader,
                                                                                     batch_size, optimizer, criterion)

        print('Train accuracy ', train_acc_val)

        test_loss_val, test_acc_val, test_label, test_pred = test_nnet_with_data(nnet, test_data_loader, batch_size,
                                                                                 optimizer, criterion)

        print('Test accuracy ', test_acc_val)

        elapsed = time.perf_counter() - start_time
        print('Time took: ', elapsed)

        # write the output file

        write_nn_output(output_dir=output_dir,
                        epoch_num=epoch+1,
                        train_acc=train_acc_val,
                        test_acc=test_acc_val,
                        train_loss=train_loss_val,
                        test_loss=test_loss_val,
                        init_file=False)

        nn_module_dict['nnet'] = nnet
        nn_module_dict['epoch_num'] = epoch + 1
        nn_module_dict['train_acc_value'] = train_acc_val
        nn_module_dict['test_acc_value'] = test_acc_val
        nn_module_dict['train_loss_value'] = train_loss_val
        nn_module_dict['test_loss_value'] = test_loss_val
        nn_module_dict['train_label_list'] = train_label
        nn_module_dict['train_pred_list'] = train_pred
        nn_module_dict['test_label_list'] = test_label
        nn_module_dict['test_pred_list'] = test_pred

        nn_module_dict_fname = 'nn_module_epoch_' + str(epoch+1) + '.pk'
        with open(os.path.join(output_dir, nn_module_dict_fname), 'wb') as pk_file:
            pk.dump(nn_module_dict, pk_file)

    print('Total time: ', time.perf_counter() - ref_time1)



if __name__ == '__main__':

    data_directory_path = '/Users/smd4193/Documents/deep_learning_prot_struct/PDBstructure/res_center_data'
    output_dirpath = data_directory_path + '/nn_output'

    # train_dataset, test_dataset = gen_train_test_dataset(input_data_fpath, input_label_fpath)
    train_dataset, test_dataset = gen_train_test_dataset_v2(data_directory_path)

    print('Train data size ', len(train_dataset))
    print('Test data size ', len(test_dataset))

    NN_module(train_dataset=train_dataset,
                            test_dataset=test_dataset,
                            batch_size=2,
                            loader_num_worker=4,
                            epoch_size=2,
                            output_dir=output_dirpath)