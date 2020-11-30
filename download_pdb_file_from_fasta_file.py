import os
import time
import json
import Bio
from Bio.PDB import PDBList


def make_directory(dirpath):

    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

    return dirpath


def save_dict_to_json(data_dict, fname, output_dirpath):

    outdirpath = make_directory(output_dirpath)
    json_fpath = os.path.join(outdirpath, fname+'.json')
    with open(json_fpath, 'w') as outfile:
        json.dump(data_dict, outfile)



def download_pdbs(fasta_fpath, pdb_outdirpath, json_dirpath):

    outdirpath = make_directory(pdb_outdirpath)

    pdb_id_dict = dict()
    pdb_id_dict['ids'] = []

    with open(fasta_fpath, 'r') as fastaf:
        num = 0
        fasta = fastaf.read().splitlines()
        for line in fasta:
            if line.startswith('>'):
                pdb_id = line.split(':')[0][1:]
                if pdb_id not in pdb_id_dict['ids']:
                    num += 1
                    pdb_id_dict['ids'].append(pdb_id)
                    download_pdbs_from_pdb_id(pdb_id, outdirpath)

    pdb_id_dict['total'] = num

    save_dict_to_json(pdb_id_dict, 'fasta_monomer_prots_pdbs', json_dirpath)

    return pdb_id_dict


def download_pdbs_from_pdb_id(pdb_id, save_dirpath):
    start_time = time.time()
    pdb_list = PDBList()
    pdb_list.retrieve_pdb_file(pdb_id, pdir=save_dirpath)
    print('--- %s seconds ----' % (time.time() - start_time))
    return pdb_list



if __name__=='__main__':

    fasta_fpath = '/Users/smd4193/Documents/deep_learning_prot_struct/fasta_monomer_prots.txt'

    outdirpath = os.path.join(os.path.split(fasta_fpath)[0], 'fasta_monomer_prots_pdbs')

    pdb_id_dict = download_pdbs(fasta_fpath, outdirpath, os.path.split(fasta_fpath)[0])