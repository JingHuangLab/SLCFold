#!/usr/bin/python

"""
Fitness score for new fold 
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np

from pandas import read_csv
from multiprocessing import Pool
from Bio.PDB import PDBParser, Selection
from Bio.PDB.DSSP import DSSP, ss_to_index

TMALIGN='/home/xiety/software/TMalign/TMalign'
TMSCORE='/home/xiety/software/TMscore/TMscore'
PDBDIR='../PDBs'
NUMPROCESS=4
def score(args):
    weight = 0.5
    if C_3d > 0.5:
        weight = S_3d / 0.5
    else:
        weight = 1

    score_ = 0.5 * weight * S_3d * ( 1 + (C_max_all - C_3d) * (1 - 0.5 * (C_2d_ave + C_1d_ave)))
    return score_

def TMalign(f1,f2, ifTMscore=False):
    TMscore = -1
    if ifTMscore:
        grep_string = '"TM-score    ="'
        result = os.popen(TMSCORE + " " + f1 + ' ' + f2 +
                      '|grep ' + grep_string + "|awk \'{ print $3 }'")
    else:
        grep_string = '"if normalized by length of Chain_2"'
        result = os.popen(TMALIGN + " -byresi 0 -a T  " + f2 + ' ' + f1 +
                          '|grep ' + grep_string + "|awk \'{ print $2 }'")
    TMscore = result.read()
    if len(TMscore) > 0:
        TMscore = float(TMscore)
    else:
        TMscore = -1
    return TMscore

def comp_prot2PDB(f):
    TMscores = []
    PDBs = os.listdir(PDBDIR)
    for PDB in PDBs:
        PDB_path = os.path.join(PDBDIR, PDB)
        TMscore = TMalign(f, PDB_path)
        print("Compare " + str(f) + " " + str(PDB_path) + ' ' + str(TMscore))
        TMscores.append(TMscore)

    max_TMscore = max(TMscores)
    index = TMscores.index(max_TMscore)
    most_similiar_PDB = PDBs[index]
    return max_TMscore, most_similiar_PDB
    
def get_TMscore_pre_exp(args, protein):
    prot_dir = os.path.join(args.input, protein) 
    pdb_paths = []
    for f in os.listdir(prot_dir):
        if f.endswith(".pdb"):
            f_path = os.path.join(prot_dir, f)
            pdb_paths.append(f_path)
            # compare the pdb file to all PDB structures

    #max_TMscore, most_similiar_PDB = comp_prot2PDB(f_path)
    with Pool(processes=NUMPROCESS) as pool:
        data = pool.map(comp_prot2PDB, pdb_paths)
    max_TMscores = [i[0] for i in data]
    max_TMscore = max(max_TMscores)
    index = max_TMscores.index(max_TMscore)
    most_similiar_PDB = data[index][1]
    return max_TMscore, most_similiar_PDB

def get_TMscore_pre_pre(args, protein):
    prot_dir = os.path.join(args.input, protein) 
    pdb_paths = []
    for f in os.listdir(prot_dir):
        if f.endswith(".pdb"):
            f_path = os.path.join(prot_dir, f)
            pdb_paths.append(f_path)

    max_TMscore = 0
    len_pdbs = len(pdb_paths)
    for i in range(len_pdbs-1):
        for j in range(i+1, len_pdbs):
            TMscore = TMalign(pdb_paths[i], pdb_paths[j], ifTMscore=True)
            print("Compare " + str(pdb_paths[i]) + " " + str(pdb_paths[j]) + ' ' + str(TMscore))
            TMscore = float(TMscore)
            if TMscore > max_TMscore:
                max_TMscore = TMscore

    return max_TMscore

def get_precision(contacts, topN, res, contact_gap=6, cutoff = 8):
    n_gap_contact= 0
    n_hit = 0

    for contact in contacts:
        index1 = int(contact[0] - 1)
        index2 = int(contact[1] - 1)

        if abs(index1 - index2) > contact_gap:
            n_gap_contact = n_gap_contact + 1
            dis = res[index1]['CA'] - res[index2]['CA']
            
            if dis < cutoff:
                n_hit = n_hit + 1

        if n_gap_contact > topN:
            break
    prec = n_hit * 1.0 / topN

    return prec

def get_residues(name, path):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(name, path)
    residues = Selection.unfold_entities(struct, "R") 
    return residues

def get_contact_consist(args, protein):
    contact_path = None
    prot_dir = os.path.join(args.input, protein) 
    pdb_paths = []
    for f in os.listdir(prot_dir):
        if f.endswith(".contact"):
            contact_path = os.path.join(prot_dir,f)
        if f.endswith(".pdb"):
            pdb_paths.append(os.path.join(prot_dir,f))
    if contact_path is None:
        return 0
    else:
        contacts = np.genfromtxt(contact_path, delimiter=' ')
        precisions = []
        for pdb_path in pdb_paths:
            res = get_residues(protein, pdb_path)
            topN = len(res)
            prec = get_precision(contacts, topN, res)
            print("Contact consistency for " + pdb_path + ": " + str(prec))
            precisions.append(prec)

        mean_prec = np.mean(precisions)
        return mean_prec
def read_ss(path):
    if os.path.exists(path):
        ss = np.genfromtxt(path, comments='#')[:, -3:]
        return ss
    else:
        return []

def  ss_predict_struct(protein, path):
        p = PDBParser(QUIET=True)
        structure = p.get_structure(protein, path)
        model = structure[0]
        dssp = DSSP(model, path, dssp='dssp')
        ss = []
        for i in dssp:
            ss_temp = i[2]
            if i[2] in ['G','H','I']:
                ss_temp = ss_to_index('H')
            elif i[2] in ['E','B']:
                ss_temp = ss_to_index('E')
            elif i[2] in ['S','T','-']:
                ss_temp = ss_to_index('C')
            else:
                print(ss_temp + ' is not belong to 8 classes in dssp')
                exit(1)
                
            ss.append(ss_temp)        
        return ss

def ss_consistency(ss_predict, ss_predict_struct_):
    """
    calculate consistency of ss from struct and predicted ss
    """
    n_seq = len(ss_predict)
    score = 0
    for i in range(n_seq):
        #score_temp = ss_predict[i, ss_predict_struct_[i]] 
        score_temp = abs(2 * ss_predict[i, ss_predict_struct_[i]] - np.sum(ss_predict[i]))
        score = score + score_temp
    return score / n_seq

def get_ss_consist(args, protein):
    ss3_path = None
    prot_dir = os.path.join(args.input, protein) 
    pdb_paths = []
    for f in os.listdir(prot_dir):
        if f.endswith(".ss3"):
            ss3_path = os.path.join(prot_dir,f)
        if f.endswith(".pdb"):
            pdb_paths.append(os.path.join(prot_dir,f))
    if ss3_path is None:
        return 0
    else:
        ss_predict = read_ss(ss3_path)

        scores = []
        for pdb_path in pdb_paths:
            ss_predict_struct_ = ss_predict_struct(protein, pdb_path)
            score = None
            if len(ss_predict) != 0 and len(ss_predict_struct_) != 0:
                score = ss_consistency(ss_predict, ss_predict_struct_)
            else:
                score = 0
            scores.append(score)
            print("Secondary structure consistency for " + pdb_path + ": " + str(score))

        ave_score = np.mean(scores)
        return ave_score
        # return scores[0] # in our protocol, only trRosetta predicted structures used (with suffix: partial.pdb) for secondary structure consistency, as RaptorX-3DModeling adopts the secondary structures generated by DeepCNF. 

def get_weight(S_3d, C_3d):
    if C_3d > 0.5:
        weight = S_3d / 0.5
        return weight
    else:
        return 1

def term_collection(protein, args):
    S_3d, most_similiar_PDB = get_TMscore_pre_exp(args,protein) # compare all predicted structures to all pdb structures
    C_3d = get_TMscore_pre_pre(args,protein) # compare pairwise predicted structures
    C_2d_ave = get_contact_consist(args,protein) # get the consistency of concat map and predicted 3d structures
    C_1d_ave = get_ss_consist(args,protein) # get the consistency of secondary structures and predicted 3d structures

    weight = get_weight(S_3d, C_3d)
    return S_3d, C_3d, C_2d_ave, C_1d_ave, weight

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input', help='directory of target proteins')
    parser.add_argument('--term', action='store_true', help='calculate each term of the fitness score')
    parser.add_argument('--score', action='store_true', help='calculate fitness score')
    parser.add_argument('--fixscore', action='store_true', help='calculate fitness score')

    args = parser.parse_args()

    if args.term:
        df = pd.DataFrame(columns=['protein', 'S_3d', 'C_3d', 'C_2d_ave', 'C_1d_ave', 'weight'])
        for dir_ in os.listdir(args.input):
            S_3d, C_3d, C_2d_ave, C_1d_ave, weight = term_collection(dir_, args)
            df = df.append({'protein': dir_,
                            'S_3d':S_3d,
                            'C_3d':C_3d,
                            'C_2d_ave':C_2d_ave,
                            'C_1d_ave':C_1d_ave,
                            'weight':weight}, 
                            ignore_index=True)
        df.to_csv("../output/score_items.csv")
    if args.score:
        data = pd.read_csv("../output/score_items.csv", index_col=0)
        C_max_all = data.max()['C_3d']
        # C_max_all = 0.9583 # the maximum TMscore of predicted structure for all SLCs
        print("C_max_all: " + str(C_max_all))
        data['score'] = 1 / 2 * data['weight'] * data['S_3d'] * ( 1 + ( C_max_all - data['C_3d']) * ( 1 - 1 / 2 * ( data['C_2d_ave'] + data['C_1d_ave'])))
        # data['score'] = 1 / 2 * data['weight'] * data['S_3d'] * ( 1 + ( C_max_all - data['C_3d']) * ( 1 - 2 / 3 * ( data['C_2d_ave'] + 1 / 2 * data['C_1d_ave']))) # the first version fitness score

        data.to_csv("../output/score.csv", index=False)
