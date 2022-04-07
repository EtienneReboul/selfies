# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 16:13:51 2020

@author: jacqu

Download moses datasets and saves them 
"""

import moses
import pandas as pd
import os


def download_moses():
    script_dir = os.path.dirname(os.path.realpath(__file__))

    print('>>> Loading data from moses')
    train = moses.get_dataset('train')
    test = moses.get_dataset('test')
    test_scaffolds = moses.get_dataset('test_scaffolds')

    train = pd.DataFrame(train).rename(columns={0: 'smiles'})
    train.index.name='idx'
    test = pd.DataFrame(test).rename(columns={0: 'smiles'})
    test.index.name='idx'
    scaf = pd.DataFrame(test_scaffolds).rename(columns={0: 'smiles'})
    scaf.index.name='idx'

    print(scaf.head())

    print('>>> Saving data to csv files in ./data')
    train.to_csv(os.path.join(script_dir, 'datasets/moses_train.csv'))
    test.to_csv(os.path.join(script_dir, 'datasets/moses_test.csv'))
    scaf.to_csv(os.path.join(script_dir, 'datasets/moses_test_scaffolds.csv'))


download_moses()
