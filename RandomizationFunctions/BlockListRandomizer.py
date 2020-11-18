# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 13:38:03 2020

@author: Julius Kricheldorff MSc.

@affilation: University of Oldenburg

@contact: julius@kricheldorff.de

Pseudo-Randomization Algorithm for the adaptive control project. The pseudo-
randomization algorithm is a very rough approximation (i.e. my implementation)
based on the algorithm idea by  Maarten van Casteren and Mathew Davis. (See van 
Casteren, M., & Davis, M. H. (2006). Mix, a program for pseudorandomization. 
Behavior research methods, 38(4), 584-589.). The code for the backtracking 
was copied from Ray Toal (https://cs.lmu.edu/~ray/notes/backtracking/).
"""

import random
import pandas as pd
import numpy as np
from time import time
from math import floor

class BlockListRandomizer(object):
    def __init__(self, Trial_Nr = 268, Block_Nr = 4, Ind_prop = 0.8,
                 Inducer_Nr = 188, IWSPC_feature = 'distance', size_L = 300,
                 size_M = 260):
        self.trlnum = Trial_Nr
        self.Block_Nr = Block_Nr
        self.inducer_nr = Inducer_Nr
        self.diag_nr = self.trlnum - self.inducer_nr
        if IWSPC_feature == 'distance':
            # the manipulated features will be either distances one and two
            self.IWSPC = IWSPC_feature 
        elif IWSPC_feature == 'size':
            # the manipulated feature will be low and high numbers respectively
            self.IWSPC = IWSPC_feature 
        self.IWSPC = IWSPC_feature 
        self.Order = [['MI', 'MC'], ['IWPC_1', 'IWPC_2']]
        random.shuffle(self.Order)
        self.Correct = ['Left', 'Right']
        self.Congruence = ['congruent', 'incongruent']
        self.item_pool_dist_1 = [(1,2), (3,4), (6,7), (8,9)]
        self.item_pool_dist_2 = [(1,3), (2,4), (6,8), (7,9)]
        # self.item_pool_dist_1 = ['A1', 'B1', 'C1', 'D1']
        # self.item_pool_dist_2 = ['A2', 'B2', 'C2', 'D2']
        self.item_number =len(self.item_pool_dist_1)+len(self.item_pool_dist_2)
        self.inducer_rep = int(Inducer_Nr/(self.item_number/2))
        self.diagnostic_rep = int((Trial_Nr-Inducer_Nr)/(self.item_number/2))
       
    def create_LWPC_sets(self, blockType):
        try:
            if blockType is 'MI':
                main  = 'incongruent'
                minor = 'congruent'
            elif blockType is 'MC':
                main = 'congruent'
                minor = 'incongruent'
        except NameError:
            print('No blockType has been specified. Please use either',
                                 'the string MI for a mainly incongruent- or', 
                                 ' MC for mainly congruent block.')
        # function that creates the basic LWPC List layout
        # First we shuffle the two item lists so that items are randomly assigned
        # to inducer and diagnostic items respectively.
        D1 = list(self.item_pool_dist_1)
        D2 = list(self.item_pool_dist_2)
        random.shuffle(D1)
        random.shuffle(D2)
        # Assign the randomized items to diagnostic and inducer item lists
        LWPC_ind = (D1[:2] + D2[:2])*self.inducer_rep
        LWPC_diag = (D1[2:] + D2[2:])*self.diagnostic_rep
        Item_list =  LWPC_ind + LWPC_diag 
        # Item Type - Inducer Diagnostic Item
        Item_type = list(['Inducer'])*self.inducer_rep*int(self.item_number/2) + list(['Diagnostic'])*self.diagnostic_rep*int(self.item_number/2)
        # Split the sets in equally sized Blocks
        Block_Nr = [1]*int(self.inducer_nr/2) + [2]*int(self.inducer_nr/2) + [1]*int(self.diag_nr/2) + [2]*int(self.diag_nr/2)
        # Now throw together an empty data frame with columns correct response, 
        # congruence, ID and block
        LWPC_list = pd.DataFrame({'Item': Item_list, 'Trl_type': Item_type, 
                               'Correct': [None]*len(Item_list),
                                 'Congruency': [None]*len(Item_list),
                                 'Analysis_type': [None]*len(Item_list),
                                 'Block': Block_Nr,
                                 'ID': range(len(Item_list))})
        LWPC_list.Analysis_type = blockType
        # now assgin biased congruency conditions to the inducer items
        # First shuffle correct so left and right are not biased across blocks
        
        
        '''
        Note: this is probably the ugliest bit of coding I ever did. Sorry for 
        the eyesore! 
        '''
        correct_vals = self.Correct*2
        correct_vals2 = self.Correct
        # First we add Proportions of congruent and incongruent items for each of
        # the inducer items with congruent and incongruent choices being balanced
        # by presentation side (Correct)
        for block in set(LWPC_list.Block):
            count = 0
            random.shuffle(correct_vals2)
            for ind_item in set(LWPC_ind):
                # get a subset of the data for each item
                Subset = LWPC_list.index[(LWPC_list.Item == ind_item) & (LWPC_list.Block == block)]
                # then add correct and congruency values to the dataframe
                for num, rep in enumerate(Subset):
                    if num%2 == 0:
                        LWPC_list.loc[rep, 'Correct'] = self.Correct[0]
                    if num%2 != 0:
                        LWPC_list.loc[rep, 'Correct'] = self.Correct[1]
                    if num%5 == 0:
                        LWPC_list.loc[rep, 'Congruency'] = minor
                    if num%5 != 0:
                        LWPC_list.loc[rep, 'Congruency'] = main
                        # if we cannot assign correct and incorrect equally because
                        # an item is not presented an even number of times (must 
                        # happen for two items our chosen trial numbers). The last item is
                        # chosen randomly between correct and incorrect to make 
                        # sure the block is overall still balanced. Hard coded and
                        # ugly I know.
                    if num == len(Subset)-1 and len(Subset)%2!=0:
                        print('count is  ', count)
                        LWPC_list.loc[rep, 'Correct'] = correct_vals2[count]
                        count += 1
        # now assgin unbiased congruency conditions to the diagnostic items
        # again shuffle correct and congruence so left and right and congruent and incongruent 
        # are not (SLIGHTLY) biased across blocks and participants for the diagnostic items
        congruent_vals = self.Congruence*2
        random.shuffle(correct_vals)
        random.shuffle(congruent_vals)
        count = 0
        for block in set(LWPC_list.Block):
            for item in set(LWPC_diag):
                Subset = LWPC_list.index[(LWPC_list.Item == item) & (LWPC_list.Block == block)]
                for num, rep in enumerate(Subset):
                    # we 
                    if num < floor(self.diagnostic_rep/4):
                        LWPC_list.loc[rep, 'Correct'] = self.Correct[0]
                    if num >= floor(self.diagnostic_rep/4) and num < floor(self.diagnostic_rep/2):
                        LWPC_list.loc[rep, 'Correct'] = self.Correct[1]
                    if num%2 == 0:
                        LWPC_list.loc[rep, 'Congruency'] = self.Congruence[0]
                    if num%2 == 1:
                        LWPC_list.loc[rep, 'Congruency'] = self.Congruence[1]
                    if num == floor(self.diagnostic_rep/2):
                        # for the last index we draw randomly
                        LWPC_list.loc[rep, 'Congruency'] = congruent_vals[count]
                        LWPC_list.loc[rep, 'Correct'] = correct_vals[count]
                        count += 1
        # Next we return the datasets as individual blocks
        Block_A = LWPC_list[LWPC_list.Block == 1]
        Block_B = LWPC_list[LWPC_list.Block == 2]
        # next we first shuffle the dataset, reset the indices and set the ID 
        # variable to correspond to the new indexes (this variable is used by 
        #the pseudorandomization function)
        
        Block_A = Block_A.sample(frac=1)
        Block_A = Block_A.reset_index(drop=True)
        Block_A['ID'] = Block_A.index
        Block_B = Block_B.sample(frac=1)
        Block_B = Block_B.reset_index(drop=True)
        Block_B['ID'] = Block_B.index
        
        return Block_A, Block_B
    
    
    def create_IWPC_Blocks(self, IWSPC_feature):
        try:
            if IWSPC_feature is 'size':
                small_num = self.item_pool_dist_1[:2] + self.item_pool_dist_2[:2]
                large_num = self.item_pool_dist_1[2:] + self.item_pool_dist_2[2:]
            elif IWSPC_feature is 'distance':
                small_num = self.item_pool_dist_1
                large_num = self.item_pool_dist_2
        except NameError:
            print('No blockType has been specified. Please use either',
                                 'the string distance for to manipulate the ',
                                 'numerical distance between comparisons as a feature', 
                                 ' or size to manipulate comparisons with small',
                                 ' numbers (2-4) and large numbers (6-8) as a feature.')
        # Shuffle the comparisons before we randomly determine which items are
        #  inducer and which are diagnostic items
        random.shuffle(small_num)
        random.shuffle(large_num)
        # randomly determine which type of item will be manipualted as congruent
        # and which one as incongruent
        conditions = 'small', 'large'
        main_incon = random.sample(conditions,1)
        if main_incon == ['small']:
            main_con = ['large']
            s_minor = 'congruent'
            s_major = 'incongruent'
            l_minor = 'incongruent'
            l_major = 'congruent'
        else:
            main_con = ['small']
            s_minor = 'incongruent'
            s_major = 'congruent'
            l_minor = 'congruent'
            l_major = 'incongruent'
        # Assign the randomized items to diagnostic and inducer item lists
        IWPC_ind = (small_num[:2] + large_num[:2])*self.inducer_rep
        IWPC_diag = (small_num[2:] + large_num[2:])*self.diagnostic_rep
        Item_list =  IWPC_ind + IWPC_diag 
        # Item Type - Inducer Diagnostic Item
        Item_type = list(['Inducer'])*self.inducer_rep*int(self.item_number/2) + list(['Diagnostic'])*self.diagnostic_rep*int(self.item_number/2)
        # Split the sets in equally sized Blocks
        Block_Nr = [1]*int(self.inducer_nr/2) + [2]*int(self.inducer_nr/2) + [1]*int(self.diag_nr/2) + [2]*int(self.diag_nr/2)
        # Now throw together an empty data frame with columns correct response, 
        # congruence, ID and block
        IWPC_list = pd.DataFrame({'Item': Item_list, 'Trl_type': Item_type, 
                               'Correct': [None]*len(Item_list),
                                 'Congruency': [None]*len(Item_list),
                                 'Analysis_type': [None]*len(Item_list),
                                 'Block': Block_Nr,
                                 'ID': range(len(Item_list))})
        # now assgin biased congruency conditions to the inducer items
        # First shuffle correct so left and right are not biased across blocks
        '''
        Note: Ugly Code Part 2
        '''
        correct_vals = self.Correct*2
        correct_vals2 = self.Correct
        # First we add Proportions of congruent and incongruent items for each of
        # the inducer items with congruent and incongruent choices being balanced
        # by presentation side (Correct)
        for block in set(IWPC_list.Block):
            count = 0
            random.shuffle(correct_vals2)
            for ind_item in set(IWPC_ind):
                # get a subset of the data for each item
                Subset = IWPC_list.index[(IWPC_list.Item == ind_item) & (IWPC_list.Block == block)]
                # then add correct and congruency values to the dataframe
                for num, rep in enumerate(Subset):
                    if num%2 == 0:
                        IWPC_list.loc[rep, 'Correct'] = self.Correct[0]
                    if num%2 != 0:
                        IWPC_list.loc[rep, 'Correct'] = self.Correct[1]
                        # if we cannot assign correct and incorrect equally because
                        # an item is not presented an even number of times (must 
                        # happen for two items our chosen trial numbers). The last item is
                        # chosen randomly between correct and incorrect to make 
                        # sure the block is overall still balanced. Hard coded and
                        # ugly I know.
                    if num == len(Subset)-1 and len(Subset)%2!=0:
                        print('count is  ', count)
                        IWPC_list.loc[rep, 'Correct'] = correct_vals2[count]
                        count += 1
                    if ind_item in small_num:
                        if num%5 == 0:
                                IWPC_list.loc[rep, 'Congruency'] = s_minor
                        if num%5 != 0:
                                IWPC_list.loc[rep, 'Congruency'] = s_major
                        if main_con == ['small']: 
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_con'
                        else:
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_incon'
                    if ind_item in large_num:
                        if num%5 == 0:
                                IWPC_list.loc[rep, 'Congruency'] = l_minor
                        if num%5 != 0:
                                IWPC_list.loc[rep, 'Congruency'] = l_major
                        print(main_con)
                        if main_con == ['large']: 
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_con'
                        else:
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_incon'
        # now assgin unbiased congruency conditions to the diagnostic items
        # again shuffle correct and congruence so left and right and congruent and incongruent 
        # are not (SLIGHTLY) biased across blocks and participants for the diagnostic items
        congruent_vals = self.Congruence*2
        random.shuffle(correct_vals)
        random.shuffle(congruent_vals)
        count = 0
        for block in set(IWPC_list.Block):
            for item in set(IWPC_diag):
                Subset = IWPC_list.index[(IWPC_list.Item == item) & (IWPC_list.Block == block)]
                for num, rep in enumerate(Subset):
                    # we 
                    if num < floor(self.diagnostic_rep/4):
                        IWPC_list.loc[rep, 'Correct'] = self.Correct[0]
                    if num >= floor(self.diagnostic_rep/4) and num < floor(self.diagnostic_rep/2):
                        IWPC_list.loc[rep, 'Correct'] = self.Correct[1]
                    if num%2 == 0:
                        IWPC_list.loc[rep, 'Congruency'] = self.Congruence[0]
                    if num%2 == 1:
                        IWPC_list.loc[rep, 'Congruency'] = self.Congruence[1]
                    if num == floor(self.diagnostic_rep/2):
                        # for the last index we draw randomly
                        IWPC_list.loc[rep, 'Congruency'] = congruent_vals[count]
                        IWPC_list.loc[rep, 'Correct'] = correct_vals[count]
                        count += 1
                    print(main_con)
                    if item in small_num:
                        if main_con == ['small']: 
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_con'
                        else:
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_incon'
                    if item in large_num:
                        if main_con == ['large']: 
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_con'
                        else:
                            IWPC_list.loc[rep, 'Analysis_type'] = 'main_incon'
        # Next we return the datasets as individual blocks
        Block_A = IWPC_list[IWPC_list.Block == 1]
        Block_B = IWPC_list[IWPC_list.Block == 2]
        # next we shuffle the datasets, reset the indices and set the ID variable to correspond to the 
        # new indexes (this variable is used by the pseudorandomization function)
        Block_A = Block_A.sample(frac=1)
        Block_A = Block_A.reset_index(drop=True)
        Block_A['ID'] = Block_A.index
        Block_B = Block_B.sample(frac=1)
        Block_B = Block_B.reset_index(drop=True)
        Block_B['ID'] = Block_B.index
        # Since we need four blocks we simply copy the two we have already created
        Block_C = Block_A.copy()
        Block_C['Block'] = 3
        Block_C = Block_C.sample(frac=1)
        Block_C = Block_C.reset_index(drop=True)
        Block_C['ID'] = Block_C.index
        Block_D = Block_B.copy()
        Block_D['Block'] = 4
        Block_D = Block_D.sample(frac=1)
        Block_D = Block_D.reset_index(drop=True)
        Block_D['ID'] = Block_D.index
        
        return Block_A, Block_B, Block_C, Block_D
    
  
