# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 12:29:59 2020

@author: Julius Kricheldorff MSc.

@affilation: University of Oldenburg

@contact: julius@kricheldorff.de
"""


import random
import pandas as pd
import numpy as np
from time import time

class List_repair(object):
    def __init__(self, start_backtracking = 0.95, shuffle_limit = 50):
        self.backtrack_lim = start_backtracking
        self.shuff_lim = shuffle_limit

    def constrain_checks(self, target_list, start_ind):
        '''
        Function which returns a set of indices where each constrain was violated first.
    
        Parameters
        ----------
        target_list : Pandas data frame - 
            The data frame we would like to check.
        
        start_ind : integer -
            The index from which we would like to start our checks
        
        Returns
        -------
        violations : set - 
            Set of first indexes where a violation of the constrains 
            occured.
    
            '''
        # diagnostic items should not follow each other
        crit_diag = self.check_repeats(1, target_list['Trl_type'], item=['Diagnostic'], start_ind=start_ind)
        # the same item should not be shown more than twice in a row
        crit_item = self.check_repeats(2, target_list['Item'], item=None, start_ind=start_ind)
        # inducer items should not follow more than 4 times in a row
        crit_ind = self.check_repeats(4, target_list['Trl_type'], item=['Inducer'], start_ind=start_ind)
        # the same response side should not follow more than three times in a row
        crit_side = self.check_repeats(3, target_list['Correct'], item=None, start_ind=start_ind)
        # the same congruency should not follow more than three times in a row
        crit_con = self.check_repeats(3, target_list['Congruency'], item=None, start_ind=start_ind)
        # turns to one only if all our conditions are met
        crit_total = (crit_diag, crit_item, crit_ind, crit_side, crit_con)
        # return unique problematic indices
        violations = set(crit_total)
        return violations
        
    def check_repeats(self, numberRepeats, dframe, item, start_ind):
        '''
        Function outputs the index where a specified condition is not met.
        
        Parameters
        ----------
        numberRepeats : int
            Number of repeats that are allowed.
            
        dframe : pandas series
            Data series that we check for violations
        item : string or None
            check for repeats - None for any item - or a string with a specific 
            item
        start_ind : int
            index from which we start to search

        Returns
        -------
        int
            Returns either either the first index where a violdation occured or
            the length of the data series - as a reference that no violation was
            detected.
        '''
        # if the start index is smaller than the number of repeats we start at 
        # the number of repeats as our start index
        if start_ind <= numberRepeats+1:
            start = numberRepeats
        # else we select an earlier staring index (not sure why anymore) - but 
        # I guess I do not like making mistakes - better safe than sorry ;P. 
        else:
            start = start_ind - numberRepeats
        for i in range(start, len(dframe)):
            # get a subset including the preceed n (numberRepeats) indexes
            subset = dframe[i-numberRepeats:i+1]
            if item == None:
                # so if we just want to look at any generic number of reating
                # elements
                if len(subset.unique()) == 1:
                    # write down the index where the condition was violated. 
                    # We only return the last index (reps before are fine) 
                    indices =  i
                    return indices
            else:
                # if we want to check for one specific element
                if sum(subset.isin(item)) == numberRepeats+1:
                    indices = i
                    return indices
        # now return all unique indices that included violations as tuple
        return len(dframe) #tuple(set(indices))
    
   
     
    def repair(self, indices, frame):
        '''
        Function that takes the data frame and performs two kinds of repair 
        options in order to adhere to the constrains of the pseudo-randomization. 
        The first repair option is shuffling the list from the last index where 
        our constrains were not met. When shuffling improves the our list 
        (our defined constrains are met for more indexes) we move on the next 
        faulty index and continue in this fashion. Typically for the last indexes 
        we cannot improve by just using shuffling. For those indexes we use 
        backtracking - going through the list again and see if we can replace 
        the last faulty indices recursively with items earlier in the data frame 
        - without violating our constrains until the data frame adheres to all constrains. 
        If this fails the list reshuffled again and the process starts from scratch. 
        The algorithm implemented here is based on the idea described by Maarten 
        van Casteren and Mathew Davis. (See van Casteren, M., & Davis, M. H. 
        (2006). Mix, a program for pseudorandomization. Behavior research methods, 
        38(4), 584-589.). The code for the backtracking was adapted from a script 
        by Ray Toal (https://cs.lmu.edu/~ray/notes/backtracking/).

        Parameters
        ----------
        indices : Set of indices that violate each of our specified conditions.
        
        frame : The data frame we operate on.
        
        Returns
        -------
        frame : returns the now updated data frame

        '''
        # we first identify the first index where our restictions were not met
        min_1 = min(indices)
        # we reset the counter that keeps track of how often we shuffle 
        reset_counter = 0
        # var
        total = len(frame)
        while min_1 != len(frame):
             
             # next we randomly shuffle all rows from the part where our conditions aren't met yet
             frame.loc[min_1:] = np.random.permutation(frame.loc[min_1:])
             # we have to reset the indexes so that our search algorithm can find the correct locations 
             # of the current
             frame = frame.reset_index(drop=True)
             # check if we actually improved our list
             indices = self.constrain_checks(frame, min_1)
             try:
                 if min(indices) <= min_1:
                     reset_counter +=1
                     # if shuffling did not bring an improvement go back to the last state
                     frame = frame.sort_values('ID')
                     frame = frame.reset_index(drop=True)
                     # if we cannot improve by shuffling we resort to backtracking
                     if reset_counter > self.shuff_lim:
                         # backtracking is relatively slow - so we only use it if our data frame is sufficiently shuffled already
                         if round(min_1/total, 1) >= self.backtrack_lim:
                             # perform backtracking
                             solution = self.repair_reverse_backtrack(min_1, frame) 
                             # update the dataframe with the solution provided by backtracking
                             frame = self.update_backtracked_indices(frame , solution)
                             # perform the condition check again
                             indices = self.constrain_checks(frame, 0)
                             frame['ID'] = frame.index
                             min_1 = min(indices)
                         else:
                             # if we get stuck early restart
                             reset_counter = 0
                             min_1 = 0
                         
                 else:
                     reset_counter = 0
                     frame['ID'] = frame.index
                     min_1 = min(indices)
             except(ValueError):
                  print('randomization successfully!')
        return frame
    
    
    def repair_reverse_backtrack(self, position, frame):
        '''
        Code is adapted and was originally written by Ray Toal in his CSMI 282 course - see:
            https://cs.lmu.edu/~ray/notes/backtracking/. Implements backtracking 
            for the last few indices that could not be properly assigned 
            by the random shuffling solution. 

        Parameters
        ----------
        position : Index where shuffling failed.
        
        frame : The data frame we operate on.
        
        Returns : returns the solution of the backtracking function, which is a 
        list of list indices, the first position indicating the index that is 
        to be relaced, the second being the suitable replacement.
        -------
        TYPE
            DESCRIPTION.

        '''
        print('called again, position is ', position)
        # create a list were we will safe the swaps we will be doing (in order)
        solution =  list()
        def backtrck(best_pos,frame):
            '''
            

            Parameters
            ----------
            best_pos : The current index until which our assignmend conditions are fullfilled.
    
            frame : The data frame we operate on.

            Returns : The above described solution list. 
            -------
            TYPE
                DESCRIPTION.

            '''
            # try all available items in the dataframe
            for i in range(0,len(frame)): 
                # swap the olde and the new index
                updatedFrame = self.swap_element(frame, best_pos, i)
                # determine which items violate our conditions (on the whole list)
                index = self.constrain_checks(updatedFrame, 0)
                # determine the smallest index a.k.a. hopefully the new best solution 
                new_best = min(index)
                # check if we already reached our criteria or our new best 
                # solution is better than the old one
                if new_best >= len(frame) or new_best > best_pos:
                    # safe the new best solution
                    solution.append([best_pos, i])
                    # end the the function when done or proceed recursively
                    if new_best >= len(frame) or backtrck(new_best, updatedFrame):
                        # if all worked out fine we return the updated frame and
                        # the final index
                        return solution
                # if we fail to solve the problem we start over
            return False  
        return backtrck(position, frame)

    def swap_element(self, frame, problem_pos, swap_pos):
        '''
        Function that swaps indices in the data frame

        Parameters
        ----------
        frame : Data frame we operate on
        
        problem_pos : The first index for the swap (the index we want to replace)
        
        swap_pos : The second index for the swap (the candidate index for replacement)

        Returns
        -------
        data : The data frame with the indexes swapped.

        '''
        #get a copy
        data = frame.copy()
        # save the two items 
        old, new = data.iloc[problem_pos].copy(), data.iloc[swap_pos].copy()
        #replace the two items with each other
        data.iloc[problem_pos], data.iloc[swap_pos] = new, old
        # reset indexes
        data = data.reset_index(drop=True)
        return data
    
    def update_backtracked_indices(self, frame, ind_list):
        '''
        Function that updates the returned indices from the backtracking function

        Parameters
        ----------
        frame : The data frame we performed backtracking on
        ind_list : The solution that the backtracking function provided

        Returns
        -------
        final_frame : The updated datafram

        '''
        final_frame = frame
        for pair in ind_list:
            final_frame = self.swap_element(final_frame, pair[0], pair[1])
        return final_frame

## Example script

# Get the dataframe
randomize = BlockListRandomizer()
my_list = randomize.create_LWPC_sets()
my_list = my_list.sample(frac=1)
my_list = my_list.reset_index(drop=True)
target_list = my_list

# import repair functions
repair = List_repair()
indices = repair.constrain_checks(my_list, 0)
min_1 = min(indices)

# repair the data frame
fina_1 = repair.repair(indices,  my_list)
repair.constrain_checks(fina_1 , 0)