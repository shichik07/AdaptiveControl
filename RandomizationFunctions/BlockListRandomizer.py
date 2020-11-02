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

class BlockListRandomizer(object):
    def __init__(self, Trial_Nr = 336, Block_Nr = 4, Ind_prop = 0.8,
                 Inducer_Nr = 236, IWSPC_feature = 'distance', size_L = 300,
                 size_M = 260):
        self.trlnum = Trial_Nr
        self.Block_Nr = Block_Nr
        if IWSPC_feature == 'distance':
            # the manipulated features will be either distances one and two
            self.IWSPC = IWSPC_feature 
        elif IWSPC_feature == 'size':
            # the manipulated feature will be low and high numbers respectively
            self.IWSPC = IWSPC_feature 
        self.IWSPC = IWSPC_feature 
        self.Order = [['MI', 'MC'], ['IWPC_1', 'IWPC_2']]
        random.shuffle(self.Order)
        # self.item_pool_dist_1 = [[1,2], [3,4], [6,7], [8,9]]
        # self.item_pool_dist_2 = [[1,3], [2,4], [6,8], [7,9]]
        self.item_pool_dist_1 = ['A1', 'B1', 'C1', 'D1']
        self.item_pool_dist_2 = ['A2', 'B2', 'C2', 'D2']
        self.item_number =len(self.item_pool_dist_1)+len(self.item_pool_dist_2)
        self.inducer_rep = int(Inducer_Nr/(self.item_number/2))
        self.diagnostic_rep = int((Trial_Nr-Inducer_Nr)/(self.item_number/2))
       
    def create_LWPC_sets(self):
        # function that creates the basic LWPC List layout
        D1 = list(self.item_pool_dist_1)
        D2 = list(self.item_pool_dist_2)
        random.shuffle(D1)
        random.shuffle(D2)
        # Randomly 
        LWPC_ind = (D1[:2] + D2[:2])*self.inducer_rep
        LWPC_diag = (D1[2:] + D2[2:])*self.diagnostic_rep
        Item_list =  LWPC_ind + LWPC_diag 
        Item_type = list(['Inducer'])*self.inducer_rep*4 + list(['Diagnostic'])*self.diagnostic_rep*4
        # add equal Nr of trials that are congruent, incongruent * Left or right sided
        # answers - first for the inducer items (surely this can be done simpler)
        if self.inducer_rep%2 != 0:
            # if we cannot balance the responses completly across sets - simply
            # allocate the last three randomly (while maintaining responses and
            # congruency balanced)
            add_side = ['Left','Right']*2
            add_congr = ['Congruent','Incongruent']*2
            random.shuffle(add_side)
            random.shuffle(add_congr)
            CR_ind = ((['Left'])*4 + (['Right'])*4)*int(self.inducer_rep//2) + add_side
            Congr_ind = (['Congruent'])*int(self.inducer_rep//2)*4  + (['Incongruent'])*int(self.inducer_rep//2)*4  + add_congr
        else:
            # in case items are balanced
            CR_ind = ((['Left'])*4 + (['Right'])*4) * int(self.inducer_rep//2)
            Congr_ind = (['Congruent'])*int(self.inducer_rep//2)*4  + (['Incongruent'])*int(self.inducer_rep//2)*4 
        # Repeat the same exercise now for the diagnostic items
        if self.diagnostic_rep%2 != 0:
            # if we cannot balance the responses completly across sets - simply
            # allocate the last three randomly (while maintaining responses and
            # congruency balanced)
            add_side = ['Left','Right']*2
            add_congr = ['Congruent','Incongruent']*2
            random.shuffle(add_side)
            random.shuffle(add_congr)
            CR_diag = ((['Left'])*4 + (['Right'])*4)*int(self.diagnostic_rep//2) + add_side
            Congr_diag = (['Congruent'])*int(self.diagnostic_rep//2)*4  + (['Incongruent'])*int(self.diagnostic_rep//2)*4  + add_congr
        else:
            # in case items are balanced
            CR_diag = ((['Left'])*4 + (['Right'])*4) * int(self.inducer_rep//2)
            Congr_diag = (['Congruent'])*int(self.inducer_rep//2)*4  + (['Incongruent'])*int(self.inducer_rep//2)*4 
        Congruence = Congr_ind + Congr_diag
        Correct_Resp = CR_ind + CR_diag
        # Construct a Pandas data frame from it
        LWPC_list = pd.DataFrame({'Item': Item_list, 'Trl_type': Item_type, 
                               'Correct': Correct_Resp,
                                 'Congruency': Congruence,
                                 'ID': range(len(Item_list))})
        # shuffle the rows until all of our conditions are met 
        # i = 0
        # while constrain_checks(LWPC_list) == 0:
        #     LWPC_list = LWPC_list.sample(frac=1)
        return LWPC_list
    constrain_checks
   
    def constrain_checks(target_list, start_ind):
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
        crit_diag = check_repeats(1, target_list['Trl_type'], item=['Diagnostic'], start_ind=start_ind)
        # the same item should not be shown more than twice in a row
        crit_item = check_repeats(2, target_list['Item'], item=None, start_ind=start_ind)
        # inducer items should not follow more than 4 times in a row
        crit_ind = check_repeats(4, target_list['Trl_type'], item=['Inducer'], start_ind=start_ind)
        # the same response side should not follow more than three times in a row
        crit_side = check_repeats(3, target_list['Correct'], item=None, start_ind=start_ind)
        # the same congruency should not follow more than three times in a row
        crit_con = check_repeats(3, target_list['Congruency'], item=None, start_ind=start_ind)
        # turns to one only if all our conditions are met
        crit_total = (crit_diag, crit_item, crit_ind, crit_side, crit_con)
        # return unique problematic indices
        violations = set(crit_total)
        return violations
        
    def check_repeats(numberRepeats, dframe, item, start_ind):
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
        return len(frame) #tuple(set(indices))
    
   
     
    def repair(indices, frame):
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
       
        while min_1 != len(frame):
             
             # next we randomly shuffle all rows from the part where our conditions aren't met yet
             frame.loc[min_1:] = np.random.permutation(frame.loc[min_1:])
             # we have to reset the indexes so that our search algorithm can find the correct locations 
             # of the current
             frame = frame.reset_index(drop=True)
             # check if we actually improved our list
             indices = constrain_checks(frame, min_1)
             try:
                 if min(indices) <= min_1:
                     reset_counter +=1
                     # if shuffling did not bring an improvement go back to the last state
                     frame = frame.sort_values('ID')
                     frame = frame.reset_index(drop=True)
                     # if we cannot improve by shuffling we resort to backtracking
                     if reset_counter > 50:
                         # backtracking is relatively slow - so we only use it if our data frame is sufficiently shuffled already
                         if round(min_1/total, 1) >= 95:
                             # perform backtracking
                             solution = repair_reverse_backtrack(min_1, frame) 
                             # update the dataframe with the solution provided by backtracking
                             frame = update_backtracked_indices(frame , solution)
                             # perform the condition check again
                             indices = constrain_checks(frame, 0)
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
    
    
    def repair_reverse_backtrack(position, frame):
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
                updatedFrame = swap_element(frame, best_pos, i)
                # determine which items violate our conditions (on the whole list)
                index = constrain_checks(updatedFrame, 0)
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

    def swap_element(frame, problem_pos, swap_pos):
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
    
    def update_backtracked_indices(frame, ind_list):
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
            final_frame = swap_element(final_frame, pair[0], pair[1])
        return final_frame
        
             
    def repair_backtrack(curent_min, frame):
        criterion = 0
        repeats = 0
        # since we do not get any further by shuffling we will try to swap 
        # random rows from prior to the sticking point
        
        while criterion == 0:
            repeats += 1
            print('Call Backtrack', repeats, 'times.')
            lucky_nr = random.choice(range(curent_min))
            old, new = frame.iloc[curent_min].copy(), frame.iloc[lucky_nr].copy()
            frame.iloc[curent_min], frame.iloc[lucky_nr] = new, old
            # next check if we actually made things better or worse
            # reset indexes
            frame = frame.reset_index(drop=True)
            # check if we actually improved our list from the beginning
            indices = constrain_checks(frame, 0)
            try:
                #check if our criterion has already been reached
                if min(indices) <= curent_min:
                    # if shuffling did not bring an improvement go back to the last state
                    frame = frame.sort_values('ID')
                    frame = frame.reset_index(drop=True)
                    # and start again
                    if repeats == 50:
                        print('End backtrack.')
                        # if we seemingly cannot find a candidate replacement the whole 
                        # process starts anew
                        return frame, curent_min -50
                else:
                    criterion = 1
                    frame['ID'] = frame.index
                    min_1 = min(indices)
                    return frame, min_1 
            except(ValueError):
                print('The list has successfully been randomized!')
                return frame, None
    
    def repair_backtrack1(frame, min_1):
        """Finds a solution to a backtracking problem.
        
        Code is adapted and was originally written by Ray Toal in his CSMI 282 course - see:
            https://cs.lmu.edu/~ray/notes/backtracking/
    
        values     -- a sequence of values to try, in order. For a map coloring
                      problem, this may be a list of colors, such as ['red',
                      'green', 'yellow', 'purple']
        safe_up_to -- a function with two arguments, solution and position, that
                      returns whether the values assigned to slots 0..pos in
                      the solution list, satisfy the problem constraints.
        size       -- the total number of “slots” you are trying to fill
        frame      -- the data frame that we want to work with
        min_1      -- the current smallest instance until which the constrains
                      data frame work
    
        Return the solution as a list of values.
        """
        size = len(frame) - min_1
        solution = [None] * size
        values = frame['ID'].loc[min_1:]
        print('Call repair_backtracking.')
        def extend_solution(position, values):
            print('We are at index ', position, ' and we have a total of ', len(values), ' to go.')
            for value in values:
                solution[position] = value
                old, new = frame.iloc[min_1 + position].copy(), frame[frame['ID'] == value].copy()
                frame.iloc[min_1 + position], frame.iloc[new.index[0]] = new.squeeze(), old
                indices = constrain_checks(frame, min_1 + position)
                try: 
                    if min(indices) >= min_1 + position:
                        print('added the value', value)
                        # get the index of the value we want to drop from the list
                        value_idx = values[values == value].index[0]
                        # drop the value from the list
                        new_values = values.drop(value_idx)
                        # while there are still spots to fill, pass the new list to the function
                        if position >= size-1 or extend_solution(position+1, new_values):
                            return solution, frame
                except(ValueError):
                    return solution, frame
            print('Take a step back?')
            return None, frame

        return extend_solution(0, values)
  
# Playing_field
randomize = BlockListRandomizer()
my_list = randomize.create_LWPC_sets()
my_list = my_list.sample(frac=1)
my_list = my_list.reset_index(drop=True)
target_list = my_list
indices = constrain_checks(my_list, 0)
min_1 = min(indices)
frame = my_list
frame.loc[min_1:] = np.random.permutation(frame.loc[min_1:])
# we have to reset the indexes so that our search algorithm can find the correct locations 
# of the current
frame['ID'] = frame.index
indices = constrain_checks(frame, 0)
min_1 = min(indices)


b = np.random.permutation(frame.loc[min_1:])
a = frame.loc[min_1:]
c = b.loc[3:]
my_list.iloc[3:7] = my_list.iloc[133:137].sample(frac=1)
my_list.loc[3:7] = np.random.permutation(my_list.loc[233:237])


randomize = BlockListRandomizer()
my_list = randomize.create_LWPC_sets()
my_list = my_list.sample(frac=1)
my_list = my_list.reset_index(drop=True)
dm = convert.from_pandas(my_list)
print(dm)


# specify constraints
ef = Enforce(dm)
ef.add_constraint(MaxRep, cols=[dm.Item], maxrep=1)
ef.add_constraint(MaxRep, cols=[dm.Congruency], maxrep=2)
ef.add_constraint(MaxRep, cols=[dm.Correct], maxrep=2)
ef.add_constraint(MaxRep, cols=[dm.Trl_type], maxrep=2)
# Enforce the constraints
dm = ef.enforce()
# See the resulting DataFrame and a report of how long the enforcement took.
print(dm)
print(ef.report)
## Working Backtracker
import random

vals = [1,2,2,2,4,5]*60
random.shuffle(vals)
res = [None]* len(vals)
total = len(vals)
emergency = 0
    
def solve(values, results, position, total):
    print(results)
    #print('the length of all remaining vlaues is', len(values))
    #print(position)
    print('The length of the results is', len(results))
    
    for ind, integer in enumerate(values):
        results[position] = integer
        if results[position-1] != results[position]:
            if position >= total-1 or solve(remove_element(values, ind), results, position+1, total):
                return results
   
    return None

def remove_element(t_list, ind):
        new_list = t_list[:ind] + t_list[ind+1:]
        return new_list                         

    
final =  solve(vals, res, 0, total) 
