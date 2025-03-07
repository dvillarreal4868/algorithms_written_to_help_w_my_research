#! /usr/bin/env python

# This python script contains a function where a takes an individual RNA sequence
# It then performs a either a dinucleotide and codon or only codon preserving shuffle of the sequence
# The codon_swap and codon_dinuc_swap code was originally written by Ofer Kimchi

# Change occurances to map_occurances
# translation needs to be better written

##################################################################################################################################

# Imports
import numpy as np
import copy
from time import time
from collections import Counter
import random
import coding_notebook

##################################################################################################################################

def generate_codon_table():
    '''
    Generate a codon table.
    Type: dictionary.
    
    '''
    
    codon_table = dict()  # for each aa, what codons code for it?
    codon_table['F'] = ['UUU', 'UUC']
    codon_table['L'] = ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG']
    codon_table['I'] = ['AUU', 'AUC', 'AUA']
    codon_table['M'] = ['AUG']
    codon_table['V'] = ['GUU', 'GUC', 'GUA', 'GUG']
    codon_table['S'] = ['UCU', 'UCC', 'UCA', 'UCG']
    codon_table['P'] = ['CCU', 'CCC', 'CCA', 'CCG']
    codon_table['T'] = ['ACU', 'ACC', 'ACA', 'ACG']
    codon_table['A'] = ['GCU', 'GCC', 'GCA', 'GCG']
    codon_table['Y'] = ['UAU', 'UAC']
    codon_table['X'] = ['UAA', 'UAG', 'UGA']  # stop
    codon_table['H'] = ['CAU', 'CAC']
    codon_table['Q'] = ['CAA', 'CAG']
    codon_table['N'] = ['AAU', 'AAC']
    codon_table['K'] = ['AAA', 'AAG']
    codon_table['D'] = ['GAU', 'GAC']
    codon_table['E'] = ['GAA', 'GAG']
    codon_table['C'] = ['UGU', 'UGC']
    codon_table['W'] = ['UGG']
    codon_table['R'] = ['CGU', 'CGC', 'CGA', 'CGG']
    codon_table['S'] += ['AGU', 'AGC']
    codon_table['R'] += ['AGA', 'AGG']
    codon_table['G'] = ['GGU', 'GGC', 'GGA', 'GGG']
    
    return codon_table

##################################################################################################################################

def generate_reverse_codon_table():
    '''
    Generate reverse codon table.
    Type: Dictionary 
    '''
    codon_table = generate_codon_table()
    reverse_codon_table = dict()  # for each codon, what aa does it code for?
    num_codons_counter = 0
    for amino_acid in codon_table.keys():
        for codon in codon_table[amino_acid]:
            num_codons_counter += 1
            reverse_codon_table[codon] = amino_acid
            
    return reverse_codon_table

##################################################################################################################################

def codon_swap(sequence, num_random_sequences=10):
    '''
    Perform a shuffle on to your sequence.
    Returns a list of strings.
    
    import numpy as np
    import copy
    from time import time
    
    '''
    
    len_RNA = len(sequence)
    num_codons = len_RNA // 3
    
    codon_table = generate_codon_table()
    reverse_codon_table = generate_reverse_codon_table()

    nts = 'ACGU'
    amino_acids = [aa for aa in codon_table.keys()]
    
    # Convert the RNA sequence to the protein it codes for
    aa_sequence = '' 
    for codon_counter in range(num_codons):
        codon = sequence[3 * codon_counter : 3 * (codon_counter + 1)]
        aa_sequence += reverse_codon_table[codon]
        
    # For each codon, how many times does it appear in the sequence?
    codons_in_sequence = dict()  

    # First initialize the table
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                codons_in_sequence[nt1 + nt2 + nt3] = 0

    # Then populate it
    for codon_counter in range(num_codons):
        codon = sequence[3 * codon_counter : 3 * (codon_counter + 1)]
        codons_in_sequence[codon] += 1

    # For each amino acid, make a list of all codons (with repeats) that the
    # sequence has coding for that amino acid

    codons_in_sequence_per_aa = dict()

    # Initialize the table
    for aa in amino_acids:
        codons_in_sequence_per_aa[aa] = []

    # Populate table
    for codon_counter in range(num_codons):
        codon = sequence[3 * codon_counter : 3 * (codon_counter + 1)]
        aa = reverse_codon_table[codon]
        codons_in_sequence_per_aa[aa] += [codon]

    # Create random RNA coding for the same protein with same codon usage
    random_sequences = []
    for _ in range(num_random_sequences):
        random_sequence = ''

        # Shuffle the order of codons for each aa to create new random sequence
        codons_in_sequence_per_aa_copy = copy.deepcopy(codons_in_sequence_per_aa)
        for aa in amino_acids:
            np.random.shuffle(codons_in_sequence_per_aa_copy[aa])

      # Generate the random sequence
        for codon_counter in range(num_codons):
            # for each aa, take a random codon from the list of codons used in the 
            # sequence to code for that aa, and remove that codon from the list.
            aa = aa_sequence[codon_counter] 
            random_sequence += codons_in_sequence_per_aa_copy[aa].pop()

        random_sequences += [random_sequence]

    return random_sequences

##################################################################################################################################

def occurrences(string, sub):
    '''
    Function for indexing overlapping occurrences of a substring
    
    '''
    count = start = 0
    lst = []
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
            lst.append(start)
        else:
            return lst
        
##################################################################################################################################

def translation(coding_sequence):
    '''
    Translate your RNA sequence
    '''
    len_RNA = len(coding_sequence)
    num_codons = len_RNA // 3
    codon_table = generate_codon_table()
    reverse_codon_table = generate_reverse_codon_table()
    
    aa_sequence = ''  # Convert the RNA sequence to the protein it codes for

    for codon_counter in range(num_codons):
        codon = coding_sequence[3 * codon_counter : 3 * (codon_counter + 1)]
        aa_sequence += reverse_codon_table[codon]
        
    return aa_sequence

##################################################################################################################################

def codon_pairwise_swap(sequence, num_random_sequences=10, num_swaps=5):
    '''
    Perform a pairwise shuffle on to your sequence.
    Returns a list of strings.
    
    import numpy as np
    import copy
    from time import time
    
    '''

    RNA_sequence = copy.deepcopy(sequence)

    # Start of Ofer's code
    len_RNA = len(sequence)
    num_codons = len_RNA // 3

    codon_table = generate_codon_table()
    reverse_codon_table = generate_reverse_codon_table()

    nts = 'ACGU'
    amino_acids = [aa for aa in codon_table.keys()]

    # Convert the RNA sequence to the protein it codes for
    aa_sequence = ''
    for codon_counter in range(num_codons):
        codon = sequence[3 * codon_counter : 3 * (codon_counter + 1)]
        aa_sequence += reverse_codon_table[codon]

    # End of Ofer's code

    # Make a dictionary where the keys are the amino acids and the items are the list of map_occurances
    # Do not add amino acid W or M and do not add any amino acid that shows up less than twice
    # Do not add list of codons that only contains one version
    aa_occurance_dict = {aa : list(np.subtract(occurrences(aa_sequence, aa), 1)) for aa in codon_table.keys()\
                         if (len(list(np.subtract(occurrences(aa_sequence, aa), 1))) > 1) and\
                         (aa != 'W') and (aa != 'M')}
    temp_dict = {key : [sequence[3 * i : 3 * (i + 1)] for i in lst if len(lst) != Counter(lst)[lst[0]]]\
                 for key, lst in aa_occurance_dict.items()}  
    aa_codon_dict = {key : lst for key, lst in temp_dict.items() if len(lst) != Counter(lst)[lst[0]]}

    sequences = []
    while len(sequences) < num_random_sequences:
        # Perform some number of pairwise swaps
        cntr = 0
        while cntr < num_swaps:
            # Choose a random codon and count how many times it appears in the sequence
            aa_to_change = random.choice(list(aa_codon_dict.keys())) # Choosing a random amino acid to change
            lst_of_occurances = aa_occurance_dict[aa_to_change]

            # Make a list of codons based off occurances
            lst_of_codons = [sequence[3 * i : 3 * (i + 1)] for i in lst_of_occurances]

            # Perform pairwise swap of codons
            i1, i2 = random.sample(lst_of_occurances, 2)

            if i1 > i2: # Orientate in correct order
                i1, i2 = i2, i1

            codon_1 = sequence[3 * i1 : 3 * (i1 + 1)]
            codon_2 = sequence[3 * i2 : 3 * (i2 + 1)]

            # If codon_1 equals codon_2, then choose two more codons
            while codon_1 == codon_2:     
                i1, i2 = random.sample(lst_of_occurances, 2)

                if i1 > i2: # Orientate in correct order
                    i1, i2 = i2, i1    

                codon_1 = sequence[3 * i1 : 3 * (i1 + 1)]
                codon_2 = sequence[3 * i2 : 3 * (i2 + 1)]


            sequence = sequence[0 : 3 * i1] + codon_2 +\
            sequence[3 * (i1 + 1) : 3 * i2] + codon_1 + sequence[3 * (i2 + 1) : len(sequence)]

            cntr += 1
        sequences.append(sequence)
        sequence = RNA_sequence

    return sequences

##################################################################################################################################

def get_codons_in_sequence_dict(RNA_sequence_input):
    '''
    For each codon, how many times does it appear in the sequence? 

    '''
    len_RNA = len(RNA_sequence_input)
    num_codons = len_RNA // 3
    nts = 'ACGU'

    codons_in_sequence = dict()  

    # First initialize the table
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                codons_in_sequence[nt1 + nt2 + nt3] = 0

    # Then populate it
    for codon_counter in range(num_codons):
        codon = RNA_sequence_input[3 * codon_counter : 3 * (codon_counter + 1)]
        codons_in_sequence[codon] += 1

    return(codons_in_sequence)

##################################################################################################################################

def get_dinucleotides_in_sequence_dict(RNA_sequence_input):
    '''
    For each dinucleotide, how many times does it appear in the sequence?
    Only need to consider dinucleotides bridging codons, since 
    if codons are conserved, all other dinucleotides will be also

    '''
    len_RNA = len(RNA_sequence_input)
    num_codons = len_RNA // 3
    nts = 'ACGU'

    dinucleotides_in_sequence = dict()  

  # First initialize the table
    for nt1 in nts:
        for nt2 in nts:
            dinucleotides_in_sequence[nt1 + nt2] = 0

  # Then populate it
    for dinucleotide_counter in range(num_codons - 1):
        dinucleotide = RNA_sequence_input[3 * dinucleotide_counter + 2 : 3 * (dinucleotide_counter + 1) + 1]
        dinucleotides_in_sequence[dinucleotide] += 1

    return(dinucleotides_in_sequence)

##################################################################################################################################

def codon_dinuc_swap(RNA, num_random_sequences=10, num_swaps=5):
    '''
    Perform a pairwise shuffle on to your sequence.
    Returns a list of strings.

    '''
    
    codon_table = generate_codon_table()
    reverse_codon_table = generate_reverse_codon_table()

    RNA_sequence = RNA

    reverse_codon_table = dict()  # for each codon, what aa does it code for?
    num_codons_counter = 0
    for amino_acid in codon_table.keys():
      for codon in codon_table[amino_acid]:
        num_codons_counter += 1
        reverse_codon_table[codon] = amino_acid

    nts = 'ACGU'
    nts_array = ['A', 'C', 'G', 'U']
    amino_acids = [aa for aa in codon_table.keys()]


    len_RNA = len(RNA_sequence)
    num_codons = len_RNA // 3

    # Convert the RNA sequence to the protein it codes for
    aa_sequence = ''
    for codon_counter in range(num_codons):
      codon = RNA_sequence[3 * codon_counter : 3 * (codon_counter + 1)]
      aa_sequence += reverse_codon_table[codon]


    codons_in_sequence = get_codons_in_sequence_dict(RNA_sequence)

    # For each amino acid, make a list of all codons (with repeats) that the
    # sequence has coding for that amino acid

    codons_in_sequence_per_aa = dict()

    # Initialize the table
    for aa in amino_acids:
      codons_in_sequence_per_aa[aa] = []

    # Populate table
    for codon_counter in range(num_codons):
      codon = RNA_sequence[3 * codon_counter : 3 * (codon_counter + 1)]
      aa = reverse_codon_table[codon]
      codons_in_sequence_per_aa[aa] += [codon]

    # Print table
    # codons_in_sequence_per_aa

    dinucleotides_in_sequence = get_dinucleotides_in_sequence_dict(RNA_sequence)

    # Keep codon usage and dinucleotides constant by randomly switching codons that are followed by the same first nt

    # For each aa, for each nt it can start with, and for each subsequent nt (i.e. 
    # first nt of following codon), keep track of the places where each combination
    # of these appears. These codons can then be swapped straightforwardly while 
    # maintaining constant codon & dinucleotide contents.

    codons_with_same_subsequent_nt_per_aa = dict()

    for e, aa in enumerate(aa_sequence[:-1]):
      next_nt = RNA_sequence[3 * (e + 1)]
      first_nt = RNA_sequence[3 * e]
      if aa + first_nt + next_nt in codons_with_same_subsequent_nt_per_aa.keys():
        codons_with_same_subsequent_nt_per_aa[aa + first_nt + next_nt] += [e]
      else:
        codons_with_same_subsequent_nt_per_aa[aa + first_nt + next_nt] = [e]

    # codons_with_same_subsequent_nt_per_aa  # prints out the table

    random_sequences_swap = []  # an empty list to store the random sequences we generate

    # keep track of which codons were swapped for each sequence for testing purposes
    # swapped_codons = [[] for _ in range(num_random_sequences)]   # comment out when not testing

    for rand_sequence_counter in range(num_random_sequences):
      random_sequence = list(copy.deepcopy(RNA_sequence))  # since we need to do assignment to string

      for e in range(num_swaps):
        # Pick a random codon
        random_codon_pos = np.random.randint(1, len_RNA//3 - 1)

        # Find its corresponding key in codons_with_same_subsequent_nt_per_aa
        for afl, pos_list in codons_with_same_subsequent_nt_per_aa.items():
          # Loop through each dictionary element to find the one whose value is a
          # list that includes random_codon_pos
          if random_codon_pos in pos_list:
              afl_of_random_codon = afl  # afl = aa + first_nt + last_nt is the key
              break  # don't need to continue looping after we've found it

        # If there are other codons with this key, randomly switch it with one of those
        if len(codons_with_same_subsequent_nt_per_aa[afl_of_random_codon]) > 1:
          # Otherwise, there's nothing to swap it with

          # Get the position to swap it with      
          codon_position_to_swap = np.random.choice(
              codons_with_same_subsequent_nt_per_aa[afl_of_random_codon])

          # Make the swap
          random_codon = copy.copy(
              random_sequence[3 * random_codon_pos : 3 * (random_codon_pos + 1)])
          codon_to_swap = copy.copy(
              random_sequence[3 * codon_position_to_swap : 3 * (codon_position_to_swap + 1)])
          random_sequence[3 * random_codon_pos : 3 * (random_codon_pos + 1)] = codon_to_swap
          random_sequence[3 * codon_position_to_swap : 3 * (codon_position_to_swap + 1)] = random_codon

          # Comment out the following line when not testing, to save time & space
          # swapped_codons[rand_sequence_counter] += [(random_codon_pos, codon_position_to_swap)]
      random_sequences_swap += [''.join(random_sequence)]

    # Check that codon usage and dinucleotide content are unchanged:
    for e, random_sequence in enumerate(random_sequences_swap):
      if not (codons_in_sequence == 
            get_codons_in_sequence_dict(random_sequence)):
        print('We have a codon usage error in sequence # ' + str(e))

      if not (dinucleotides_in_sequence == 
            get_dinucleotides_in_sequence_dict(random_sequence)):
        print('We have a dinucleotide usage error in sequence # ' + str(e))

    return random_sequences_swap

##################################################################################################################################

def nt_pairwise_swap(sequence, num_random_sequences=10, num_swaps=5):
    '''
    Perform a pairwise shuffle on to your sequence.
    Returns a list of strings.
    
    import numpy as np
    import copy
    from time import time
    
    '''
    nts = 'ACGU'

    RNA_sequence = copy.deepcopy(sequence)

    len_RNA = len(sequence)

    lst = list(np.arange(len_RNA))

    sequences = []
    while len(sequences) < num_random_sequences:
        # Perform some number of pairwise swaps
        for cntr in range(num_swaps):
            # Randomly choose two nucleotides to swap 
            i1, i2 = random.sample(lst, 2)

            nt_1 = sequence[i1]
            nt_2 = sequence[i2]
            # If nt_1 equals nt_2, then choose two more codons
            if nt_1 == nt_2:
                while nt_1 == nt_2:
                    i1, i2 = random.sample(lst, 2)
                    nt_1 = sequence[i1]
                    nt_2 = sequence[i2]
            temp = list(sequence)
            temp[i1] = nt_2
            temp[i2] = nt_1
            sequence = ''.join(temp)


        sequences.append(sequence)
        sequence = copy.deepcopy(RNA_sequence)

    return sequences

##################################################################################################################################

def dinuc_pairwise_swap(RNA_sequence, num_random_sequences=10, num_swaps=5):
    '''
    Perform pairwise dinucleotide swapping
    '''
    nts = sorted(set('ACGU'))
    RNA_sequence_dinuc = coding_notebook.dinuc_count(RNA_sequence)
    
    lst_of_sequences = []

    for dummy_var_1 in range(num_random_sequences):

        sequence = RNA_sequence

        for dummy_var_2 in range(num_swaps):

            # Grab 4 nts at random
            sub_seq_1 = ''.join([random.sample(nts, 1)[0] for _ in range(4)])
            dinuc_pair_1 = sub_seq_1[1:3] # Dinucleotide Pair to swap

            # Define which Dinucleotide pair to swap with
            dinuc_pair_2 = ''.join([random.sample(nts, 1)[0] for _ in range(2)])

            # Do not have dinucleotides match
            if dinuc_pair_1 == dinuc_pair_2:
                while dinuc_pair_1 == dinuc_pair_2:
                    dinuc_pair_2 = ''.join([random.sample(nts, 1)[0] for _ in range(2)])


            sub_seq_2 = sub_seq_1[0] + dinuc_pair_2 + sub_seq_1[3]

            # Find all the occurances for both sub_sequences
            occurances_1 = coding_notebook.map_of_occurrences(sequence, sub_seq_1)
            occurances_2 = coding_notebook.map_of_occurrences(sequence, sub_seq_2)

            # If occurances list is empty
            if (len(occurances_1) == 0) | (len(occurances_2) == 0):
                while (len(occurances_1) == 0) | (len(occurances_2) == 0):
                    # Grab 4 nts at random
                    sub_seq_1 = ''.join([random.sample(nts, 1)[0] for _ in range(4)])
                    dinuc_pair_1 = sub_seq_1[1:3] # Dinucleotide Pair to swap

                    # Define which Dinucleotide pair to swap with
                    dinuc_pair_2 = ''.join([random.sample(nts, 1)[0] for _ in range(2)])

                    # Do not have dinucleotides match
                    if dinuc_pair_1 == dinuc_pair_2:
                        while dinuc_pair_1 == dinuc_pair_2:
                            dinuc_pair_2 = ''.join([random.sample(nts, 1)[0] for _ in range(2)])


                    sub_seq_2 = sub_seq_1[0] + dinuc_pair_2 + sub_seq_1[3]

                    # Find all the occurances for both sub_sequences
                    occurances_1 = coding_notebook.map_of_occurrences(sequence, sub_seq_1)
                    occurances_2 = coding_notebook.map_of_occurrences(sequence, sub_seq_2)


            # Select two occurances at random
            i1 = random.sample(occurances_1, 1)[0]
            i2 = random.sample(occurances_2, 1)[0]

            # Do not select overlapping sub-sequences
            if (i1 - 4 < i2) & (i2 < i1 + 4):
                while (i1 - 4 < i2) & (i2 < i1 + 4):
                    i1 = random.sample(occurances_1, 1)[0]
                    i2 = random.sample(occurances_2, 1)[0]

            # Orientate in correct order
            if i1 > i2:
                i1, i2 = i2, i1
                sub_seq_1, sub_seq_2 = sub_seq_2, sub_seq_1

            # Perform swap
            temp = sequence[0 : i1] + sub_seq_2 + sequence[i1 + 4 : i2] + sub_seq_1 + sequence[i2 + 4 : len(sequence)]
            
            # Check if Dinucleotide Counts are equal
            temp_dinuc = coding_notebook.dinuc_count(temp)

            # Print Error
            if RNA_sequence_dinuc != temp_dinuc:
                print('Dinucleotide frequencies do not match')
                break

            sequence = temp
        lst_of_sequences.append(sequence)
        
    return lst_of_sequences

##################################################################################################################################

def swap(sequence, swap_type='codon_dinuc', num_random_sequences=10, num_swaps=5):
    '''
    '''
    
    if swap_type == 'codon_dinuc':
        seqs = codon_dinuc_swap(sequence, num_random_sequences, num_swaps)
        
    elif swap_type == 'codon_pairwise':
        seqs = codon_pairwise_swap(sequence, num_random_sequences, num_swaps)
    
    elif swap_type == 'codon':
        seqs = codon_swap(sequence, num_random_sequences)
        
    elif swap_type == 'nt_pairwise':
        seqs = nt_pairwise_swap(sequence, num_random_sequences, num_swaps)
        
    elif swap_type == 'dinuc_pairwise':
        seqs = dinuc_pairwise_swap(sequence, num_random_sequences, num_swaps)
        
    else:
        print('Error')
    
    return seqs

##################################################################################################################################
