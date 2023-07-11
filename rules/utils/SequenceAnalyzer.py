import numpy as np
import pandas as pd
import pacmap as pm
import pathlib
import sklearn
from Bio import SeqIO
from Bio.Seq import Seq

class SequenceAnalyzer:
    """
    Organizes many sequences into simple strings,  integer- and one-hot-encoded arrays, and pd.DataFrame formats, and enables
    interconversion between these different formats, vectorized mutation analysis operations,
    and facilitates linkage between sequences and associated data such as fitness scores

    args:
        reference_fasta:    str, fasta file containing two or three sequences. The first is the sequence used for alignment,
                                the second is the sequence used for nucleotide-level analysis, the third is the nucleotide sequence
                                of an associated protein sequence (length%3 = 0)
        genotypesCSV:       str, .csv file name that contains a NT_substitutions, NT_insertions, and NT_deletions columns
        exclude_indels:     bool, whether to remove any sequences with indels from the dataset

    """

    def __init__(self, reference_fasta=None, genotypesCSV=None, exclude_indels=False, MSA_fasta=None):
        self.characters = {'NT':'ATGC-', 'AA':'ACDEFGHIKLMNPQRSTVWY*-'}

        # two different kinds of initialization, one for genotypes CSV input and one for MSA fasta input
        if MSA_fasta is None:
            if reference_fasta is None:
                raise ValueError("reference fasta file must be input if MSA fasta is not used for initialization")
            if genotypesCSV is None:
                raise ValueError("genotypes CSV file must be input if MSA fasta is not used for initialization")
            
            self.assign_genotypes_df(genotypesCSV, exclude_indels=exclude_indels)
            self.do_AA_analysis = False
            if 'AA_substitutions_nonsynonymous' in self.genotypes.columns:
                self.do_AA_analysis = True
            self.assign_reference(reference_fasta, include_AA=self.do_AA_analysis )
            self.integer_matrix = {'NT':self.seq_array_from_genotypes('NT')}
            if self.do_AA_analysis:
                self.integer_matrix['AA'] = self.seq_array_from_genotypes('AA')
        
        # TODO: add MSA initialization
        else:
            self.genotypes = None

        # convert integer array(s) to onehot array(s)
        self.onehot_matrix = {}
        for array_type in self.integer_matrix:
            self.onehot_matrix[array_type] = SequenceAnalyzer.integer_to_onehot_matrix(self.integer_matrix[array_type], len(self.characters[array_type]) )
        
        # make all indices the most recently selected indices
        key = list(self.integer_matrix.keys())[0]
        self.selected_idx = np.arange(self.integer_matrix[key].shape[0])
        self.selected_data = {}

    def assign_reference(self, reference_fasta, include_AA):
        """ stores sequences from reference fasta as SequenceEncoder objects

        args:
            reference_fasta:    str, fasta file containing two or three sequences. The first is the sequence used for alignment,
                                    the second is the sequence used for nucleotide-level analysis, the third is the nucleotide sequence
                                    of an associated protein sequence (length%3 = 0)
            include_AA:         bool, Whether to include the sequence corresponding to AA sequence. If True, will store the 2nd sequence in the fasta file as NT
                                    and the 3rd sequence in the fasta file as AA. Otherwise will only store the 2nd sequence in the fasta file
        """
        reference_sequences = [ str(s.seq).upper() for s in list(SeqIO.parse(reference_fasta, 'fasta')) ]
        self.ref_seq = {'NT':SequenceEncoder(reference_sequences[1], self.characters['NT'])}
        if include_AA:
            AA_ref = Seq.translate(reference_sequences[2])
            self.ref_seq['AA'] = SequenceEncoder(AA_ref, self.characters['AA'])

    def assign_genotypes_df(self, genotypesCSV, exclude_indels):
        """ store a genotypes csv file as a pd.DataFrame in self,
            first removing any sequences with indels if necessary
        params:
            genotypesCSV:       str, path to csv file containing genotypes
            exclude_indels:     bool, whether to exclude sequences with indels
        """
        genotypes = pd.read_csv(genotypesCSV)
        if exclude_indels:
            genotypes = genotypes[genotypes['NT_insertions'].isna()]
            genotypes = genotypes[genotypes['NT_deletions'].isna()]
        genotypes = genotypes.reset_index().drop('index', axis='columns')
        self.genotypes = genotypes

    def seq_array_from_genotypes(self, NTorAA):
        """
        args:
            NTorAA:             string, either 'NT' or 'AA'
        
        returns an array of integer-encoded sequences of
        shape (N,L) where:
                N = len(substitutionsList = len(deletionsList),
                L = len(reference_sequence)
        """
        if NTorAA=='NT':
            subs = 'NT_substitutions'
            dels = 'NT_deletions'
        if NTorAA=='AA':
            subs = 'AA_substitutions_nonsynonymous'
            dels = 'AA_deletions'
        genotypes_cols = self.genotypes[[subs, dels]]

        # make an array of integer encoded genotypes of shape (N,L)
        reference = self.ref_seq[NTorAA]
        integer_list = []
        for subs,dels in genotypes_cols.itertuples(index=False, name=None):
            integer = reference.genotype_modify(subs,dels)
            integer_list.append(integer)
        integer_matrix = np.array(integer_list)

        return integer_matrix

    @staticmethod
    def integer_to_onehot_matrix(integer_matrix, num_characters):
        """
        converts a 2D matrix of integer-encoded sequences of shape (N,L)
        into a 3D matrix of onehot-encoded sequences of shape (N,L,C) where:
            N = integer_matrix.shape[0]
            L = integer_matrix.shape[1]
            C = num_characters

        args:
            integer_matrix:     np.array, 2D matrix of integer-encoded sequences
            num_characters:     int, size to make the 3rd dimension of the output
                                    3D array. Must be at least 1 more than the maximum
                                    value in integer_matrix (because of 0-indexing)
        """
        onehot_matrix = np.zeros((integer_matrix.shape[0], integer_matrix.shape[1], num_characters), dtype=np.int8)
        onehot_matrix[np.arange(integer_matrix.shape[0])[:, np.newaxis], np.arange(integer_matrix.shape[1]), integer_matrix] = 1
        return onehot_matrix
    
    def select(self, idx=None):
        """
        returns integer matrices, onehot matrices, and, if it exists, the genotypes pd.DataFrame
        as a single dictionary, applies a selection to these if indices are provided, and
        saves the indices to self.selected_idx for future retrieval

        args:
            idx:    list-like, a list of integers referring to indexes within the integer matrix,
                        onehot matrix, and, if it exists, genotypes pd.DataFrame. If not provided,
                        no subselection is applied to the matrices or dataframe.
            new:    bool, whether to save the selection as a new selection and wipe associated data
                        or to keep the selection as the current selection and retain the associated data
        
        returns:
            dictionary, {'integer':selection from each integer 2D matrix in self.integer_matrix dict,
                        'onehot':selection from each onehot 3D matrix in self.onehot_matrix dict,
                        (optional)'df':subset of pd.DataFrame self.genotypes}
        """
        if idx is None:
            key = list(self.integer_matrix.keys())[0]
            idx = np.arange(self.integer_matrix[key].shape[0])
        else:
            idx = np.sort(np.array(idx))

        integer_selected = {key:self.integer_matrix[key][idx,:] for key in self.integer_matrix}
        onehot_selected = {key:self.onehot_matrix[key][idx,:,:] for key in self.onehot_matrix}
        out_dict = {'integer':integer_selected, 'onehot':onehot_selected}
        if self.genotypes is not None:
            out_dict['df'] = self.genotypes.iloc[idx]

        # if selection is new, assign as most recent selection and wipe associated data
        if not np.array_equal(idx, self.selected_idx):
            self.selected_idx = idx
            self.selected_data = {}
        return out_dict
    
    def get_selection(self):
        """
        returns the most recently retrieved selection indices of integer matrices, onehot matrices,
        and, if it exists, the genotypes pd.DataFrame.
        If additional data about the selection exists, adds this to new columns in the dataframe
        
        returns:
            dictionary, {'integer':selection from each integer 2D matrix in self.integer_matrix dict,
                        'onehot':selection from each onehot 3D matrix in self.onehot_matrix dict,
                        (optional)'df':subset of pd.DataFrame self.genotypes}
        """
        selected_dict = self.select(idx=self.selected_idx)
        for key in self.selected_data:
            selected_dict['df'][key] = self.selected_data[key]
        return selected_dict
    
    def downsample(self, size):
        """
        downsamples the dataset to a given size by randomly selecting a subset of the data,
            then returns the selected data via the select method. If the requested size is
            larger than the dataset, returns the entire dataset via the select method.

        args:
            size:   int, number of sequences to sample from the dataset
        """
        key = list(self.integer_matrix.keys())[0]
        if size > self.integer_matrix[key].shape[0]:
            idx = np.arange(self.integer_matrix[key].shape[0])
        else:
            idx = np.sort(np.random.choice(np.arange(self.integer_matrix[key].shape[0]), size=size, replace=False))
            
        return self.select(idx=idx)

    @staticmethod
    def pairwise_hamming_distance_matrix(seqArray):
        """
        Given an array of shape (N,L) containing N integer-encoded sequences of length L,
        computes the pairwise hamming distance of all pairs of sequences and outputs these as
        a matrix of shape (N,N) containing all pairwise hamming distances
        """
        # sklearn pairwise distances returns distances as a fraction of the maximum distance,
        #   so multiplying by the maximum distance, then rounding, and converting to int
        return np.rint(sklearn.metrics.pairwise_distances(seqArray,metric='hamming')*seqArray.shape[1]).astype(int)

    @staticmethod
    def bincount_2D(matrix):
        """
        matrix:     2d array

        given a matrix (2d array) of positive integer values of shape (M,N),
            returns bincounts for the matrix as a 2d array of shape
            (M, matrix.max()+1). Essentially a vectorized version of
            np.bincount applied to each element along the 0th dimension of a 2d array

        taken from https://stackoverflow.com/questions/40591754/vectorizing-numpy-bincount
        """
        M = matrix.shape[0]
        maxVal = matrix.max()+1
        matrix = matrix.transpose()

        arr = matrix + (maxVal*np.arange(M))
        return np.bincount(arr.ravel(), minlength=M*maxVal).reshape(M,-1)


    def HD_matrix_and_dist(self, NTorAA, downsample=False):
        """
        NTorAA:         string, 'NT' or 'AA', determines whether to calculate
                            NT hamming distance or AA hamming distance
        downsample:     False or int, if int, after removal of genotypes with insertions or deletions,
                            the number of sequences will be reduced down to this number
        """

        if downsample:
            self.downsample(downsample)

        HD_matrix = self.pairwise_hamming_distance_matrix(self.get_selection()['integer'][NTorAA])

        triangle = np.tril(HD_matrix)            # matrix is symmetric along diagonal, so zero out hamming distances from upper triangle
        HD_bincount = self.bincount_2D(triangle)
        HD_bincount[:,0] = np.subtract(HD_bincount[:,0], np.arange(HD_bincount.shape[0],0,-1))     # remove 0 counts resulting from zeroing out the upper triangle
        counts = self.get_selection()['df']['count'].to_numpy().reshape(len(self.get_selection()['df']),1)
        HD_bincount = np.multiply(HD_bincount, counts)                                            # multiply hamming distance bincounts for each sequence by the counts for each sequence
        same_pairs = np.multiply(counts, (counts-1)) / 2
        HD_bincount[:,0] = np.add(HD_bincount[:,0], same_pairs[:,0])                                                            # add appropriate amount of 0 counts for sequences with count>1 according to the formula pairs=count*(count-1)/2
        HD_dist = HD_bincount.sum(0)
        
        return HD_matrix, HD_dist
    
    def write_fasta(self, filename, NTorAA, idx=None, unique_only=True):
        """
        writes a fasta file of the nucleotide or protein sequences in the dataset

        args:
            filename:           string, path to the output file
            NTorAA:             string, either 'NT' or 'AA'
            idx:                list-like, a list of integers referring to indices of sequences in the 
                                    integer matrix to be written. If not provided,
                                    all sequences are written
            unique_only:        bool, whether to write each unique sequence only once or to write each sequence
                                    the number of times it appears in the dataset
        """
        if idx is None:
            idx = np.arange(self.integer_matrix.shape[0])

        if self.genotypes is not None: # add genotype ID to fasta headers
            ID_list = self.genotypes['genotype_ID'].values

        pathlib.Path(filename).parent.absolute().mkdir(parents=True, exist_ok=True)
        
        with open(filename, 'w') as f:
            for i in idx:
                ID_line = f'>index:{i}'
                if self.genotypes is not None:
                    ID_line += f'_genotype-ID:{ID_list[i]}'
                integer_encoded_seq = self.integer_matrix[NTorAA][i,:]
                sequence_line = ''.join( self.ref_seq[NTorAA].decode_integer(integer_encoded_seq) )

                count = 1
                if ( not unique_only ) and ( self.genotypes is not None ):
                    count = self.genotypes['count'][i]
                
                for _ in range(count):
                    f.write(f'{ID_line}\n{sequence_line}\n')


    def get_count(self, idx=None):
        """
        returns the number of sequences in the dataset

        args:
            idx:    list-like, a list of integers referring to indexes within the integer matrix,
                        onehot matrix, and, if it exists, genotypes pd.DataFrame. If not provided,
                        no subselection is applied to the matrices or dataframe.
        
        returns:
            int, number of sequences in the dataset
        """
        if idx is None:
            idx = np.arange(self.integer_matrix.shape[0])
        subset = self.select(idx)

        if 'df' in subset:  # if genotypes DF is present, use the count column
            return subset['df']['count'].sum()
        else:               # otherwise, just take the length of an integer matrix (all are the same length)
            key = list(subset['integer'].keys())[0]
            return subset['integer'][key].shape[0]
    
    def get_mutations_of_interest(self, NTorAA, muts, max_groups, idx=None):
        """
        generate a list of the NT or AA Mutations Of Interest as a list of strings containing comma separated mutations
            
        parameters:
            NTorAA (str):            "NT" or "AA"
            muts (str):             comma separated list of NT or AA mutations of interest (must match NTorAA)
                                        can be a full mutation, or just the position, or just a
                                        mutated AA/NT. Underscore can be used in place of a position number
                                        to include any position. Input should be 1-index examples:
                                            A100G, specific A100G mutation
                                            150 or 150_, mutation to position 150
                                            G_, mutation to any G
                                            _C, any mutation that produces a C
                                            A_C, any A->C mutation
            max_groups (int):        maximum number of unique groups to include, decided based on frequency
                                            the least common groups will be combined under the name 'other',
                                            and sequences without any of the listed mutations will be given the name 'none'
            idx (list(ints)):         indices of selected sequences to use. All sequences are used if not provided

        returns:
            list of strings, each string is a comma separated list of mutations of interest, the most frequently
                appearing combinations (n=max_groups) are specified, all others are combined under the name 'other'.
                Length of output list is equal to len(idx) if idx is provided, otherwise onehot_matrix.shape[0]
        """
        column_name = f'{NTorAA}_muts_of_interest'

        # if no mutations of interest are specified, return a list of 'none' strings
        if muts == '':
            out_length = len(idx) if idx is not None else self.integer_matrix[NTorAA].shape[0]
            out = ['none']*out_length
            self.selected_data[column_name] = out
            return out

        onehot_matrix = self.select(idx)['onehot'][NTorAA]
        N,L,C = onehot_matrix.shape
        all_L_idx = np.arange(L)
        all_C_idx = np.arange(C)
        
        muts_list = muts.replace(' ','').split(',')
        number_strings = [str(x) for x in range(0,10)]
        allowed_characters = ['_', *self.characters[NTorAA], *number_strings]
        
        # build a 2D array of shape (L,C) of onehot encoded mutations of interest
        #    to be broadcast and multiplied by the 3D onehot array of shape (N,L,C)
        #    to find genotypes that have one or more of the mutations in muts
        MOIs_2D_arr = np.zeros((L,C), dtype=np.int8)
        
        for mut in muts_list: # parse input mutations

            if any([s not in allowed_characters for s in mut]):
                raise ValueError(f'Invalid characters in provided mutation {mut}. Please check your input.')
            
            wt_input = mut[0] if mut[0] in self.characters[NTorAA] else ''
            mut_input = mut[-1] if mut[-1] in self.characters[NTorAA] else ''
            position = mut.strip(wt_input+mut_input)

            C_idx = all_C_idx
            if mut_input:
                # set mut idx to all character positions that are mut_input
                C_idx = self.characters[NTorAA].find(mut_input[-1])

            if position == '_':
                if wt_input == mut_input == '':
                    raise ValueError(f'Underscore must be used with a wt or mut input. Please check your input.')
                
                L_idx = all_L_idx
                if wt_input:
                    # set wt idx to all ref seq positions that are wt_input
                    L_idx = np.where(self.ref_seq[NTorAA].array == wt_input)[0]

            else:
                L_idx = int(position)-1
                wt = self.ref_seq[NTorAA].string[L_idx]
                if wt != wt_input:
                    print(f'[WARNING] User provided wt nucleotide/residue for {NTorAA} mutation {mut} does not match the wt nucleotide/residue {wt}. Please double check indexing.')
                
            MOIs_2D_arr[L_idx,C_idx] = 1
            
        # zero out all non-mutations by subtracting out the onehot encoded wild type sequence
        #    then converting any -1s to 0
        MOIs_2D_arr = np.subtract( MOIs_2D_arr, self.ref_seq[NTorAA].onehot )
        MOIs_2D_arr[MOIs_2D_arr<0] = 0

        # perform matrix multiplication between the 2D onehot array representing all mutations 
        #    of interest and the 3D onehot array representing all sequences. Resulting 3D array
        #    will have ones at positions where a sequence contained a mutation of interest, 0s elsewhere
        MOIs_2D_arr = (MOIs_2D_arr==1) # mask non-1 values as False so these positions are not multiplied. positions in resulting array are 0. ~10-fold speedup in einsum
        all_MOIs_3D_arr = np.einsum('NLC,LC->NLC', onehot_matrix, MOIs_2D_arr)
        
        # need to find all unique combinations of mutations of interest but np.unique can be sped up ~10-fold by
        #   getting rid of values that we don't care about, so we make a new array
        #   in which only mutation positions are present and call np.unique on that

        # list to store subset of all_MOIs_3D_arr
        all_MOIs_subset = []

        for (i,j), value in np.ndenumerate(MOIs_2D_arr):
            if value == 1:
                all_MOIs_subset.append(all_MOIs_3D_arr[:,i,j])
        all_MOIs_subset = np.stack(all_MOIs_subset, axis=1)
        _, idx_first, idx_original, counts = np.unique(all_MOIs_subset, axis=0, return_index=True, return_inverse=True, return_counts=True)

        # retrieve the original 2D arrays representing all unique combinations of mutations of interest that are observed,
        #   sort by count, iterate through the most frequently observed, and label with human readable strings 
        unique_MOI_combo_arrs = all_MOIs_3D_arr[idx_first,:,:]
        sorted_idx = np.argsort(-counts)
        MOI_combo_strings_sorted = []
        none_found = False # flag for doing a search for the MOI_combo without any mutations of interest

        for combo_idx, unique_idx in enumerate(sorted_idx):
            
            # if no more groups need to be made, fill remaining values with 'other'
            if (combo_idx+1) > max_groups: # convert to 1-indexed
                remaining_genotypes = len(sorted_idx) - len(MOI_combo_strings_sorted)
                MOI_combo_strings_sorted.extend(['other']*(remaining_genotypes))
                break

            combo_arr = unique_MOI_combo_arrs[unique_idx,:,:]
            mut_indices = np.argwhere(combo_arr == 1)
            if len(mut_indices) > 0:
                muts = []
                for L_idx, C_idx in mut_indices:
                    wt = self.ref_seq[NTorAA].string[L_idx]
                    mut = self.characters[NTorAA][C_idx]
                    muts.append(wt + str(L_idx+1) + mut)
                mut_string = ', '.join(muts)
            else:
                max_groups += 1
                none_found = True
                mut_string = 'none'

            MOI_combo_strings_sorted.append(mut_string)

        MOI_combo_strings_sorted = np.array(MOI_combo_strings_sorted)
        MOI_combo_strings_unsorted = MOI_combo_strings_sorted[np.argsort(sorted_idx)] # convert to original order of unique arrays with reverse indexing

        # search for unique array without any MOIs if it hasn't already been found
        if not none_found:
            none_arr = np.zeros( (len(self.ref_seq[NTorAA].string), len(self.characters[NTorAA])) )
            none_arr_matches = np.any(np.all( unique_MOI_combo_arrs == none_arr, axis=(1,2)), axis=0)
            none_idx = np.argmax(none_arr_matches) if none_arr_matches.any() else None
            if none_idx!=None:
                MOI_combo_strings_unsorted[none_idx] = 'none'

        MOI_combo_strings_original = MOI_combo_strings_unsorted[idx_original] # convert to original order of all arrays/genotypes with reverse indexing

        self.selected_data[column_name] = MOI_combo_strings_original
        return MOI_combo_strings_original
    
    def aggregate_identities(self, NTorAA, idx=None):
        """
        Generate a numpy array of counts of all NT or AA identities in a selection by summing
            the onehot sequence array along axis 0
        
        Parameters:
            NTorAA (str):           'NT' or 'AA' to indicate nucleotide or amino acid mutation level aggregation
            idx (np.array):         1d array of indices of a selection. If None, all sequences are used
        
        Returns:
            np.array: A 2D array of shape L,C that tabulates counts for all possible mutations
        """

        selected = self.select(idx)
        onehot_matrix = selected['onehot'][NTorAA]
        if 'df' in selected:
            counts = selected['df']['count'].values.reshape(-1,1,1)
            onehot_matrix = onehot_matrix * counts
        aggregated_identities = np.sum(onehot_matrix, axis=0, dtype=np.int32)

        return aggregated_identities
    
    def aggregate_mutations(self, NTorAA, idx=None):
        """
        Generate a tidy format pd.DataFrame of counts of all NT or AA mutations 
        
        Parameters:
            NTorAA (str):           'NT' or 'AA' to indicate nucleotide or amino acid mutation level aggregation
            idx (np.array):         1d array of indices of a selection. If None, all sequences are used
        
        Returns:
            pd.DataFrame: A tidy-formatted dataframe that tabulates counts for all possible mutations
        """

        aggregated_identities = self.aggregate_identities(NTorAA, idx=idx)
        total_seqs = self.get_count(idx=idx)
        
        rows = []
        chars = self.characters[NTorAA]

        # loop through all possible mutations and generate a dataframe of all, including those with count of 0
        for index, count in np.ndenumerate(aggregated_identities):
            posi, char_idx = index
            wt = self.ref_seq[NTorAA].string[posi]
            posi += 1
            mut = chars[char_idx]

            # don't add a row for the wildtype sequence identities
            if wt != mut:
                rows.append([wt, posi, mut, count])
                
        df = pd.DataFrame(rows, columns=['wt', 'position', 'mutation', 'total_count'])
        df['proportion_of_seqs'] = df['total_count']/total_seqs

        return df
    
    def get_consensus(self, NTorAA, idx=None, write_to='', append=False, name=''):
        """
        Get the consensus sequence for a selection of sequences
        
        Parameters:
            NTorAA (str):           'NT' or 'AA' to indicate nucleotide or amino acid mutation level aggregation
            idx (np.array):         1d array of indices of a selection. If None, all sequences are used
            write_to (str):         If not empty, write the consensus sequence to the given .fasta file
            append (bool):          If True, append to the given .fasta file instead of overwriting it.
                                        Does nothing if string is empty
            name (str):             Name to prepend to the consensus sequence in the .fasta file
        
        Returns:
            str: The consensus sequence
        """

        aggregated_identities = self.aggregate_identities(NTorAA, idx=idx)
        integer_encoded_consensus = np.argmax(aggregated_identities, axis=1)
        consensus = ''.join(SequenceEncoder.convert_array(integer_encoded_consensus, self.ref_seq[NTorAA].decoder_dict))

        if write_to:
            pathlib.Path(write_to).parent.absolute().mkdir(parents=True, exist_ok=True)
            total_seqs = self.get_count(idx=idx)
            write_mode = 'a+' if append else 'w'
            with open(write_to, write_mode, encoding="utf-8") as f:
                if name:
                    name = name + '_'
                f.write(f'>{name}{NTorAA}_consensus_sequence_from_{total_seqs}_sequences\n{consensus}\n')

        return consensus
    
    def dimension_reduction(self, NTorAA, idx=None, onehot=True, encoding=None, dimensions=2):
        """ 
        Perform dimension reduction (DR) on the sequences of a particular level (NT or AA), and return the result as a DataFrame
        
        Parameters:
            NTorAA (str):           'NT' or 'AA' to indicate DR applied to nucleotide or amino acid sequences
            onehot (bool):          whether to use onehot-encoding of the sequence matrix or not. Must be False for a non-integer encoding.
                                        Should be set to true unless the encoding being used has meaningful ordering
                                        (e.g. defining some amino acid property on a scale)
            encoding (dict):        encoding dictionary to use. Keys correspond to the original sequence, values correspond to the final value.
                                        If None, the standard integer encoding is used
            dimensions (int):       number of dimensions to reduce to
        
        Returns:
            pd.DataFrame:           DataFrame of shape (n_sequences, dimensions) containing the DR output
        """
        if idx is None:
            idx = np.arange(self.integer_matrix[NTorAA].shape[0])
        
        # behavior depends on both onehot and encoding. if both are default, then just use the stored onehot matrix instead of rebuilding it
        if encoding is None:
            if onehot:
                matrix = self.select(idx)['onehot'][NTorAA]
            else: 
                matrix = self.select(idx)['integer'][NTorAA]

        else:
            matrix = self.select(idx)['integer'][NTorAA]
            # convert encoding to an encoding that can be used directly on the stored integer matrix
            letter_to_int = self.ref_seq[NTorAA].encoder_dict
            encoding = {letter_to_int[k]:v for k,v in encoding.items()}
            matrix = SequenceEncoder.convert_array(matrix, encoding)

            if onehot: # convert to onehot encoding
                # this is not a good assumption, but is just a maximum. if there are fewer values, then the blank values will be removed prior to DR
                num_characters = len(encoding)

                matrix = SequenceAnalyzer.integer_to_onehot_matrix(matrix, num_characters)            

        if onehot:
            # flatten onehot matrix, remove axis 1 positions that are all 0s. accomplishes a similar task to the 'else' block but is faster and only works for onehot encoding
            N,L,C = matrix.shape
            matrix = matrix.reshape((N,L*C))
            cols_with_nonzero = np.where(np.any(matrix, axis=0))[0]
            matrix = matrix[:, cols_with_nonzero]
        else:
            # remove values that are the same across axis 1 (i.e. all sequences have the same value for that position)
            same_values = np.all(matrix == matrix[0, :], axis=0)
            different_values = np.logical_not(same_values)
            matrix = matrix[:, different_values]
            
        # initializing the pacmap instance
        embedding = pm.PaCMAP(n_components=dimensions, MN_ratio=0.5, FP_ratio=2.0)

        # fit the data to provided number of dimensions then make into a dataframe with # of columns = dimensions
        reduced = embedding.fit_transform(matrix, init="pca")
        columns = [f'{NTorAA}_PaCMAP{i}' for i in range(1,dimensions+1)]
        reduced_DF = pd.DataFrame(reduced, columns=columns, index=self.genotypes.index)
        return reduced_DF
    
    def assign_dimension_reduction(self, NTorAA, onehot=True, encoding=None, dimensions=2):
        """
        Perform dimension reduction (DR) on the sequences of a particular level (NT or AA), and assign as new columns
            to self.genotypes
        
        Parameters:
            NTorAA (str):           'NT' or 'AA' to indicate DR applied to nucleotide or amino acid sequences
            onehot (bool):          whether to use onehot-encoding of the sequence matrix or not. Must be False for a non-integer encoding.
                                        Should be set to true unless the encoding being used has meaningful ordering
                                        (e.g. defining some amino acid property on a scale)
            encoding (dict):        encoding dictionary to use. Keys correspond to the original sequence, values correspond to the final value.
                                        If None, the standard integer encoding is used
            dimensions (int):       number of dimensions to reduce to
        
        """
        reduced_DF = self.dimension_reduction(NTorAA, onehot=onehot, encoding=encoding, dimensions=dimensions)
        self.genotypes = pd.merge(self.genotypes, reduced_DF, how='left', left_index=True, right_index=True)


class SequenceEncoder:
    """ encodes a sequence into string, integer, and onehot arrays
    args:
        sequence (str): sequence to be encoded
        characters (str): characters to be used for encoding.
            Positions of characters within string is used for encoding
            Allows for consistency in encoding among different input sequences
    
    """
    def __init__(self, sequence, characters):

        self.string = sequence.upper()
        
        #get sequence into an array
        self.array = np.array(list(self.string))

        # create an encoder dict for consistent encoding
        self.encoder_dict = {char:index for index,char in enumerate(characters.upper())}
        self.decoder_dict = {index:char for index,char in enumerate(characters.upper())}
        
        # create string, integer- and onehot-encoded arrays for the sequence and add to self
        self.assign_integer()
        self.assign_onehot()

    @staticmethod
    def convert_array(arr, converter_dict):
        """ converts an array of strings or integers into an array of the same shape with converted values using the provided dictionary.
                converter_dict keys correspond to values of the input array and values correspond to values of the output array.
        args:
            arr_1D (np.array): 1D array of values to be converted
            dict (dict): dictionary to use for conversion

        returns:
            np.array: 1D array of converted values
        
        example:
            arr_1D = np.array(['A', 'T', 'G', 'C'])
            dict = {'A':0, 'T':1, 'C':2, 'G':3}
            convert_array(arr_1D, dict)
            >>> array([0, 1, 3, 2])
        
        """

        output_dtype = np.array(list(converter_dict.values())).dtype
        
        # Create an empty output array with the same shape as the input array and the inferred data type
        output_array = np.empty_like(arr, dtype=output_dtype)
        
        # Iterate through the converter dictionary and apply the conversion
        for key, value in converter_dict.items():
            output_array[np.where(arr == key)] = value

        return output_array

    def encode_integer(self, sequence_array):
        """ encodes an array of strings into an array of integers using the encoder_dict """
        return self.convert_array(sequence_array, self.encoder_dict)
    
    def assign_integer(self):
        """ assign the integer-encoded sequence to the sequence object """
        self.integer = self.encode_integer(self.array)
    
    def decode_integer(self, integer_array):
        """ decodes an array of integers into an array of strings using the decoder_dict """
        return self.convert_array(integer_array, self.decoder_dict)

    def integer_to_onehot(self, integer):
        """ converts a 1D integer encoded sequence to a 2D one-hot encoded sequence """
        onehot = np.zeros((len(integer), len(self.encoder_dict)), dtype=np.int8)
        onehot[np.arange(len(integer)), integer] = 1
        return onehot
    
    def assign_onehot(self):
        """ assign the one-hot-encoded sequence to the sequence object """
        self.onehot = self.integer_to_onehot(self.integer)
    
    def genotype_modify(self, substitutions, deletions):
        """
        modify an integer-encoded sequence according to provided substitutions and deletions
                args:
                    substitutions    string, comma separated list of 1-index mutations, ex: "A234C, T301G"
                                        first and last characters must be present in `characters`
                    deletions        string, comma separated list of 1-index nt deletions, ex: "200del3, 206del1"
                returns:
                    integer          modified integer-encoded sequence
        """
        integer = self.integer.copy()

        if type(substitutions) == str:
            for mut in substitutions.replace(' ','').split(','):
                posi = int(mut[1:-1])-1 # convert to 0-index
                mutNT = mut[-1]
                integer[posi] = self.encoder_dict[mutNT]
        if type(deletions) == str:
            for deletion in deletions.replace(' ','').split(','):
                posi, number = [int(x) for x in deletion.split('del')]
                posi = posi-1 # convert to 0-index
                integer[posi:posi+number] = self.encoder_dict['-']
        
        return integer
    