import pandas as pd
import hhelper
import htf
import numpy as np

class Data:

    def __init__(self, dir_data, dir_doid,max_depth):
        self.data_mit = pd.read_csv(dir_data)
        self.doid_labels = pd.read_pickle(dir_doid)
        self.max_depth = max_depth

    def filter_data(self, tpm_cutoff=1):
        # idx for exons with TPM >1
        # TODO: investigate with 1 is a good cutoff
        rm_me = list(self.data_mit.max(axis=1) > tpm_cutoff)
        # set exon names as index
        self.data_mit.set_index('exons', inplace=True, drop=False)
        # select exons with TPM >1
        self.data_mit = self.data_mit.loc[rm_me]
        # select exons with valid labeling
        self.data_mit = self.data_mit[~pd.isna(self.data_mit['labels']) == True]
        # send to preprocess (drops zeros, ambigious exons (exons annotated with multiple diseases)
        # and does max-min scaling)
        self.data_mit = htf.pre_process(self.data_mit)
        # convert label string to list
        self.data_mit['labels'] = self.data_mit['labels'].apply(hhelper.str_to_list)
        # list of lenghts (i.e. depth) of label list per exon
        labeling_depth = list(map(len, self.data_mit['labels'].values))
        # idx for exons with depth >= XX
        min_length = list(map(lambda x: x >= self.max_depth + 1, labeling_depth))
        # select exons with min depth
        self.data_mit = self.data_mit[min_length]

    def build_tree(self):
        df_tree = hhelper.build_tree(self.data_mit, self.doid_labels)

        # this tree has to be pruned now to the desired subtree
        # 'disease' is the only entry on lvl 0
        # drop all lvl 1 nodes that are not XX and 'disease'
        use_this_tree = df_tree.drop(df_tree[(df_tree['label'] != 'disease of anatomical entity') &
                                             (df_tree['lvl'] <= 1)].index)

        # now that we removed some of the top nodes, let's remove all their children
        # starting one lvl below new root node to not remove root node
        for lvl in range(2, self.max_depth + 1):
            # remove parent less nodes
            use_this_tree.drop(use_this_tree[~use_this_tree['parent'].isin(use_this_tree['id'].values) &
                                             (use_this_tree['lvl'] == lvl)].index, inplace=True)

        # drop all nodes below lvl 1 that have less than XX members
        use_this_tree.drop(use_this_tree[(use_this_tree['count_all'] < 50) &
                                         (use_this_tree['lvl'] > 1)].index, inplace=True)
        # drop all nodes below max depth
        use_this_tree.drop(use_this_tree[use_this_tree['lvl'] > self.max_depth].index, inplace=True)

        # remove leaf nodes w/o siblings
        number_of_sib = []
        for label in use_this_tree['label'].values:
            # get int value of parent for given label and the number of siblings for label
            parent = use_this_tree.loc[use_this_tree['label'] == label, 'parent'].values[0]
            number_of_sib.append(use_this_tree[use_this_tree['parent'] == parent].shape[0])

        # number_of_sib if 1 it has 0 siblings, the 1 come from counting the node itself
        use_this_tree['number_of_sib'] = number_of_sib
        use_this_tree.drop(use_this_tree[(use_this_tree['number_of_sib'] < 2) &
                                         (use_this_tree['lvl'] == self.max_depth)].index, inplace=True)

        # remove children less nodes from bottom up
        for lvl in range(self.max_depth - 1, 1, -1):
            print(lvl)
            use_this_tree.drop(use_this_tree[~use_this_tree['id'].isin(use_this_tree['parent'].values) &
                                             (use_this_tree['lvl'] == lvl)].index, inplace=True)

        # reindex df such that idx corresponds to axis=1 cooridnate in label matrix
        # UWAGA:the root node (idx=0) is not a label and labels are therefor 1 indexed
        use_this_tree = use_this_tree.sort_values('lvl')
        use_this_tree['idx'] = range(use_this_tree.shape[0])
        use_this_tree.set_index('idx', inplace=True)

        # now make sure we only have data that is labeled for this tree
        # get all the leaf nodes
        leaf_nodes = use_this_tree.loc[use_this_tree['lvl'] == self.max_depth, 'label'].values
        # get the leaf node label for all data
        list_leaf_data = [x[4] for x in self.data_mit['labels'].values]
        # does the exon labeling contain one of the leaf nodes?
        bool_keep_data = [True if x in leaf_nodes else False for x in list_leaf_data]
        self.data_mit = self.data_mit[bool_keep_data]

        return use_this_tree

    def label_matrix(self, use_this_tree):
        label_matrix = np.zeros((self.data_mit.shape[0], use_this_tree.shape[0] - 1))

        # TODO: find a way to avoid at least one loop

        # loop through each exon and assign labels to label_matrix
        for i, list_labels in enumerate(self.data_mit['labels']):
            # loop through the list starting on the first lvl after root node
            for label in list_labels[2:5]:
                # get label idx
                idx = use_this_tree[use_this_tree['label'] == label].index[0]
                # add to idx-1 (shift due to root node not beign a class)
                label_matrix[i, idx - 1] = 1

        # drop meta data and hstack to datamatrix
        # for_np = self.data_mit.drop(['labels', 'COID', 'phenString', 'DOID', 'exons'], axis=1)
        # np_data_mit = np.hstack([for_np, label_matrix])
        return label_matrix