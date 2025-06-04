import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import copy
import sys
import random

def remove_dups(df, graph_type='undirected'):
    '''
    Inputs: df (dataframe) containing at least two columns named 'source' and 'target'. An optional column can be 'edge_type'.
            Depending on the presence of 'edge_type' a pair, e.g., (a,b) or a triplet, e.g., (a,b,C) can represent an edge.
            Here, a = source node, b = target node, C = edge type.
    graph_type (str): 'directed' or 'undirected'

    Function: Remove redundant/repeated edges.
    - When 'edge_type' is present then keep only one instance of (a, b, C). Note that, (a, b, C) and (a, b, D) can
     remain in the final dataset where D is another edge_type.
    - When 'edge_type' is absent then keep only one instance of (a,b).
    - If it is an undirected graph, then keep either (a, b, C) or (b, a, C) (and (a, b) or (b, a)).
    '''

    if 'edge_type' in df.columns:
        df.drop_duplicates(subset=['source', 'target', 'edge_type'], keep='first', inplace=True)
    else:
        df.drop_duplicates(subset=['source', 'target'], keep='first', inplace=True)

    if graph_type == 'undirected':
        df['sorted_col1'] = np.minimum(df['source'], df['target'])
        df['sorted_col2'] = np.maximum(df['source'], df['target'])

        # Drop duplicates based on the sorted columns
        df = df.drop_duplicates(subset=['sorted_col1', 'sorted_col2']).drop(
            columns=['sorted_col1', 'sorted_col2'])
    return df


def split_list(d_list, n_folds, seed):
    '''
    :param d_list: list of nodes, edges, or edge types.
    :param n_folds: In how many folds to split the values in the list.
    :return: a dict where key=fold number, value=list of values in that fold.
    '''
    #split in n_folds
    random.Random(seed).shuffle(d_list)

    split_size= int(len(d_list) / n_folds)
    folds= {i: d_list[split_size * i: split_size * (i + 1)] for i in range(n_folds-1)}
    folds[n_folds-1] = d_list[split_size*(n_folds-1):]

    return folds

def verify_split(df, train_idx, test_idx, split_type):
    '''
    :return: Given df, train idx and test_idx, verify that the train and test
    data have been created properly according to the split type.'''
    train_df = df.iloc[train_idx]
    test_df = df.iloc[test_idx]

    if split_type=='random':
        if 'edge_type' in df.columns:
            test_triplets = set(zip(test_df['source'],test_df['target'],test_df['edge_type']))
            train_triplets = set(zip(train_df['source'],train_df['target'],train_df['edge_type']))
            n_common = len(test_triplets.intersection(train_triplets))

        else:
            test_edges = set(zip(test_df['source'], test_df['target']))
            train_edges = set(zip(train_df['source'], train_df['target']))
            n_common = len(test_edges.intersection(train_edges))
        assert n_common == 0, print(f'error in {split_type} split')


    if split_type=='edge':
        test_edges = set(zip(test_df['source'],test_df['target']))
        train_edges = set(zip(train_df['source'],train_df['target']))
        n_common = len(test_edges.intersection(train_edges))
        assert n_common == 0, print(f'error in {split_type} based split')


    if split_type=='node':
        test_nodes = set(test_df['source']).union(set(test_df['target']))
        train_nodes = set(train_df['source']).union(set(train_df['target']))
        n_common = len(test_nodes.intersection(train_nodes))
        assert n_common == 0, print(f'error in {split_type} based split')


    if split_type=='edge_type':
        test_edge_type = set(test_df['edge_type'])
        train_edge_type = set(train_df['edge_type'])
        n_common = len(test_edge_type.intersection(train_edge_type))
        assert n_common == 0, print(f'error in {split_type} based split')


def split_random_train_test(df, test_frac, seed=None):
    df['ID'] = list(range(len(df)))
    if 'edge_type' in df.columns:
        edges = list(set(zip(df['source'], df['target'], df['edge_type'])))
    else:
        edges = list(zip(df['source'], df['target']))
    # as I used set() to get unique edges, we will lose the initial order of edges. Hence, we need to sort the edges. Otherwise the same edges can appear in different
    # order in different runs, thus using random_state=int(e.g., 42) will not ensure the reproducible split each time on the same dataset.
    edges.sort()
    #split edges into train and test set
    train_edges, test_edges = train_test_split(edges, test_size=test_frac, random_state=seed)
    if 'edge_type' in df.columns:
        test_idx = list(df[df[['source','target','edge_type']].apply(tuple, axis=1).isin(test_edges)]['ID'])
    else:
        test_idx = list(df[df[['source','target']].apply(tuple, axis=1).isin(test_edges)]['ID'])

    # train_idx = list(df['ID'].drop(test_idx))
    train_idx = list(df['ID'][~df['ID'].isin(test_idx)])
    df.drop(columns='ID', inplace=True)

    return train_idx, test_idx

def split_edge_train_test(df, test_frac, seed=None):
    ''' Input: df = a dataframe with atleast two columns ['source','target'].
        Function: edge appearing in train (irrespective of the edge type) will not appear in test. '''

    df['ID'] = list(range(len(df)))

    edges = list(set(zip(df['source'], df['target'])))  #make sure that edges are unique.
    # as I used set() to get unique edges, we will lose the initial order of edges. Hence, we need to sort the edges. Otherwise the same edges can appear in different
    # order in different runs, thus using random_state=int(e.g., 42) will not ensure the reproducible split each time on the same dataset.
    edges.sort()
    #split edges into train and test set
    train_edges, test_edges = train_test_split(edges, test_size=test_frac, random_state=seed)

    test_idx = list(df[df[['source','target']].apply(tuple, axis=1).isin(test_edges)]['ID'])
    # train_idx = list(df['ID'].drop(test_idx))
    train_idx = list(df['ID'][~df['ID'].isin(test_idx)])
    df.drop(columns='ID', inplace=True)

    return train_idx, test_idx



def split_node_train_test(df, test_frac, seed=None):

    '''
    Function: Nodes appearing in training data will not appear in any edge in the test.
    '''
    df['ID'] = list(range(len(df)))

    nodes = list(set(df['source']).union(set(df['target'])))
    # as I used set() to get unique nodes, we will lose the initial order of nodes. Hence, we need to sort the nodes. Otherwise the same nodes can appear in different
    # order in different runs, thus using random_state=int(e.g., 42) will not ensure the reproducible split each time on the same dataset.
    nodes.sort()

    count=0
    while True:
        count+=1

        # split nodes into train and test set
        train_nodes, test_nodes = train_test_split(nodes, test_size=test_frac, random_state=seed)

        #both source and target has to be in test_nodes for a triplet or edge to be considered in test dataset.
        test_idx = list(df[(df['source'].isin(test_nodes)) & (df['target'].isin(test_nodes))]['ID'])

        train_idx = list(df[(df['source'].isin(train_nodes)) & (df['target'].isin(train_nodes))]['ID'])

        if (len(test_idx)>0) and (len(train_idx)>0): #make sure that we don't get empty test or train dataset
            break

        if isinstance(seed, (int)):
            print('Failed to split nodes such that train and test dataset are not empty. Try using random_state=None or np.random.RandomState or a different integer value than the current one ')
            sys.exit(1)

        if count>1000:
            print('Failed to split nodes such that train and test dataset are not empty. Tried 1000 times ')
            sys.exit(1)

    df.drop(columns='ID', inplace=True)
    return train_idx, test_idx




def split_edge_type_train_test(df, test_frac, seed=None):

    '''
    Edge type appearing in training data will not appear in the test.
    '''
    df['ID'] = list(range(len(df)))

    edge_types = list(df['edge_type'].unique())
    edge_types.sort()


    # split nodes into train and test set
    train_edge_types, test_edge_types = train_test_split(edge_types, test_size=test_frac, random_state=seed)

    # both source and target has to be in test_nodes for a triplet or edge to be considered in test dataset.
    test_idx = list(df[(df['edge_type'].isin(test_edge_types))]['ID'])
    # train_idx = list(df['ID'].drop(test_idx))
    train_idx = list(df['ID'][~df['ID'].isin(test_idx)])

    df.drop(columns='ID', inplace=True)
    return train_idx, test_idx



def split_random_cv(df, n_folds, seed=None):
    df['ID'] = list(range(len(df)))
    if 'edge_type' in df.columns:
        edges = list(set(zip(df['source'], df['target'], df['edge_type'])))  # list of tuples
    else:
        edges = list(set(zip(df['source'], df['target'])))  # list of tuples
    edges.sort()

    # prepare train and val split
    edge_folds = split_list(edges, n_folds, seed=seed)

    if 'edge_type' in df.columns:
        df_split = {i: df[df[['source', 'target','edge_type']].apply(tuple, axis=1).isin(edge_folds[i])] for i in range(n_folds)}
    else:
        df_split = {i: df[df[['source', 'target']].apply(tuple, axis=1).isin(edge_folds[i])] for i in range(n_folds)}

    train_idx = {}
    val_idx = {}
    for i in range(n_folds):
        val_idx[i] = list(df_split[i]['ID'])
        # train_idx[i] = list(df['ID'].drop(val_idx[i]))
        train_idx[i] = list(df['ID'][~df['ID'].isin(val_idx[i])])

    df.drop(columns='ID', inplace=True)

    return train_idx, val_idx


def split_edge_cv(df, n_folds, seed=None):
    '''
        Input: df = a dataframe with atleast two columns ['source','target'].
        Function: edge appearing in train (irrespective of the edge type) will not appear in test.
        Return: 
    '''
    df['ID'] = list(range(len(df)))

    edges = list(set(zip(df['source'], df['target'])))  # list of tuples
    edges.sort()

    #prepare train and val split
    edge_folds = split_list(edges, n_folds, seed=seed)

    df_split = {i: df[df[['source','target']].apply(tuple, axis=1).isin(edge_folds[i])] for i in range(n_folds)}

    train_idx = {}
    val_idx = {}
    for i in range(n_folds):
        val_idx[i] = list(df_split[i]['ID'])
        # train_idx[i] = list(df['ID'].drop(val_idx[i]))
        train_idx[i] = list(df['ID'][~df['ID'].isin(val_idx[i])])

    df.drop(columns='ID', inplace=True)

    return train_idx, val_idx



def split_node_cv(df, n_folds, seed=None):

    '''
    Nodes appearing in training data will not appear in any edge in the test.
    '''

    df['ID'] = list(range(len(df)))

    nodes = list(set(df['source']).union(set(df['target']))) # list of tuples
    nodes.sort()

    # now create 5 folds from the remaining drug and use each fold as val and the rest as train
    node_folds = split_list(nodes, n_folds, seed)

    df_split = {i: df[(df['source'].isin(node_folds[i])) & (df['target'].isin(node_folds[i]))]  for i in range(n_folds)}

    train_idx = {}
    val_idx = {}
    for i in range(n_folds):
        val_idx[i] = list(df_split[i]['ID'])
        train_df = df[~(df['source'].isin(node_folds[i]) | df['target'].isin(node_folds[i]))]
        train_idx[i] = list(train_df['ID'])

    df.drop(columns='ID', inplace=True)

    return train_idx, val_idx




def split_edge_type_cv(df, n_folds, seed=None):
    '''
    Edge type appearing in training data will not appear in the test.
    '''
    df['ID'] = list(range(len(df)))

    edge_types = list(df['edge_type'].unique())# list of strings
    edge_types.sort()

    # now create 5 folds from the remaining drug and use each fold as val and the rest as train
    edge_type_folds = split_list(edge_types, n_folds, seed)

    df_split = {i: df[df['edge_type'].isin(edge_type_folds[i])]  for i in range(n_folds)}

    train_idx = {}
    val_idx = {}
    for i in range(n_folds):
        val_idx[i] = list(df_split[i]['ID'])
        # train_idx[i] = list(df['ID'].drop(val_idx[i]))
        train_idx[i] = list(df['ID'][~df['ID'].isin(val_idx[i])])

    df.drop(columns='ID', inplace=True)

    return train_idx, val_idx

def split_train_test(df, split_type, test_frac, seed=None):
    '''
    Inputs:
        df (dataframe): contains at least two columns named 'source' and 'target'. An optional column can be 'edge_type'.
                Depending on the presence of 'edge_type' a pair, e.g., (a,b) or a triplet, e.g., (a,b,C) can represent an edge.
                Here, a = source node, b = target node, C = edge type.
                There should not be any duplicated edge. If needed, use the remove_dups() function beforehand
                to remove duplicates.
        split_type (str): The values can be: 'random', 'node', 'edge', or 'edge_type'.
        test_frac (float): Fraction of  expected test samples.
        seed (int or None):
                - int will assure the same splits each time you split the data.
                - None will give different splits at each run.
    Return:
        1. train_idx(list): df.iloc[train_idx] consists the training dataset.
        2. test_idx(list): df.iloc[test_idx] consists the test dataset.
    '''
    split_type_2_function_map = {'random': split_random_train_test ,'node': split_node_train_test, 'edge': split_edge_train_test, 'edge_type': split_edge_type_train_test}
    train_idx, test_idx = split_type_2_function_map[split_type](df, test_frac, seed=seed)
    verify_split(df, train_idx, test_idx, split_type)
    return train_idx, test_idx

def split_cv(df, split_type, n_folds, seed=None):
    '''
    Inputs:
        df (dataframe): contains at least two columns named 'source' and 'target'. An optional column can be 'edge_type'.
                Depending on the presence of 'edge_type' a pair, e.g., (a,b) or a triplet, e.g., (a,b,C) can represent an edge.
                Here, a = source node, b = target node, C = edge type.
                There should not be any duplicated edge. If needed, use the remove_dups() function beforehand
                to remove duplicates.
        split_type (str): The values can be: 'random', 'node', 'edge', or 'edge_type'.
        n_folds (int): Number of folds in cross-validation.
        seed (int or None):
                - int will assure the same splits each time you split the data.
                - None will give different splits at each run.
    Return:
        1. train_idx(dict): df.iloc[train_idx[i]] consists the training dataset for ith fold.
        2. test_idx(dict): df.iloc[test_idx[i]] consists the test dataset for ith fold.
    '''
    split_type_2_function_map = {'random': split_random_cv ,'node': split_node_cv, 'edge': split_edge_cv, 'edge_type': split_edge_type_cv}
    train_idx, test_idx = split_type_2_function_map[split_type](df, n_folds, seed=seed)

    for i in range(n_folds):
        verify_split(df, train_idx[i], test_idx[i], split_type)
    return train_idx, test_idx

def generate_negative_samples_CNNC(df, graph_type='directed', anchor='source', seed=None):
    '''
    Parameters:
        graph_type: 'undirected', 'directed'
        anchor: 'source', 'target', 'both'

    If graph_type=='directed' and anchor=='source':
        - without edge_type: For any positive edge (a, b) create a negative edge (a, c) such that (a, c) was not present in the set of positive edges.
    '''
    
    #works for directed, anchor based graph without edgetype.
    if df.shape[1]>2:
        exit('Error: Not implemented for extra information on edges except for source and target. ')

    all_sampled_sources = []
    all_sampled_targets = []

    if (graph_type=='directed') and (anchor=='source'):
        init_sample_space = set(df['target'].unique())
        # print(df.head())
        # source_wise_targets = df.groupby('source').agg(target_list= ('target', lambda x:set(x)), count = ('target', 'size')).reset_index()
        source_wise_targets = df.groupby('source')['target'].agg([('target_list', lambda x:set(x)), ('count', 'size')]).reset_index() # TODO beeline compatibility for named aggregations not in place for pandas < 0.25
        source_wise_targets['target_list'] = source_wise_targets['target_list'].apply(lambda x: sorted(init_sample_space.difference(x)))
        # print(source_wise_targets)
        # print(source_wise_targets["count"].sum())

        for i, row in source_wise_targets.iterrows():
            sample_count = row['count']
            all_sampled_targets.extend(list(random.Random(seed).choices(row['target_list'], k=sample_count)))
            all_sampled_sources.extend([row['source']]*sample_count)
        negative_df = pd.DataFrame({'source': all_sampled_sources, 'target': all_sampled_targets})

        return negative_df
    else:
        exit('Error: current code only works for graph_type= directed, anchor=source.')







