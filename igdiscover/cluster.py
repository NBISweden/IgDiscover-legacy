import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy
from .utils import distances, iterative_consensus


def inner_nodes(root):
	"""
	Return a list of all inner nodes of the tree
	"""
	if root.is_leaf():
		return []
	return inner_nodes(root.left) + [root] + inner_nodes(root.right)


def collect_ids(root):
	"""
	Return a list of ids of all leaves of the given tree
	"""
	if root.is_leaf():
		return [root.id]
	return collect_ids(root.left) + collect_ids(root.right)


def cluster_sequences(sequences, minsize=5):
	"""
	Cluster the given sequences into groups of similar sequences.

	Return a triple that contains a pandas.DataFrame with the edit distances,
	the linkage result, and a list that maps sequence ids to their cluster id.
	If an entry is zero in that list, it means that the sequence is not part of
	a cluster.
	"""
	matrix = distances(sequences)
	linkage = hierarchy.linkage(distance.squareform(matrix), method='average')
	# Linkage columns are:
	# 0, 1: merged clusters, 2: distance, 3: number of nodes in cluster
	inner = inner_nodes(hierarchy.to_tree(linkage))
	prev = linkage[:,2].max()  # highest distance
	clusters = [0] * len(sequences)
	cl = 1
	for n in inner:
		if prev/n.dist < 0.8 and n.left.count >= minsize and n.right.count >= minsize:
			for id in collect_ids(n.left):
				# Do not overwrite previously assigned ids
				if clusters[id] == 0:
					clusters[id] = cl
			cl += 1
		prev = n.dist

	return pd.DataFrame(matrix), linkage, clusters
