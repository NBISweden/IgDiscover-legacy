import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from .utils import distances


def all_nodes(root):
	"""Return a list of all nodes of the tree, from left to right (iterative implementation)."""
	result = []
	path = [None]
	node = root
	while node is not None:
		if node.left is not None:
			path.append(node)
			node = node.left
		elif node.right is not None:
			result.append(node)
			node = node.right
		else:
			result.append(node)
			node = path.pop()
			if node is not None:
				result.append(node)
				node = node.right

	return result


def inner_nodes(root):
	"""
	Return a list of all inner nodes of the tree, from left to right.
	"""
	return [node for node in all_nodes(root) if not node.is_leaf()]


def collect_ids(root):
	"""
	Return a list of ids of all leaves of the given tree
	"""
	return [node.id for node in all_nodes(root) if node.is_leaf()]


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
	prev = linkage[:, 2].max()  # highest distance
	clusters = [0] * len(sequences)
	cl = 1
	for n in inner:
		if n.dist > 0 and prev / n.dist < 0.8 \
				and n.left.count >= minsize and n.right.count >= minsize:
			for id in collect_ids(n.left):
				# Do not overwrite previously assigned ids
				if clusters[id] == 0:
					clusters[id] = cl
			cl += 1
		prev = n.dist
	# At the end of the above loop, we have not processed the rightmost
	# subtree. In our experiments, it never contains true novel sequences,
	# so we omit it.

	return pd.DataFrame(matrix), linkage, clusters
