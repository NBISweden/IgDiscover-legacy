from igdiscover.cluster import inner_nodes


def inner_nodes_recursive(root):
	"""
	Return a list of all inner nodes of the tree, from left to right.
	"""
	if root.is_leaf():
		return []
	return inner_nodes_recursive(root.left) + [root] + inner_nodes_recursive(root.right)


def collect_ids_recursive(root):
	"""
	Return a list of ids of all leaves of the given tree
	"""
	if root.is_leaf():
		return [root.id]
	return collect_ids_recursive(root.left) + collect_ids_recursive(root.right)


def test_inner_nodes():
	class Node:
		def __init__(self, value, left, right):
			self.value = value
			self.left = left
			self.right = right

		def is_leaf(self):
			return self.left is None and self.right is None

		def __repr__(self):
			# return 'Node({!r}, {!r}, {!r})'.format(self.value, self.left, self.right)
			return 'Node({!r}, ...)'.format(self.value)

	empty_tree = None
	assert inner_nodes(empty_tree) == []

	leaf = Node(0, None, None)
	assert inner_nodes(leaf) == []

	tree = Node(1, leaf, leaf)
	values = [ n.value for n in inner_nodes(tree) ]
	assert values == [1]

	tree = Node(1, Node(2, Node(3, leaf, leaf), leaf), leaf)
	values = [ n.value for n in inner_nodes(tree) ]
	assert values == [3, 2, 1]

	tree = Node(1, Node(2, Node(3, leaf, Node(4, leaf, leaf)), leaf), leaf)
	values = [ n.value for n in inner_nodes(tree) ]
	assert values == [3, 4, 2, 1]

	tree = Node(1,
		leaf,
		Node(2,
			leaf,
			Node(3,
	            Node(4,
	                leaf,
	                Node(5, leaf, leaf)
	            ),
	            Node(6, leaf, leaf)
			)
		)
	)
	values = [ n.value for n in inner_nodes(tree) ]
	assert values == [1, 2, 4, 5, 3, 6]
