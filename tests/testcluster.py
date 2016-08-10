from igdiscover.cluster import inner_nodes


class Node:
	def __init__(self, value, left, right):
		self.value = value
		self.left = left
		self.right = right

	def is_leaf(self):
		return self.left is None and self.right is None


def test_inner_nodes():
	leaf = Node(0, None, None)
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
	inner = [ n.value for n in inner_nodes(tree) ]
	assert inner == [1, 2, 4, 5, 3, 6]
