"""
This trie can be used to store a set of strings and to retrieve sequences similar to a query
sequence. 
"""


class Trie:
    """
    A tree-like datastructure for storing strings. It supports queries similar to a set() of
    strings. This particular implementation also supports searching for strings by Hamming distance.
    """
    # Children
    A = None
    C = None
    G = None
    T = None
    # The name attribute is set to the string that is spelled from the root of the tree to this
    # node - but only for leaf nodes.
    name = None

    def __init__(self, iterable=None):
        if iterable is not None:
            for it in iterable:
                self.add(it)

    def add(self, s: str):
        """Add a string to this trie. If it already exists, the trie remains unchanged."""
        self._insert(s, leaf_name=s)

    def _insert(self, s: str, leaf_name: str):
        """Recursive insert, called by add()"""
        if len(s) == 0:
            # This needs to become a leaf node
            assert self.name is None or self.name == leaf_name
            self.name = leaf_name
        else:
            subtrie = getattr(self, s[0])
            if subtrie is None:
                subtrie = Trie()
                setattr(self, s[0], subtrie)
            subtrie._insert(s[1:], leaf_name)

    def __repr__(self):
        parts = []
        for c in list('ACGT'):
            if getattr(self, c) is not None:
                if getattr(self, c) is None:
                    parts.append(c)
                else:
                    parts.append(c + ':' + repr(getattr(self, c)))
        if self.name is not None:
            parts.append('Leaf({})'.format(self.name))
        if parts:
            return '{' + ', '.join(parts) + '}'
        else:
            return '*'

    def find_node(self, s: str):
        """
        Search for s. If s is a prefix of any string in the trie, then return the node
        corresponding to that prefix. The node can be an internal or a leaf node.
        A leaf node is a node in which name attribute is not None. If the returned node
        is a leaf node, then s is in the trie. If it is an internal node (that is,
        the name attribute is None), then there exists a string in the trie that has
        s as a proper prefix.
        
        If no such string (with s as prefix) exists, then None is returned.
        """
        if len(s) == 0:
            return self
        subtrie = getattr(self, s[0])
        if subtrie is None:
            return None
        else:
            return subtrie.find_node(s[1:])

    def __contains__(self, s: str):
        """Return whether s is in the trie"""
        node = self.find_node(s)
        return node is not None and node.name == s

    def count_nodes(self, internal=True):
        """Return number of nodes in this trie. internal nodes are counted if "internal" is True"""
        n = 1 if internal or self.name is not None else 0
        for c in list('ACGT'):
            subtrie = getattr(self, c)
            if subtrie is not None:
                n += len(subtrie)
        return n

    def __len__(self):
        """Return the number of unique strings in this trie"""
        return self.count_nodes(internal=False)

    def has_similar(self, s, mismatches):
        """
        Return whether a string exists in this trie that has a Hamming distance
        of at most 'mismatches' to s.
        """
        if len(s) == 0:
            return self.name is not None
        # As a runtime heuristic, descend into the the mismatches==0 subtrie first
        # if possible.
        subtrie = getattr(self, s[0])
        if subtrie is not None and subtrie.has_similar(s[1:], mismatches):
            return True
        if mismatches == 0:
            return False
        # The above did not work - try all three possible substitutions
        # and recursively check subtries
        for c in list('ACGT'):
            if c == s[0]:
                continue
            subtrie = getattr(self, c)
            if subtrie is None:
                continue
            if subtrie.has_similar(s[1:], mismatches - 1):
                return True
        return False

    def find_all_similar(self, s, mismatches):
        """
        Yield all strings in this trie that have a Hamming distance of at most
        'mismatches' to s.
        """
        # This routine is similar to has_similar, but since we are interested
        # in all similar strings, the optimization of going into the matching
        # subtrie first does not help.
        if len(s) == 0:
            if self.name is not None:
                yield self.name
            return

        # The code below is an optimized version of the following:
        #
        # for c in list('ACGT'):
        #   if c != s[0] and mismatches == 0:
        #     continue
        #   subtrie = getattr(self, c)
        #   if subtrie is not None:
        #     yield from subtrie.find_all_similar(s[1:], mismatches - int(c != s[0]))

        if mismatches > 0:
            s0 = s[0]
            s1 = s[1:]
            subtrie = self.A
            if subtrie is not None:
                yield from subtrie.find_all_similar(s1, mismatches - (int(s0 != 'A')))
            subtrie = self.C
            if subtrie is not None:
                yield from subtrie.find_all_similar(s1, mismatches - (int(s0 != 'C')))
            subtrie = self.G
            if subtrie is not None:
                yield from subtrie.find_all_similar(s1, mismatches - (int(s0 != 'G')))
            subtrie = self.T
            if subtrie is not None:
                yield from subtrie.find_all_similar(s1, mismatches - (int(s0 != 'T')))
        else:
            subtrie = self
            for c in s:
                subtrie = getattr(subtrie, c)
                if subtrie is None:
                    break
            else:
                if subtrie.name is not None:
                    yield subtrie.name
