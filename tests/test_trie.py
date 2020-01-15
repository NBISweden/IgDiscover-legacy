import random
from igdiscover.trie import Trie
from igdiscover.cli.group import hamming_neighbors
from tinyalign import hamming_distance


def random_nt(length):
    return ''.join(random.choice('ACGT') for _ in range(length))


class TestTrie:
    LENGTHS = (0, 1, 2, 3, 4, 5, 6, 10, 12, 15, 20)

    def setup(self):
        self.strings = set()
        self.trie = Trie()
        # Create a set and a Trie both containing the same strings
        for _ in range(80):
            for length in self.LENGTHS:
                s = random_nt(length)
                self.strings.add(s)
                self.trie.add(s)

    def test_empty_string(self):
        trie = Trie()
        trie.add('')
        assert '' in trie
        assert 'A' not in trie

    def test_contains(self):
        for s in self.strings:
            assert s in self.trie

        for length in self.LENGTHS:
            for _ in range(min(100, 4**length)):
                s = random_nt(length)
                assert (s in self.strings) == (s in self.trie)

    def test_len(self):
        assert len(self.strings) == len(self.trie)

    def naive_has_similar(self, s, distance):
        for t in self.strings:
            if len(t) != len(s):
                continue
            if hamming_distance(t, s) <= distance:
                return True
        return False

    def naive_find_all_similar(self, s, distance):
        for t in self.strings:
            if len(t) != len(s):
                continue
            if hamming_distance(t, s) <= distance:
                yield t

    def test_has_similar(self):
        for s in self.strings:
            assert self.trie.has_similar(s, 0)
            assert self.trie.has_similar(s, 1)
            assert self.trie.has_similar(s, 2)

        for base in self.strings:
            for modified in hamming_neighbors(base):
                assert self.trie.has_similar(modified, 1)

        for errors in range(4):
            for length in self.LENGTHS:
                for _ in range(min(100, 4**length)):
                    s = random_nt(length)
                    assert self.naive_has_similar(s, errors) == self.trie.has_similar(s, errors)

    def test_find_all_similar(self):
        t = Trie()
        t.add('ACGT')
        result = list(t.find_all_similar('ACGG', 1))
        assert set(result) == frozenset(('ACGT',))

        for s in self.strings:
            result = list(self.trie.find_all_similar(s, 0))
            assert result == [s]

            for errors in range(1, 4):
                result = list(self.trie.find_all_similar(s, errors))
                assert s in result

        for base in self.strings:
            for modified in hamming_neighbors(base):
                assert base in self.trie.find_all_similar(modified, 1)

        for errors in range(1, 3):
            for length in self.LENGTHS:
                for _ in range(min(100, 4**length)):
                    s = random_nt(length)
                    expected = set(self.naive_find_all_similar(s, errors))
                    assert expected == set(self.trie.find_all_similar(s, errors))


def main():
    import sys
    n = int(sys.argv[1])
    dist = int(sys.argv[2])

    strings = [random_nt(13) for _ in range(n)]
    print('created random sequences')

    if sys.argv[3] == 'x':
        def naive_has_similar(t):
            for s in strings:
                if hamming_distance(s, t) <= dist:
                    return True
            return False

        hs = 0
        for s in strings:
            hs += int(naive_has_similar(s[1:] + s[0]))
        print('hs:', hs, 'out of', len(strings))

    else:
        trie = Trie()
        for s in strings:
            trie.add(s)
        print('created trie')
        hs = 0
        for s in strings:
            hs += int(trie.has_similar(s[1:] + s[0], dist))

        print('hs:', hs, 'out of', len(strings))


if __name__ == '__main__':
    main()
