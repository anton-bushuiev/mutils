import math

from mutils.mutations import MutationSpace


def test_size():
    # Simple case where numbers of candidates are equal on each position
    dct = {
        'YC17': 'AG',
        'TA20': 'AD',
        'AC18': 'AG',
        'AC11': 'NE',
    }
    space = MutationSpace(dct)
    for d in range(1, len(dct)+1):
        expected = math.comb(len(dct), d) * len(list(dct.values())[0])**d
        real = space.size(d=d)
        assert expected == real

    # "Difficult case"
    dct = {
        'YC17': 'AG',
        'TA20': 'ADG',
        'AC11': 'NE',
        'AC18': ''
    }
    space = MutationSpace(dct)
    assert 7 == space.size(d=1)
    assert 16 == space.size(d=2)
    assert 12 == space.size(d=3)
    assert 0 == space.size(d=4)

    # Zero and wt work
    assert 0 == space.size(d=0)
    assert 1 == space.size(d=0, wt=True)
    assert 13 == space.size(d=3, wt=True)

    # Total works
    assert 7 + 16 + 12 + 0 == space.size()


def test_construct():
    dct = {
        'YC17': 'AG',
        'TA20': 'AD',
        'AC18': 'AGE',
        'AC11': 'NE',
        'AE1': ''
    }
    space = MutationSpace(dct)

    # Size is correct
    assert space.size() == len(space.construct())
    for d in range(1, space.n_pos + 1):
        assert space.size(d=d) == len(space.construct(d=d))

    # Format is correct
    assert [] == space.construct(d=0)
    for d in range(1, space.n_pos + 1):
        mutations = space.construct(d=d)
        for mut in mutations:
            assert d == mut.count(',') + 1


def test__dict_to_list():
    retval = MutationSpace._dict_to_list({
        'YC17': 'AG',
        'TA20': 'A',
        '5C': ''
    })
    expected = [['YC17A', 'YC17G'], ['TA20A'], []]
    assert expected == retval
