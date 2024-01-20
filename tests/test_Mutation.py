import pytest

from mutils.mutations import Mutation


def test_constructor():
    # Different input formats work
    mut_from_str = Mutation('YC17T,TA20A')
    mut_from_list = Mutation(muttuple=['YC17T', 'TA20A'])
    mut_from_tuple = Mutation(muttuple=('YC17T', 'TA20A'))
    assert mut_from_str.__dict__ == mut_from_list.__dict__ ==\
           mut_from_tuple.__dict__

    # Empty mutation is OK
    mut_empty = Mutation()
    assert mut_empty.muttuple == ()

    # Both formats lead to error
    with pytest.raises(ValueError):
        Mutation('YC17T,TA20A', ('YC17T', 'TA20A'))


def test_str():
    mut = 'YC17T,TA20A'
    assert mut == str(Mutation(mut))


def test_revert():
    assert 'TC17Y' == str(Mutation('YC17T').revert())
    assert 'TC17Y,AA20T' == str(Mutation('YC17T,TA20A').revert())


def test_rename_chains():
    assert 'YA17T' == str(Mutation('YC17T').rename_chains({'C': 'A'}))
    assert 'YA17T,TD20A' == str(Mutation('YC17T,TA20A').
                                rename_chains({'C': 'A', 'A': 'D'}))
