import pytest

from mutils.mutations import Mutation


def test_str():
    mut = 'YC17T,TA20A'
    assert mut == str(Mutation.from_str(mut))

    mut_empty = Mutation()
    assert str(mut_empty) == ''


def test_revert():
    assert 'TC17Y' == str(Mutation.from_str('YC17T').revert())

    mut = Mutation.from_str('TC17Y,AA20T')
    assert mut == mut.revert().revert()

