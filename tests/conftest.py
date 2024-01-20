import pytest

from mutils.data import load_SKEMPI2


@pytest.fixture(scope='session')
def skempi2():
    return load_SKEMPI2()[0]


@pytest.fixture(scope='session')
def skempi2_multimeric_partners(skempi2):
    return skempi2[skempi2['Partner 1'].str.len() > 1]


@pytest.fixture(scope='session')
def skempi2_big_multimeric_partners(skempi2):
    return skempi2[skempi2['Partner 1'].str.len() > 2]
