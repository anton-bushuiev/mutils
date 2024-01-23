from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='mutils',
    version='0.1',
    packages=find_packages(),
    package_data={'mutils': ['data/*']},
    install_requires=required,
)
