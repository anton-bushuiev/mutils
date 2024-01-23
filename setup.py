from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='mutils',
    version='0.1',
    packages=find_packages(),
    package_data={'': ['data/*']},
    include_package_data=True,
    install_requires=required,
    # package_dir={"": ""}
)
