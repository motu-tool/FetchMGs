from setuptools import setup

setup(
    name='fetchmgs',
    version='1.2',
    description='FetchMGs extracts the 40 marker genes from genomes and metagenomes in an easy and accurate manner.',
    url='https://github.com/SushiLab/fetchmgs_dev',
    author='Chris Field, Shinichi Sunagawa, Daniel R. Mende',
    author_email='fieldc@ethz.ch',
    license='GPL-3.0',
    packages=setuptools.find_packages(),
    package_data={'fetchmgs': ['data/*']},
    entry_points = {
        'console_scripts': ['fetchmgs=fetchmgs.__main__:main'],
    }
)

