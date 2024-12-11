from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
install_requires = ['biopython ==1.84',
                    'tqdm',
                    'pyrodigal ==3.6.3',
                    'pyhmmer ==0.10.15',
                    'psutil']
long_description = read('README.md')

setup(
    name='FetchMGs',
    version='2.0.0',
    description='FetchMGs extracts the 40 marker genes from genomes and metagenomes in an easy and accurate manner.',
    url='https://github.com/motu-tool/FetchMGs',
    author='Hans-Joachim Ruscheweyh, Chris Field, Shinichi Sunagawa, Daniel R. Mende',
    author_email='hansr@ethz.ch',
    license='GPL-3.0',
    include_package_data=True,
    install_requires=install_requires,
    long_description=long_description,
    packages=['fetchmgs'],
    download_url = "https://github.com/motu-tool/FetchMGs/archive/refs/tags/2.0.0.tar.gz",
    entry_points = {
        'console_scripts': ['fetchMGs=fetchmgs.fetchmgs:main'],
        },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        ]
    
)
