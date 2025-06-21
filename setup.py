from setuptools import setup


install_requires = ['biopython',
                    'tqdm',
                    'pyrodigal ==3.6.3',
                    'pyhmmer ==0.10.15',
                    'psutil']


with open("README.md", "r") as fh:
    long_desc = fh.read()
    long_desc = "\n".join(long_desc.split("\n")[1:])

setup(
    name='fetchMGs',
    version='2.1.0',
    description='FetchMGs extracts the 40 marker genes from genomes and metagenomes in an easy and accurate manner.',
    url='https://github.com/motu-tool/FetchMGs',
    author='Hans-Joachim Ruscheweyh, Chris Field, Shinichi Sunagawa, Daniel R. Mende',
    author_email='hansr@ethz.ch',
    license='GPL-3.0',
    include_package_data=True,
    install_requires=install_requires,
    long_description=long_desc,
    long_description_content_type = "text/markdown",
    packages=['fetchmgs'],
    download_url = "https://github.com/motu-tool/FetchMGs/archive/refs/tags/2.1.0.tar.gz",
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
