from setuptools import setup, find_packages

setup(
    name='vira-av',
    version='1.0.0',
    author='Ales Varabyou',
    author_email='ales.varabyou@jhu.edu',
    url='https://github.com/alevar/vira',
    description='VIRA: By-Reference Exon and CDS Viral Genome Annotation',
    license='GPLv3',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'intervaltree',
        'numpy',
        'pyfaidx',
        'pysam',
        'setuptools',
        'snapper'
    ],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'vira=vira.core:main',
        ],
    },
)