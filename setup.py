from setuptools import setup
from mixr.constants import VERSION
from distutils.extension import Extension


if __name__ == '__main__':
    setup(
        name='mixr',
        packages=['mixr'],
        version=VERSION,
        entry_points={
          'console_scripts': [
              'mixr = mixr.main:main'
          ]
        },
        include_package_data=True,
        zip_safe=False,
        data_files=[('resources', ['resources/pam30.txt'])],
        description='Mismatching Isoform Exon Remover',
        keywords=['DNA', 'bioinformatics', 'genomics', 'transcriptomics', 'alignment'],
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'License :: Freely Distributable',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
