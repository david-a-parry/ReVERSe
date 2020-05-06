try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name="ReVERSe",
    packages=["ReVERSe"],
    version="0.1",
    description="Rare Variant Enrichment and Recessive Segregation",
    author="David A. Parry",
    author_email="david.parry@igmm.ed.ac.uk",
    url='https://github.com/david-a-parry/ReVERSe',
    download_url='https://github.com/david-a-parry/ReVERSe/archive/0.1.tar.gz',
    license='MIT',
    install_requires=['vase>=0.2.5', 'pandas', 'xlsxwriter', 'requests',
                      'biopython', 'scipy'],
    extras_require={
        'MYGENEINFO': ['mygene'],
    },
    scripts=["bin/ReVERSe_count",
             "bin/ReVERSe_reporter",
             "bin/ReVERSe_seg"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    dependency_links=[
        'https://github.com/david-a-parry/vase/tarball/master#egg=vase-0.2.5'],
)
