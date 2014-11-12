#from distutils.core import setup
from setuptools import setup

descr = """lesion: quantitative analysis of spinal cord lesions in zebrafish.

This package contains IO, statistic, and plotting functions for this very
specific purpose. Functions will be backported to 
"""

DISTNAME            = 'lesion'
DESCRIPTION         = 'Analyse zebrafish lesion images'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'Juan Nunez-Iglesias'
MAINTAINER_EMAIL    = 'juan.n@unimelb.edu.au'
URL                 = 'https://github.com/jni/lesion'
LICENSE             = 'BSD 3-clause'
DOWNLOAD_URL        = 'https://github.com/jni/lesion'
VERSION             = '0.1-dev'
PYTHON_VERSION      = (2, 7)
INST_DEPENDENCIES   = {} 


if __name__ == '__main__': # pragma: no cover

    setup(name=DISTNAME,
        version=VERSION,
        url=URL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        license=LICENSE,
        packages=['lesion'],
        install_requires=INST_DEPENDENCIES,
        scripts=[]
    )

