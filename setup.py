
from setuptools import setup

setup(name='pppmodel',
      version = '0.1',
      description = 'model of the pentose phosphate pathway from McIntyre et al. 1989, later adapted by Berthon et al. 1993',
      url='https://github.com/QTB-HHU/ppp-modelbase',
      author = 'Tim Nies',
      author_email = 'Tim.Nies@hhu.de',
      license = 'GPL3',
      packages = ['pppmodel'],
      intsall_require=[
          'modelbase',
          'scipy',
          'numpy',
          'numdifftools'
          'matplotlib'
      ],
      zip_safe=False)