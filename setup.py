from setuptools import setup
import versioneer

def readme():
    with open('README.md') as f:
        return f.read()
    
setup(name='wavegraph',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Wavegraph -- graph-based clustering for GW search with coherent Waveburst',
      long_description=readme(),
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      keywords='gravitational waves',
      url='http://github.com/ecm0/wavegraph',
      author='Eric Chassande-Mottin',
      author_email='ecm@apc.in2p3.fr',
      license='LGPL-V3',
      packages=['wavegraph'],
      install_requires=[
          'numpy', 'lalsuite', \
          'matplotlib'
      ],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      #entry_points={
      #    'console_scripts': ['funniest-joke=funniest.command_line:main'],
      #},
      include_package_data=True,
      zip_safe=False)
