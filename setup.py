#!/usr/bin/env python

# import os
import sys
from setuptools import setup
# from pkg_resources import resource_filename
from subprocess import call as subpcall
# from setuptools import find_packages

if sys.version < "2.6.0" or sys.version > "2.8.0":
    print "Please use a Python with higher version than 2.6.0"
    sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
    sys.exit(1)
    # exit(1)

def run_cmd(command):
    subpcall(command, shell=True)

def main():
    if not float(sys.version[:3]) >= 2.5:
        sys.stderr.write(
            "CRITICAL: Python version must be greater than or equal to 2.5! python 2.7.1 or newer is recommended!\n")
        sys.exit(1)
    setup(name="Cisrome-GO-Package",
          version="1.0.0",
          description="Cistrome-GO -- functional enrichment analysis of transcription factor ChIP-seq peaks ",
          license='MIT',
          author='Changson Wan',
          author_email='wchangson@gmail.com',
          package_dir={'Cisrome-GO': 'cistromego'},
          install_requires=['xlmhg>=2.4.9', 'mne>=0.17.0'],
          packages=['cistromego'],
          scripts=['bin/cistromego'],
          package_data={'cistromego': ['data/*']},

          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Environment :: Web Environment',
              'Intended Audience :: Developers',
              'License :: OSI Approved :: Artistic License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              'Topic :: Database',
          ],
          )


if __name__ == '__main__':
    main()