from setuptools import find_packages
from setuptools import setup

classifiers = [
  "Development Status :: 3 - Alpha ",
  "Intended Audience :: Education",
  "License :: OSI Approved :: MIT License",
  "Programming Language :: Python :: 3"
]

setup(
	name='residence',
	version='1.0.0',
	description='residence behavior codes',
	author = 'Zaf4',
	author_email = 'zafermolbio@gmail.com'
	url = 'https://github.com/Zaf4/residence',
	packages = find_packages()

)