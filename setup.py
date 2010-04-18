from setuptools import setup, find_packages

from strprofiles import __version__ as VERSION


setup(
    name='strprofiles',
    version=VERSION,
    description='An example github python project.',
    author='Martin Budden',
    author_email='',
    platforms='Posix; MacOS X; Windows',
    scripts=['example_script'],
    packages=find_packages(exclude=['test']),
    install_requires=[
        'example_require'],
    include_package_data=True,
    zip_safe=False)
