from setuptools import setup, find_packages

setup(
    name='fluxpal',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'astropy',
        'seaborn',
        'numpy'
    ],
    author='Bella Longo',
    author_email='bellalongo.mail@gmail.com',
    description='Given a FITS spectrum file, this script fits a Voigt profile to each emission line, facilitating flux calculation and data validity assessment. It generates plots for each line and outputs comprehensive results in both FITS and ECSV formats. These files contain all necessary information for analysis or reproduction.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/bellalongo/fluxpal',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],
)
