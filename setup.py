import setuptools


setuptools.setup(
    name="LiftoffTools",
    version="0.1.0",
    author="Alaina Shumate",
    author_email="ashumat2@jhmi.edu",
    description="A tool for comparing annotations across genome assemblies",
    url="https://github.com/ashumate/Planet",
    install_requires=['nltk==3.6.7','numpy==1.21.1', 'parasail==1.2.4', 'gffutils==0.10.1', 'pyfaidx==0.5.8', 'biopython==1.76'],
    python_requires='>=3.6',
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': ['liftofftools = liftofftools.liftofftools:main'], },
)
