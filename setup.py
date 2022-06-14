import setuptools


setuptools.setup(
    name="LiftoffTools",
    version="0.2.0",
    author="Alaina Shumate",
    author_email="ashumat2@jhmi.edu",
    description="A tool for comparing annotations across genome assemblies",
    url="https://github.com/ashumate/Planet",
    install_requires=[],
    python_requires='>=3.6',
    packages=['liftofftools'],
    entry_points={'console_scripts': ['liftofftools = liftofftools.liftofftools:main'], },
)
