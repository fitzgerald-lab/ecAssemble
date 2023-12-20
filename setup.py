from setuptools import setup

setup(
    name='ecAssemble',
    version='0.1',
    packages=['matplotlib'],
    entry_points={
        'console_scripts': [
            'ecAssemble = ecAssemble.run_ecAssemble:main',
        ],
    },
)
