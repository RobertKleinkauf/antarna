from setuptools import setup, find_packages

setup(
    name="antarna",
    packages=['antarna'],
    install_requires=['numpy', 'scipy', 'argparse', 'networkx', 'uuid'],
    version="2.0.1",
    description="Ant colony optimized RNA secondary structure inverse folding",
    author='Robert Kleinkauf',
    author_email='robkleinkauf@gmail.com',
    url='https://github.com/RobertKleinkauf/antarna',
    include_package_data=True,
    keywords=['ant colony optimization', 'RNA inverse folding'],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Scientists',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
    ]
)