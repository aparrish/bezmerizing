from setuptools import setup
    
with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='bezmerizing',
    version='0.0.1',
    author='Allison Parrish',
    author_email='allison@decontextualize.com',
    url='https://github.com/aparrish/bezmerizing',
    description='tiny quirky library with some bezier curve functions',
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=['flat>=0.0.3', 'scipy>=1.1.0'],
    dependency_links=[],
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    platforms='any',
)
