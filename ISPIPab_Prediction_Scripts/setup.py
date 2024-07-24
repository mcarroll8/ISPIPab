from setuptools import setup

if __name__ == '__main__':
    with open("ReadMe.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()

    setup(name='ISPIP',
        version="1.15",
        description='Integrated Structure-based Protein Interface Prediction',
        url='https://github.com/eved1018/ISPIP',
        author='Evan Edelstein',
        author_email='edelsteinevan@gmail.com',
        license='MIT',
        packages=['ispip'],
        data_files=[('',["ReadMe.md"])],
        install_requires=[         
            'pandas',         
            'scikit-learn',
            'joblib',
            'argparse',
            "numpy",
            "scipy",
            "graphviz",
            "dtreeviz",
            "matplotlib"],
        long_description=long_description,
        long_description_content_type='text/markdown',
        entry_points = {
            'console_scripts': ['ispip=ispip.ispip:main'],
            
        })

# python setup.py sdist
# python3 -m twine upload dist/*