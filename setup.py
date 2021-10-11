from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(
        name="pykna",
        version='0.1',
        description='This project summarizes remote sensing data analysis tools',
        author='Kazuki Nakata',
        author_email='knakata0924@gmail.com',
        url='https://github.com/kazuki-nakata/',
        packages=find_packages(),
#        install_requires=open('requirements.txt').read().splitlines()
    )