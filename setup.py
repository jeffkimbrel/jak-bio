import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="jak-bio",
    version="0.0.1",
    author="Jeff Kimbrel",
    author_email="jakpot@gmail.com",
    description="Various omics tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['jakomics'],
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    python_requires='>=3.6',
)