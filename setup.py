from setuptools import setup, find_packages

# Read the README file
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="coastal_fim_vis",  
    version="0.1.0",  
    author="Arash Modaresi Rad",
    author_email="arash.mod.rad@noaa.gov",
    description="A package for barycentric interpolation of Schisim outputs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lynker-spatial/coastal-fim.git",  
    packages=find_packages(),  
    include_package_data=True,  # To include non-code files
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',  
    install_requires=[],  
    extras_require={
        # Optional dependencies 
        "dev": [
            "pytest",
            "flake8",
        ]
    },
)
