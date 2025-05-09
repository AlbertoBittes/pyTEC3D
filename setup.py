from setuptools import setup, find_packages
import subprocess
import os

# Define the local_packages path
local_packages_path = os.path.join("api", "local_packages")

# Ensure the local_packages directory exists
os.makedirs(local_packages_path, exist_ok=True)

# Install openpyxl into the local_packages folder
subprocess.run(["pip", "install", "openpyxl", "--target", local_packages_path])

setup(
    name="PlasmaArc",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "vtk",
    ],
    description="A package for managing VTK unstructured grids and element types.",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/your-repo",
)
