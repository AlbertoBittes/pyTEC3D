from setuptools import setup, find_packages
import subprocess
import os

# Define the local_packages path
local_packages_path = os.path.join("api", "local_packages")

# Ensure the local_packages directory exists
os.makedirs(local_packages_path, exist_ok=True)

print("Local installation of openpyxl for ParaView...\n")

# Install openpyxl into the local_packages folder
subprocess.run(["pip", "install", "openpyxl", "--target", local_packages_path])

setup(
    name="apec_m",
    version="0.1.0",
    packages=find_packages(where="api"),  # Look for packages in the 'api' directory
    package_dir={"": "api"},  # Root package directory is 'api'
    install_requires=[
        "pandas",
        "vtk",
        "numpy",
        "ansys.mapdl.reader",  # Dependency for convert_mesh.py
    ],
    entry_points={
        "console_scripts": [
            "apec=apec_m:main",  # CLI entry point for apec_m
            "convert_mesh=convert_mesh:main",  # CLI entry point for convert_mesh
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
