import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open('requirements.txt', 'r') as file:
    requirements = file.read().splitlines()

PACKAGE_DIR = {"": "src"}
PACKAGES = setuptools.find_packages(where="src")
setuptools.setup(
    name="hyp_solver",
    version="0.0.2",
    author="Dunn Kopylov",
    author_email="38dunn@gmail.com",
    description="hyp solver",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    package_dir=PACKAGE_DIR,
    packages=PACKAGES,
    python_requires=">=3.9"
)