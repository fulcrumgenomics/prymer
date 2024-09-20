# Installation and Developer's Documentation

## Recommended Installation

The package `prymer` requires installation of [Primer3](https://github.com/primer3-org/primer3) and [interactive `bwa`](https://github.com/fulcrumgenomics/bwa-aln-interactive).

To satisfy these requirements, it is recommended to install using [bioconda](https://bioconda.github.io/):

```console
mamba install -c bioconda prymer
```

## Installation for Development and Release

1. Install the environment manager [`mamba`](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
2. Install the Python build tool [`poetry`](https://python-poetry.org/docs/#installing-with-the-official-installer)
3. Create an environment with Python, [Primer3](https://github.com/primer3-org/primer3), and [interactive `bwa`](https://github.com/fulcrumgenomics/bwa-aln-interactive):

    ```console
    mamba env create -y -f prymer.yml
    ```

4. Activate the environment:

    ```console
    mamba activate prymer
    ```

5. Configure `poetry` to install into pre-existing virtual environments:

    ```console
    poetry config settings.virtualenvs.create false
    ```

6. Install `prymer` into the virtual environment:

    ```console
    poetry install
    ```

# Checking the Build

Use `poetry` to test your code.

```console
poetry run pytest
```

Note that `poetry run pytest` will run `mypy` checks, `ruff` checks, `pytest` unit tests, and will provide a unit test coverage report.

However, `pytest` will neither run the ruff formatter nor apply `ruff`'s automatic lint fixes, which can be done by calling `ruff` directly. 

```console
poetry run ruff format && poetry run ruff check --fix
```

# Building the Documentation

Use `mkdocs` to build and serve the documentation.

```console
poetry run mkdocs build && poetry run mkdocs serve
```

# Creating a Release on PyPi

1. Clone the repository recursively and ensure you are on the `main` (un-dirty) branch
2. Checkout a new branch to prepare the library for release
3. Bump the version of the library to the desired SemVer with `poetry version #.#.#`
4. Commit the version bump changes with a Git commit message like `chore(release): bump to #.#.#`
5. Push the commit to the upstream remote, open a PR, ensure tests pass, and seek reviews
6. Squash merge the PR
7. Tag the new commit on the main branch of the origin repository with the new SemVer

> [!NOTE]
> This project follows [Semantic Versioning](https://semver.org/).
> In brief:
> 
> * `MAJOR` version when you make incompatible API changes
> * `MINOR` version when you add functionality in a backwards compatible manner
> * `PATCH` version when you make backwards compatible bug fixes

GitHub Actions will take care of the remainder of the deployment and release process with:

1. Unit tests will be run for safety-sake
2. A source distribution will be built
3. Multi-arch multi-Python binary distributions will be built
4. Assets will be deployed to PyPi with the new SemVer
5. A [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/)-aware changelog will be drafted
6. A GitHub release will be created with the new SemVer and the drafted changelog

> [!IMPORTANT]
> Consider editing the changelog if there are any errors or necessary enhancements.
