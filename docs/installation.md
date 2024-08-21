# Installation


## Installing `prymer`

The installation requires three steps:

1. Install python and other dependencies with `conda`
2. Install the custom version of bwa
3. Install `prymer` with `poetry.

Install the required Python version, [`poetry`](https://github.com/python-poetry/poetry), and [`primer3](https://github.com/primer3-org/primer3) into your environment manager of choice, e.g.

```sh
$ mamba env create -y -f prymer.yml
$ conda activate prymer
```

Install the custom version of bwa:
```sh
$ git clone -b interactive_aln git@github.com:fulcrumgenomics/bwa.git
$ cd bwa
$ make -j 12
$ cp bwa ${CONDA_PREFIX}/bin
```

Note: the `virtualenvs.create false` setting in `poetry.toml` stops poetry from creating new virtual environments and forces it to use the active conda environment instead.
This can be set once per machine/user and stored in the user's poetry configuration with:

```sh
$ poetry config settings.virtualenvs.create false
```

Install the prymer with `poetry`.

```console
$ poetry install
```

## Getting Setup for Development Work

Follow the [instructions above](#installing-prymer)

```console
$ poetry install --with dev
```

## Checking the Build

Make sure that [instructions for development work](#getting-setup-for-development-work) have been followed. 

Use `poetry` to format, lint, type-check, and test your code.
Note that `poetry run pytest` will run `mypy` and `ruff` code checks in addition to `pytest` unit tests, and will provide a unit test coverage report.

```console
$ poetry run pytest 
```

However, `pytest` will neither run the ruff formatter nor apply `ruff`'s automatic lint fixes, which can be done by calling `ruff` directly. 

```console
$ poetry run ruff format && poetry run ruff check --fix
```

Static type checking is performed using `mpyp`.

```console
poetry run mypy
```

## Building the Documentation

Make sure that [instructions for development work](#getting-setup-for-development-work) have been followed.

Use `mkdocs` to build and serve the documentation.

```console
$ poetry install --with dev
$ poetry run mkdocs build
$ poetry run mkdocs serve
```

## Creating a Release on PyPi

1. Clone the repository recursively and ensure you are on the `main` (un-dirty) branch
2. Checkout a new branch to prepare the library for release
3. Bump the version of the library to the desired SemVer with `poetry version #.#.#`
4. Commit the version bump changes with a Git commit message like `chore(release): bump to #.#.#`
5. Push the commit to the upstream remote, open a PR, ensure tests pass, and seek reviews
6. Squash merge the PR
7. Tag the new commit on the main branch of the repository with the new SemVer

GitHub Actions will take care of the remainder of the deployment and release process with:

1. Unit tests will be run for safety-sake
2. A source distribution will be built
3. Many multi-arch multi-Python binary distributions will be built
4. Assets will be deployed to PyPi with the new SemVer
5. A [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/)-aware changelog will be drafted
6. A GitHub release will be created with the new SemVer and the drafted changelog

Consider editing the changelog if there are any errors or necessary enhancements.
