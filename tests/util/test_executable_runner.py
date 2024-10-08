import os
from pathlib import Path
from unittest import mock

import pytest

from prymer.util import ExecutableRunner


def test_no_command() -> None:
    with pytest.raises(ValueError, match="Invocation must not be empty"):
        ExecutableRunner(command=[])


def test_close_twice() -> None:
    exec = ExecutableRunner(command=["sleep", "5"])
    assert exec.close() is True
    assert exec.close() is False


def test_validate_executable_path_does_not_exist(tmp_path: Path) -> None:
    """
    `validate_executable_path` should raise a ValueError when the provided executable does not
    exist.
    """
    bad_path: Path = tmp_path / "nowhere"

    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable=bad_path)


def test_validate_executable_path_not_executable(tmp_path: Path) -> None:
    """
    `validate_executable_path` should raise a ValueError when the provided executable does not have
    execute permissions.
    """
    bad_path: Path = tmp_path / "not_executable"
    bad_path.touch()

    with pytest.raises(ValueError, match="is not executable"):
        ExecutableRunner.validate_executable_path(executable=bad_path)


def test_validate_executable_path(tmp_path: Path) -> None:
    """
    `validate_executable_path` should return the `yes` executable in the following scenarios:
    1. The name of the executable is passed as a string, and the executable is on the user's PATH.
    2. The absolute path to the executable is passed, either as a string or a Path.
    """
    expected_path = tmp_path / "yes"
    expected_path.touch()
    expected_path.chmod(755)

    # Clear the PATH, to override any local versions of `yes` on the user's PATH
    with mock.patch.dict(os.environ, clear=True):
        os.environ["PATH"] = str(tmp_path)

        executables: list[str | Path] = ["yes", expected_path, str(expected_path)]
        for executable in executables:
            validated_path: Path = ExecutableRunner.validate_executable_path(executable=executable)
            assert validated_path == expected_path


def test_validate_executable_path_rejects_paths() -> None:
    """
    `validate_executable_path` should not treat non-existent Path objects as valid executables.

    Specifically, if the user passes the name of an executable on the PATH as a `Path` instead of a
    string, it should be treated as a non-existent Path and a `ValueError` should be raised.
    """
    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable=Path("yes"))
