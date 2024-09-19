import os
from pathlib import Path
from tempfile import NamedTemporaryFile
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


def test_validate_executable_path_does_not_exist() -> None:
    """
    `validate_executable_path` should raise a ValueError when the provided executable does not
    exist.
    """
    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable="/path/to/nowhere")


def test_validate_executable_path_not_executable() -> None:
    """
    `validate_executable_path` should raise a ValueError when the provided executable does not have
    execute permissions.
    """
    with NamedTemporaryFile(suffix=".exe", mode="w", delete=True) as tmpfile:
        with pytest.raises(ValueError, match="is not executable"):
            ExecutableRunner.validate_executable_path(executable=tmpfile.name)


@pytest.mark.parametrize("executable", ["/usr/bin/yes", Path("/usr/bin/yes")])
def test_validate_executable_path_returns_existing_paths(executable: str | Path) -> None:
    """
    `validate_executable_path` should return the `yes` executable when the absolute path to the
    executable is passed, either as a string or a Path.
    """
    expected_path: Path = Path("/usr/bin/yes")

    # Clear the PATH, in case the user has a local version of `yes` elsewhere on their PATH
    with mock.patch.dict(os.environ, clear=True):
        os.environ["PATH"] = "/usr/bin"

        validated_path: Path = ExecutableRunner.validate_executable_path(executable=executable)
        assert validated_path == expected_path


def test_validate_executable_finds_on_path() -> None:
    """
    `validate_executable_path` should find executables on the PATH when the name of the executable
    is passed as a string.
    """
    expected_path: Path = Path("/usr/bin/yes")

    with mock.patch("shutil.which", return_value=expected_path.as_posix()):
        validated_path: Path = ExecutableRunner.validate_executable_path(executable="yes")
        assert validated_path == expected_path


def test_validate_executable_path_rejects_paths() -> None:
    """
    `validate_executable_path` should not treat non-existent Path objects as valid executables.

    Specifically, if the user passes the name of an executable on the PATH as a `Path` instead of a
    string, it should be treated as a non-existent Path and a `ValueError` should be raised.
    """
    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable=Path("yes"))


def test_validate_executable_path_new_file() -> None:
    with NamedTemporaryFile(suffix=".exe", mode="w", delete=True) as tmpfile:
        exec_str: str = tmpfile.name
        exec_path: Path = Path(exec_str)
        # not an executable
        with pytest.raises(ValueError, match="is not executable"):
            ExecutableRunner.validate_executable_path(executable=exec_str)
        # make it executable and test again
        os.chmod(exec_str, 755)
        exec_full_path: Path = exec_path.absolute()
        assert exec_full_path == ExecutableRunner.validate_executable_path(executable=exec_str)
        assert exec_full_path == ExecutableRunner.validate_executable_path(executable=exec_path)
        assert exec_full_path == ExecutableRunner.validate_executable_path(
            executable=exec_full_path
        )
