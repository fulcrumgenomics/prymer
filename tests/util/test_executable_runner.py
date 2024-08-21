import os
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest

from prymer.util import ExecutableRunner


def test_no_command() -> None:
    with pytest.raises(ValueError, match="Invocation must not be empty"):
        ExecutableRunner(command=[])


def test_close_twice() -> None:
    exec = ExecutableRunner(command=["sleep", "5"])
    assert exec.close() is True
    assert exec.close() is False


def test_validate_executable_path_does_not_eexist() -> None:
    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable="/path/to/nowhere")


def test_validate_executable_path_not_executable() -> None:
    with NamedTemporaryFile(suffix=".exe", mode="w", delete=True) as tmpfile:
        with pytest.raises(ValueError, match="is not executable"):
            ExecutableRunner.validate_executable_path(executable=tmpfile.name)


def test_validate_executable_path() -> None:
    exec = "yes"
    exec_full_str = f"/usr/bin/{exec}"
    exec_full_path = Path(exec_full_str)
    assert exec_full_path.is_absolute()

    # find it on the PATH
    assert exec_full_path == ExecutableRunner.validate_executable_path(executable=exec)
    # find it given an absolute path as a string
    assert exec_full_path == ExecutableRunner.validate_executable_path(executable=exec_full_str)
    # find it given an absolute path as a Path
    assert exec_full_path == ExecutableRunner.validate_executable_path(executable=exec_full_path)

    # do not find it on the PATH if given as a Path
    with pytest.raises(ValueError, match="Executable does not exist"):
        ExecutableRunner.validate_executable_path(executable=Path(exec))


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
