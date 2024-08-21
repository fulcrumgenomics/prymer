"""
# Base classes and methods for wrapping subprocess

This module contains a base class to facilitate wrapping subprocess and run command line tools from
Python. Methods include functions to validate executable paths as well as initiate
and interact with subprocesses. This base class implements the context manager protocol.

"""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from types import TracebackType
from typing import Optional
from typing import Self


class ExecutableRunner:
    """
    Base class for interaction with subprocess for all command-line tools. The base class supports
    use of the context management protocol and performs basic validation of executable paths.

    The constructor makes the assumption that the first path element of the command will be
    the name of the executable being invoked. The constructor initializes a subprocess with
    file handles for stdin, stdout, and stderr, each of which is opened in text mode.

    Subclasses of [`ExecutableRunner`][prymer.util.executable_runner.ExecutableRunner]
    provide additional type checking of inputs and orchestrate parsing output data from specific
    command-line tools.
    """

    __slots__ = ("_command", "_subprocess", "_name")
    _command: list[str]
    _subprocess: subprocess.Popen[str]
    _name: str

    def __init__(
        self,
        command: list[str],
        stdin: int = subprocess.PIPE,
        stdout: int = subprocess.PIPE,
        stderr: int = subprocess.PIPE,
    ) -> None:
        if len(command) == 0:
            raise ValueError(f"Invocation must not be empty, received {command}")
        self._command = command
        self._name = command[0]
        self._subprocess = subprocess.Popen(
            command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            text=True,
            bufsize=0,  # do not buffer stdin/stdout so that we can read/write immediately
        )

    def __enter__(self) -> Self:
        logging.debug(f"Initiating {self._name} with the following params: {self._command}")
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        """Gracefully terminates any running subprocesses."""
        self.close()

    @classmethod
    def validate_executable_path(cls, executable: str | Path) -> Path:
        """Validates user-provided path to an executable.

        If a string is provided, checks whether a Path representation exists. If not, uses
        shutil.which() to find the executable based on the name of the command-line tool.

         Args:
            executable: string or Path representation of executable

         Returns:
            Path: valid path to executable (if found)

         Raises:
            ValueError: if path to executable cannot be found
            ValueError: if executable is not executable
        """
        if isinstance(executable, str):
            executable = Path(executable)
            if not executable.exists() and executable.name == f"{executable}":
                retval = shutil.which(f"{executable}", mode=os.F_OK)  # check file existence
                if retval is not None:
                    executable = Path(retval)

        if not executable.exists():
            raise ValueError(f"Executable does not exist: {executable}")
        if not os.access(executable, os.X_OK):  # check file executability
            raise ValueError(f"`{executable}` is not executable: {executable}")

        return executable

    @property
    def is_alive(self) -> bool:
        """
        Check whether a shell subprocess is still alive.

        Returns:
            bool: True if process is alive, False if otherwise

        """
        return self._subprocess.poll() is None

    def close(self) -> bool:
        """
        Gracefully terminates the underlying subprocess if it is still
        running.

        Returns:
            True: if the subprocess was terminated successfully
            False: if the subprocess failed to terminate or was not already running
        """
        if self.is_alive:
            self._subprocess.terminate()
            self._subprocess.wait(timeout=10)
            if not self.is_alive:
                logging.debug("Subprocess terminated successfully.")
                return True
            else:
                logging.debug("Subprocess failed to terminate.")
                return False
        else:
            logging.debug("Subprocess is not running.")
            return False
