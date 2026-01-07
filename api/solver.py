from __future__ import annotations

import os
import select
import subprocess
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from apec_m import APEC_M


class Solver:
    """
    Thin Python helper for driving the C++ solver.
    """

    prompt: str = "apec_m>"

    def __init__(self, apec: "APEC_M", num_processors: int = 1) -> None:
        """
        Initialize a solver bound to an ``APEC_M`` instance.

        Args:
            apec: The ``APEC_M`` instance this CLI helper is associated with.
            num_processors: Default MPI rank count to use when launching the solver.
        """
        self.apec: "APEC_M" = apec
        self.process: Optional[subprocess.Popen[str]] = None
        self.num_processors: int = num_processors

    def start_process(self, num_processors: Optional[int] = None) -> None:
        """
        Launch the solver process in interactive mode.

        Args:
            num_processors: Number of MPI ranks; when >1 the command is prefixed with
                ``mpirun -np``. Defaults to the value provided at construction.
        """
        if num_processors is not None:
            self.num_processors = num_processors
        cmd_list = ["apec_m", "-i"]
        ranks = max(1, self.num_processors)
        if ranks > 1:
            cmd_list = ["mpirun", "-np", str(ranks)] + cmd_list

        self.process = subprocess.Popen(
            cmd_list,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        print(f"Process started: {' '.join(cmd_list)}")
        self._consume_output(wait_for_prompt=True)

    def exec(
        self, command: str, *, wait_for_prompt: bool = True, echo: bool = True
    ) -> str:
        """
        Send a CLI command to the running solver.

        Args:
            command: Raw CLI text such as ``-read_mesh mesh.vtu``.
            wait_for_prompt: When True, read output until the solver prompt is seen or
                stdout goes idle.
            echo: Print the solver output to stdout.

        Returns:
            Captured solver output for this command.
        """
        process = self._ensure_process()
        if not process.stdin:
            raise RuntimeError("Process stdin is not available.")

        process.stdin.write(command + "\n")
        process.stdin.flush()
        if echo:
            print(f"Command sent: {command}")

        return self._consume_output(wait_for_prompt=wait_for_prompt, echo=echo)

    def solu(
        self,
        analysis_type: str,
        *,
        mesh_file: Optional[str] = None,
        model_filename: str = "model_data.bin",
    ) -> str:
        """
        Write model data to binary, load it into the solver, and enter analysis mode.

        Args:
            analysis_type: Analysis mode requested by the solver
                (e.g., ``static``, ``transient``, ``apec_equilibrium``).
            mesh_file: Optional mesh filepath to read via ``-read_mesh`` before loading
                model data.
            model_filename: Relative filename for the generated model-data binary.
        """
        if mesh_file:
            self.exec(f"-read_mesh {mesh_file}")

        model_path = self._write_model_data(model_filename)
        self.exec(f"-load_data {model_path}")
        return self.exec(f"-mode analysis {analysis_type}")

    def solve(
        self,
        load_step_number: int,
        *,
        grid_filename: Optional[str] = None,
        time_filename: Optional[str] = None,
    ) -> str:
        """
        Write load-step binaries and solve that load step.

        Args:
            load_step_number: 1-based load step index.
            grid_filename: Relative grid-data filename to write and load.
            time_filename: Relative time-step filename to write and load.

        Returns:
            Solver output from the ``-solve`` command.
        """
        if load_step_number < 1:
            raise ValueError("load_step_number must be >= 1.")

        time_path, grid_path = self._write_load_step_binaries(
            load_step_number, grid_filename, time_filename
        )
        self.exec(f"-load_data {time_path} {grid_path}")
        return self.exec("-solve")

    def exit(self) -> str:
        """
        Gracefully stop the solver process (CLI: ``-exit``).
        """
        if not self.process:
            return ""

        try:
            output = self.exec("-exit", wait_for_prompt=False)
        finally:
            self.process.wait()
            self.process = None
        print("Process stopped.")
        return output

    def is_running(self) -> bool:
        """
        Check whether the solver process is alive.
        """
        return self.process is not None and self.process.poll() is None

    def start(self, num_processors: Optional[int] = None) -> None:
        """
        Start an interactive CLI session (``apec_m -i``).
        """
        try:
            self.start_process(num_processors=num_processors)
            print(
                "Enter commands such as -read_mesh <mesh>, -active_mesh <name>, "
                "-load_data <files>, -solve, -solve_transient, -finish, "
                "-mode <...>, -print_model_data, -clear, -help, -exit."
            )

            while self.is_running():
                command = input(f"{self.prompt} ").strip()
                if command.lower() in {"exit", "-exit", "quit"}:
                    self.exit()
                    break
                if command:
                    self.exec(command)
        except FileNotFoundError:
            print("Error: `apec_m` executable not found. Ensure it is in your PATH.")
        except KeyboardInterrupt:
            print("\nInterrupted. Exiting solver.")
            self.exit()
        except Exception as exc:  # noqa: BLE001
            print(f"An error occurred: {exc}")
            self.exit()

    def _ensure_process(self) -> subprocess.Popen[str]:
        if not self.process:
            raise RuntimeError("No process is running. Call start_process first.")
        return self.process

    def _consume_output(
        self,
        *,
        wait_for_prompt: bool = True,
        echo: bool = True,
        idle_timeout: float = 0.2,
    ) -> str:
        """
        Read solver stdout until the prompt reappears or the process exits.
        """
        process = self._ensure_process()
        if not process.stdout:
            return ""

        output: list[str] = []
        buffer = ""

        while True:
            ready, _, _ = select.select([process.stdout], [], [], idle_timeout)
            if not ready:
                if not wait_for_prompt and process.poll() is not None:
                    break
                continue

            chunk = process.stdout.read(1)
            if chunk == "":
                break

            buffer += chunk
            output.append(chunk)
            if wait_for_prompt and self.prompt in buffer:
                break

        output_text = "".join(output)
        if echo and output_text:
            print(output_text.strip("\n"))
        return output_text

    def _write_model_data(self, filename: str) -> str:
        """
        Serialize model-level data to binary and return its path.
        """
        folder = self.apec.folderpath
        if not folder:
            raise ValueError("Folder path is not set; cannot write model data.")
        model_data_folder = os.path.join(folder, "_model_data_")
        output_path = os.path.join(model_data_folder, filename)
        self.apec.write_model_data(filename=filename)
        return output_path

    def _write_load_step_binaries(
        self,
        load_step_number: int,
        grid_filename: Optional[str],
        time_filename: Optional[str],
    ) -> tuple[str, str]:
        """
        Serialize time and grid data for a load step and return their filepaths.
        """
        folder = self.apec.folderpath
        jobname = getattr(self.apec, "jobname", None)
        if not folder:
            raise ValueError("Folder path is not set; cannot write load step binaries.")
        if not jobname and (grid_filename is None or time_filename is None):
            raise ValueError("Job name is not set; cannot infer default binary paths.")

        time_folder = os.path.join(folder, "_time_")
        time_path = (
            os.path.join(time_folder, time_filename)
            if time_filename
            else os.path.join(time_folder, f"{jobname}_timestep_{load_step_number}.bin")
        )
        grid_path = (
            os.path.join(folder, grid_filename)
            if grid_filename
            else os.path.join(
                folder, jobname, "input", f"{jobname}_grid_{load_step_number}.bin"
            )
        )

        self.apec.write_time_step(load_step_number, filename=time_filename)
        self.apec.write_loads_bcs(load_step_number, filename=grid_filename)
        return time_path, grid_path
