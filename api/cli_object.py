import subprocess

class CLIObject:
    def __init__(self, apec):
        """
        CLIObject to manage CLI interaction for APEC_M.

        Parameters:
            apec: The APEC_M instance to bind CLI operations to.
        """
        self.apec = apec
        self.process = None
        self.num_processors = 1

    def start_process(self, command, num_processors=1):
        """
        Start a process and set it as an attribute. If num_processors > 1, use mpirun.

        Args:
            command (list or str): The command to execute as a list of strings or a single string.
            num_processors (int): The number of processors to use. Defaults to 1.
        """
        self.num_processors = num_processors
        command = command.split(" ") if isinstance(command, str) else command
        if self.num_processors > 1:
            command = ["mpirun", "-np", str(self.num_processors)] + command

        self.process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        print(f"Process started: {' '.join(command)}")
        while True:
            line = self.process.stdout.readline()
            if not line:
                break
            print(line.strip())
            if line.strip().endswith(">>>"):
                break

    def exec(self, command):
        """
        Send a command to the running process.

        Args:
            command (str): The command to send to the process.
        """
        if self.process and self.process.stdin:
            self.process.stdin.write(command + "\n")
            self.process.stdin.flush()
            print(f"Command sent: {command}")
        else:
            print("No process is running.")
        while True:
            line = self.process.stdout.readline()
            if not line:
                break
            print(line.strip())
            if line.strip().endswith(">>>"):
                break

    def exit(self):
        """
        Stop the running process.
        """
        if self.process:
            self.exec("-exit")
            self.process.wait()
            print("Process stopped.")
            self.process = None
        else:
            print("No process is running.")

    def start(self, num_processors=1):
        """
        Run an interactive session with the process.

        Args:
            num_processors (int): The number of processors to use. Defaults to 1.
        """
        try:
            self.start_process(["apec_m", "-i"], num_processors=num_processors)
            print("Interactive session started with `apec_m -i`.")
            print("Enter commands (e.g., -model arg, -solu, -solve arg, -finish, -exit).")
            print("Type 'exit' to terminate the session.")
            while True:
                command = input("> ")
                if command.strip().lower() == "exit":
                    self.exit()
                    break
                self.exec(command)
        except FileNotFoundError:
            print("Error: `apec_m` executable not found. Ensure it is in your PATH.")
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            if self.process and self.process.poll() is None:
                self.process.terminate()
