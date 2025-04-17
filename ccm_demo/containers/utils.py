import os
import subprocess


class ContainerRunError(Exception):
    pass


class ContainerRunner:
    def __init__(self, container, command_dict):
        """
        Run the command inside the container and return job status.
        :param command_dict:
        This needs to look like the example below:
            {"bind_mounts":[{"local":"local_dir", "target":<target_dir>},
                            etc],
             {"singularity_args": ["--cleanenv", "--nv" "etc"],
             {"run_command": ["/bin/bash", "-c"]},
            }
        :return: runs the specified command inside the container and then returns subprocess returncode
        """
        if not os.path.exists(container):
            raise FileNotFoundError(f"Container {container} not found")

        self.container = container
        self.command_dict = command_dict

    def run(self):

        binds = []
        for mount in self.command_dict["bind_mounts"]:
            if not os.path.exists(mount["local"]):
                raise NotADirectoryError(f"Bind mount {mount} not found")
            binds.append("-B {}:{}".format(os.path.abspath(mount["local"]), mount["target"]))

        singularity_args = " ".join(self.command_dict["singularity_args"])
        singularity_args = " ".join([singularity_args, binds])
        run_command = " ".join(self.command_dict["run_command"])
        command = "singularity exec {} {} {}".format(singularity_args, os.path.abspath(self.container), run_command)

        run = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.run = run

        return self

    def check_status(self):
        poll = self.run.poll()
        if poll is None:
            print("Container {} did not return a polling response. Meaning the process is still running".format(
                self.container))
        else:
            if self.run.returncode != 0:
                out, err = self.run.communicate()
                raise ContainerRunError("The process did not finish successfully. {}".format(err))
                self.out = out
                self.err = err
            else:
                out, err = self.run.communicate()
                print("Container run finished successfully")
                self.out = out
                self.err = err

        return self

    def return_results(self):
        to_ret = []
        for results in [self.out, self.err]:
            if results is not None:
                with open(results) as f:
                    f.write(results)
            else:
                to_ret.append(results)

        if len(to_ret) > 0:
            return to_ret
