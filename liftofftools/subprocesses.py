import subprocess
import sys

class SubprocessCaller():
    def __init__(self, software, subcommand, positionals, command_line_args={}):
        self.software = software
        self.subcommand = subcommand
        self.command_line_args = command_line_args
        self.positionals = positionals


    def run_command(self):
            command = self.build_command()
            try:
                subprocess.check_call(command)
            except subprocess.CalledProcessError as e:
                print(e)
                sys.exit()


    def build_command(self):
        command_line_list = [self.software]
        if self.subcommand != '':
            command_line_list += [self.subcommand]
        command_line_list += self.positionals
        for param, value in self.command_line_args.items():
            command_line_list.append(param)
            command_line_list.append(value)
        return command_line_list