#!/usr/bin/env python

## run_command_with_retries.py by Rohan Maddamsetti

import subprocess
import threading
import argparse


def run_command_with_retry(command_string, tempdir=None, max_retries=3, timeout=20):
    ## This code handles a bug in themisto build-- sometimes randomly hangs, have to delete temp files
    ## and restart and then it usually works.
    retries = 0
    while retries < max_retries:
        process = subprocess.Popen(command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        timer = threading.Timer(timeout, process.kill)  ## Kill process if it exceeds timeout

        try:
            timer.start()
            stdout, stderr = process.communicate()
        finally:
            timer.cancel()

        if process.returncode == 0:
            print("Command succeeded:", stdout.decode())
            return stdout.decode()
        else:
            print(f"*********COMMAND FAILED (attempt {retries + 1}):", stderr.decode())
            if tempdir is not None: ## remove temporary files from the failed run.
                print(f"removing {tempdir}")
                subprocess.run(f"rm -rf {tempdir}", shell=True)
                print(f"remaking {tempdir} before restarting")
                os.mkdir(tempdir)
            retries += 1
            time.sleep(0.1)  ## Small delay before retrying
    
    print("Command failed after maximum retries.")
    return


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a command with retries.",
        usage="python run_command_with_retries.py "
    )

    ## Required option for command string
    parser.add_argument(
        "cmd_string",
        type=str,
        required=True,
        help="command to run with retries"
    )

    parser.add_argument(
        "tempdir",
        type=str,
        default=None,
        help="path to temporary directory to delete if the command fails and needs to be retried."
    )

    return parser.parse_args()


def main():
    args = parse_args()
    run_command_with_retry(args.cmd_string, args.tempdir)
    return


if __name__ == "__main__":
    main()
    
