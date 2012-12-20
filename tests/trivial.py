#!/usr/bin/env python
import subprocess, sys

shellCmd = "tests/trivial"
sys.exit(subprocess.call(shellCmd, shell=True))
