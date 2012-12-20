#!/usr/bin/env python
import subprocess, sys

shellCmd = "tests/regression"
sys.exit(subprocess.call(shellCmd, shell=True))
