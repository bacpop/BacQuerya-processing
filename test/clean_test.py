"""Remove test outputs"""

import os
import sys
import shutil

def deleteDir(dirname):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)

outputDirs = [
    "test_files",
    "test_features",
    "test_isolates",
    "test_attributes"
]

for outDir in outputDirs:
    deleteDir(outDir)
