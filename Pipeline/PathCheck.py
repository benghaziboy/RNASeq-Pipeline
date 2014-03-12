import sys
import os


def check_folder(folder):
    while True:
        if os.path.isdir(os.path.dirname(folder)):
            return
        else:
            try:
                os.mkdir(os.path.dirname(folder))
            except:
                check_folder(os.path.dirname(os.path.abspath(folder)))
    return folder
