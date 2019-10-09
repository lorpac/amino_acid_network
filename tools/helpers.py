import os
import json
from shutil import copyfile

#TODO change configurationfile parameter with configuration file path and remove __file__ line. Check for scripts and notebboks.
# in the notebook I have one level less!!

def get_configuration_parameters(configuration_file_path):
    with open(configuration_file_path) as f:
        parameters = json.load(f)
    return parameters

def save_configuration_parameters(configuration_file_path, folder_path):
    configuration_file = os.path.split(configuration_file_path)[1]
    copyfile(configuration_file_path, os.path.join(folder_path, configuration_file))