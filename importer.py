''' used by the parameters class to import the static data
'''
import json
import sys


class Importer(object):

    def __init__(self, filepath):
        self.filepath = "./" + filepath + "/"

    @staticmethod
    def get_lines_from_file(filename):
        fin = open(filename)
        lines = fin.readlines()
        fin.close()
        return lines

    @staticmethod
    def parse_json(filename):
        with open(filename, 'r') as fp:
            payload = json.load(fp)
        return payload
