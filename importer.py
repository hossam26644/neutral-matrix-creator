''' used by the parameters class to import the static data
'''
import json
import sys


class Importer(object):

    def __init__(self, filepath):
        self.filepath = ""

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

    @staticmethod
    def export_json(filename, data):
        with open(filename, 'w') as f:
            json.dump(data, f)

    @staticmethod
    def export_dict_to_tsv(filename, data):
        nuclutieds = ['A', 'C', 'G', 'T']
        text = ''
        for key, values in data.items():
            text += key + '\t'
            for nuclutied in nuclutieds:
                text += str(values[nuclutied]) + '\t'
            text += '\n'

        with open(filename, 'w') as f:
            f.write(text)
