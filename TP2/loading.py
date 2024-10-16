import gzip
from os import listdir, path


def load_fasta(file_pointer):
    """ Loads a fasta formated file into a list of sequences.
    :param file_pointer: The stream of the fasta file to load.
    :return Array: An array of strings where each string is a sequence from the fasta
    """
    texts = []
    txt = []

    for line in file_pointer:
        if line[0] == '>':
            if len(txt) > 0:
                texts.append("".join(txt))
            txt = []
        else:
            txt.append(line.strip())

    if len(txt) > 0:
        texts.append("".join(txt))
    return texts



def load_directory(directory):
    """ Loads all the fasta files from a data directory into a dictionary.
    Each subdirectory in data is considered as a different sample.
    Fatsta files (even gzipped) are loaded.
    :param str directory: Path to the data directory to load.
    :return dict: A dict containing pairs (sample, [sequence list]).
    """
    print("Loading data from directory", directory)
    files = {}
    for name in listdir(directory):
        subpath = path.join(directory, name)
        # Look for sample directories
        if path.isdir(subpath):
            # Creates one list of sequence per sample
            files[name] = []
            for filename in listdir(subpath):
                # Load raw fasta files
                if filename.endswith(".fa") or filename.endswith(".fasta"):
                    with open(path.join(subpath, filename)) as fp:
                        files[name] += load_fasta(fp)
                # Load gzipped fasta files
                elif filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
                    with gzip.open(path.join(subpath, filename), 'rt') as fp:
                        files[name] += load_fasta(fp)
                        print("Loaded", filename, len(files[name]))
    
    return files


if __name__ == "__main__":
    files = load_directory("data")
    print(len(files))
