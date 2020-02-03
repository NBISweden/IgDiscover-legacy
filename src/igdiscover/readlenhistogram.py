from collections import Counter
import dnaio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def length_histogram(path):
    """Return a list of lengths """
    lengths = []
    with dnaio.open(path) as reader:
        for record in reader:
            lengths.append(len(record.sequence))
    return lengths


def plot_histogram(lengths, path, title, min_x=0, bins=50):
    """
    Plot histogram of lengths to path
    """
    lengths = np.array(lengths)
    histomax = int(max(lengths, default=100))
    fig = plt.figure(figsize=(20/2.54, 10/2.54))
    ax = fig.gca()
    ax.set_xlabel('Read length')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    _, borders, _ = ax.hist(lengths, bins=bins, range=(min_x, histomax))
    fig.set_tight_layout(True)
    fig.savefig(path)


def read_length_histogram(path, tsv_path, plot_path, bins, left, title):
    lengths = length_histogram(path)
    freqs = Counter(lengths)
    with open(tsv_path, "w") as outfile:
        print("## File:", path, file=outfile)
        print("length", "frequency", sep='\t', file=outfile)
        for length in range(0, max(freqs, default=100) + 1):
            freq = freqs[length]
            if freq > 0:
                print(length, freq, sep='\t', file=outfile)

    if plot_path:
        plot_histogram(lengths, plot_path, title, min_x=left, bins=bins)
