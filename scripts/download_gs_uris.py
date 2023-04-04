from module.GsUri import download_gs_uri
from multiprocessing import Pool

import argparse
import sys
import re
import os


def main(input_path, n_threads, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    uri_paths = list()
    with open(input_path, 'r') as file:
        for l,line in enumerate(file):
            uri_paths.append(line.strip())

    # Multithread the downloading of files
    args = [[p,output_directory] for p in uri_paths]
    with Pool(n_threads) as pool:
        download_results = pool.starmap(download_gs_uri, args)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input file containing URIs (one per line), to be downloaded"
    )

    parser.add_argument(
        "-t",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory, where files will be downloaded using their filenames on the cloud"
    )

    args = parser.parse_args()

    main(input_path=args.i, n_threads=args.t, output_directory=args.o)
