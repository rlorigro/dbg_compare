#!/usr/bin/env python
from pysam import AlignmentFile
from module.Authenticator import GoogleToken
from collections import Counter

import argparse
import random
import sys
import os


def get_remote_header(bam_path, token):
    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()
    header = AlignmentFile(bam_path,'r').header

    return header


def main(bam_path, chunk_size, n_samples, forbidden, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    token = GoogleToken()
    header = get_remote_header(bam_path=bam_path, token=token)

    contig_lengths = list()
    total_length = 0

    for r in header.references:
        if r in forbidden:
            continue
        l = int(header.get_reference_length(r))
        contig_lengths.append([r,l])
        total_length += l

    regions = list()

    contig_weights = [float(x[1])/float(total_length) for x in contig_lengths]
    contig_choices = random.choices(contig_lengths, weights=contig_weights, k=n_samples)

    print(Counter([x[0] for x in contig_choices]))

    for contig,length in contig_choices:
        start = random.randint(1,max(1,length-chunk_size+1))
        stop = min(length,start + chunk_size - 2)

        regions.append([contig,start,stop])

    output_path = os.path.join(output_directory, "intervals.bed")
    with open(output_path, 'w') as file:
        for r in regions:
            file.write('\t'.join(list(map(str,r))))
            file.write('\n')


def parse_comma_separated_string(s):
    return s.strip().split(',')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input BAMs to be chunked (comma separated list)"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="Chunk size"
    )

    parser.add_argument(
        "-n",
        required=True,
        type=int,
        help="Number of samples of `chunk_size` to extract randomly"
    )

    parser.add_argument(
        "-f",
        required=True,
        type=parse_comma_separated_string,
        help="Forbidden contigs, as a comma-separated list"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(bam_path=args.i, output_directory=args.o, chunk_size=args.c, forbidden=args.f, n_samples=args.n)
