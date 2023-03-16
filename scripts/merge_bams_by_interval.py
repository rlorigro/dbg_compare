from module.Authenticator import GoogleToken
from multiprocessing import Pool

import subprocess
import argparse
import sys
import os
import re


def merge_bams_by_interval(bam_paths, contig, start, stop, output_directory):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_filename = region_string.replace(":","_") + ".fasta"
    output_path = os.path.join(output_directory,output_filename)

    # Not attempting to share a token across threads, so each job gets a new instance
    token = GoogleToken()

    for path in bam_paths:
        get_remote_region_as_fasta(
            bam_path=path,
            contig=contig,
            start=start,
            stop=stop,
            output_path=output_path,
            token=token
        )


# Requires samtools installed!
def get_remote_region_as_fasta(bam_path, contig, start, stop, output_path, token):
    region_string = "%s:%d-%d" % (contig, start, stop)

    print(region_string)

    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()

    samtools_view_args = ["samtools", "view", "-b", "-h", "-F", "4", bam_path, region_string]
    samtools_fasta_args = ["samtools", "fasta", "-"]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_view_args)+'\n')

        p1 = subprocess.Popen(samtools_view_args, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_fasta_args, stdin=p1.stdout, stdout=file)
        p2.communicate()

    return output_path


def main(bam_paths, bed_path, output_directory, n_cores):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Contains all the necessary args to fetch a region from the BAMs
    args = list()

    with open(bed_path, 'r') as file:
        for l,line in enumerate(file):
            contig,start,stop = line.strip().split()
            start = int(start)
            stop = int(stop)

            args.append([bam_paths, contig, start, stop, output_directory])

    with Pool(processes=n_cores) as pool:
        pool.starmap(merge_bams_by_interval, args)


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bams",
        required=True,
        type=parse_comma_separated_string,
        help="Input BAMs to be chunked (comma separated list)"
    )

    parser.add_argument(
        "--bed",
        required=True,
        type=str,
        help="BED file containing intervals which are fetched from all BAMs and merged"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="Number of cores to use (for parallelizing across intervals)"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    print(args.bams)
    exit()

    main(bam_paths=args.bams, bed_path=args.bed, output_directory=args.o, n_cores=args.c)
