from module.Authenticator import GoogleToken
from multiprocessing import Pool

import subprocess
import argparse
import sys
import os
import re


def merge_bams_by_interval(bam_paths, contig, start, stop, output_directory):
    # Not attempting to share a token across threads, so each job gets a new instance
    token = GoogleToken()

    for path in bam_paths:
        result = get_remote_region_as_fasta(
            bam_path=path,
            contig=contig,
            start=start,
            stop=stop,
            output_directory=output_directory,
            token=token
        )

        if not result[0]:
            return result

    return True, None


# Requires samtools installed!
def get_remote_region_as_fasta(bam_path, contig, start, stop, output_directory, token):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_filename = region_string.replace(":","_") + ".fasta"
    output_path = os.path.join(output_directory,output_filename)

    print(region_string)

    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()

    samtools_view_args = ["samtools", "view", "-b", "-h", "-F", "4", bam_path, region_string]
    samtools_fasta_args = ["samtools", "fasta", "-"]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_view_args)+'\n')

        p1 = subprocess.Popen(samtools_view_args, stdout=subprocess.PIPE, cwd=output_directory)
        p2 = subprocess.Popen(samtools_fasta_args, stdin=p1.stdout, stdout=file, cwd=output_directory)
        p2.communicate()

    success = (p2.returncode == 0)

    if not success:
        sys.stderr.flush()

    return success, bam_path


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
        successes = pool.starmap(merge_bams_by_interval, args)

    all_succeeded = True
    for s,success in enumerate(successes):
        if not success[0]:
            sys.stderr.write("ERROR: job for region %s:%d-%d failed\n" % (args[s][1], args[s][2], args[s][3]))
            sys.stderr.write("\t%s\n" % (success[1]))
            all_succeeded = False

    if not all_succeeded:
        exit(1)


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

    main(bam_paths=args.bams, bed_path=args.bed, output_directory=args.o, n_cores=args.c)
