from module.Authenticator import GoogleToken
from multiprocessing import Pool

import subprocess
import argparse
import hashlib
import sys
import os
import re
import io


def merge_bams_by_interval(bam_paths, contig, start, stop, output_directory):
    # Not attempting to share a token across threads, so each job gets a new instance
    token = GoogleToken()

    for path in bam_paths:
        result = get_remote_region_as_bam(
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
def get_remote_region_as_bam(bam_path, contig, start, stop, output_directory, token, index=True):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_filename = os.path.basename(bam_path).split('.')[0] + "_" + region_string.replace(":","_") + ".bam"
    output_path = os.path.join(output_directory,output_filename)

    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()

    samtools_view_args = ["samtools", "view", "-b", "-h", "-F", "4", "-o", output_filename, bam_path, region_string]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_view_args)+'\n')

        p1 = subprocess.run(samtools_view_args, cwd=output_directory, check=True)

    success = (p1.returncode == 0)

    print("success: %s" % success)

    if not success:
        sys.stderr.write("ERROR: " + " ".join([bam_path, contig, start, stop]))
        sys.stderr.flush()

    if index:
        samtools_index_args = ["samtools", "index", output_filename]
        sys.stderr.write(" ".join(samtools_index_args)+'\n')
        p1 = subprocess.run(samtools_index_args, cwd=output_directory, check=True)

        success = (p1.returncode == 0)

        print("success: %s" % success)

        if not success:
            sys.stderr.write("ERROR: " + " ".join([bam_path, contig, start, stop]))
            sys.stderr.flush()

    return output_path


def get_region_coverage(bam_path, contig, start, stop, output_directory):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_path = os.path.join(output_directory, os.path.basename(bam_path))
    output_path = output_path.replace(".bam", "_coverage.csv")

    print(output_path)

    samtools_args = ["samtools", "coverage", "-r", region_string, bam_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        p1 = subprocess.run(samtools_args, stdout=file, cwd=output_directory, check=True)

    success = (p1.returncode == 0)

    if not success:
        sys.stderr.write("ERROR: " + " ".join([bam_path, contig, start, stop]))
        sys.stderr.flush()

    return output_path


def get_reads_from_bam(bam_path, output_directory):
    output_path = os.path.join(output_directory, os.path.basename(bam_path))
    output_path = output_path.replace(".bam", ".fasta")

    print(output_path)

    samtools_args = ["samtools", "fasta", bam_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        p1 = subprocess.run(samtools_args, stdout=file, cwd=output_directory, check=True)

    success = (p1.returncode == 0)

    if not success:
        sys.stderr.write("ERROR: " + " ".join([bam_path]))
        sys.stderr.flush()

    return output_path


def process_region(bam_paths, contig, start, stop, output_directory):
    # Not attempting to share a token across threads, so each job gets a new instance
    token = GoogleToken()

    region_string = "%s:%d-%d" % (contig, start, stop)
    output_subdirectory = region_string.replace(":","_")
    output_directory = os.path.join(output_directory, output_subdirectory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for path in bam_paths:
        local_bam_path = get_remote_region_as_bam(
            bam_path=path,
            contig=contig,
            start=start,
            stop=stop,
            output_directory=output_directory,
            token=token)

        coverage_csv_path = get_region_coverage(
            bam_path=local_bam_path,
            contig=contig,
            start=start,
            stop=stop,
            output_directory=output_directory)

        fasta_path = get_reads_from_bam(
            bam_path=local_bam_path,
            output_directory=output_directory)
        
        os.remove(local_bam_path)
        os.remove(local_bam_path + ".bai")
        os.remove(os.path.join(output_directory,os.path.basename(path) + ".bai"))

    return


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
        sys.stderr.write("ERROR: " + " ".join([bam_path, contig, start, stop]))
        sys.stderr.flush()

    return


def main(bam_paths, bed_path, output_directory, n_cores):
    output_directory = os.path.abspath(output_directory)

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
        results = pool.starmap(process_region, args)

    # all_succeeded = True
    # for r,result in enumerate(results):
    #     result.get()


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
