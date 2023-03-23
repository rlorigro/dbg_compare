from module.Authenticator import GoogleToken
from multiprocessing import Pool

import subprocess
import argparse
import hashlib
import tarfile
import shutil
import sys
import os
import re
import io


def merge_bams_by_interval(bam_paths, contig, start, stop, output_directory, token):

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
def get_remote_region_as_bam(bam_path, contig, start, stop, output_directory, token, index=True, region_label=False):

    region_string = "%s:%d-%d" % (contig, start, stop)

    if region_label:
        output_filename = os.path.basename(bam_path).split('.')[0] + "_" + region_string.replace(":","_") + ".bam"
    else:
        output_filename = os.path.basename(bam_path).split('.')[0] + ".bam"

    output_path = os.path.join(output_directory,output_filename)

    samtools_view_args = ["samtools", "view", "-b", "-h", "-F", "4", "-o", output_path, bam_path, region_string]

    sys.stderr.write(" ".join(samtools_view_args)+'\n')

    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()
    try:
        p1 = subprocess.run(samtools_view_args, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.write("ERROR: " + " ".join([bam_path, contig, str(start), str(stop)]) + '\n')
        sys.stderr.flush()
        return None

    if index:
        samtools_index_args = ["samtools", "index", output_path]
        sys.stderr.write(" ".join(samtools_index_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_index_args, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.write("ERROR: " + " ".join([bam_path, contig, str(start), str(stop)]) + '\n')
            sys.stderr.flush()
            return None

    return output_path


def get_region_coverage(bam_path, contig, start, stop, output_directory, token):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_path = os.path.join(output_directory, os.path.basename(bam_path))
    output_path = output_path.replace(".bam", "_coverage.tsv")

    print(output_path)

    samtools_args = ["samtools", "coverage", "-r", region_string, bam_path]

    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_args, stdout=file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.write("ERROR: " + " ".join([bam_path, contig, str(start), str(stop)]) + '\n')
            sys.stderr.flush()
            return None

    return output_path


def get_reads_from_bam(bam_path, output_directory, token):
    output_path = os.path.join(output_directory, os.path.basename(bam_path))
    output_path = output_path.replace(".bam", ".fasta")

    print(output_path)

    samtools_args = ["samtools", "fasta", bam_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.write("ERROR: " + bam_path + '\n')
            sys.stderr.flush()
            return None

    return output_path


def merge_coverages(output_directory):
    output_filename = "coverage.tsv"
    output_path = os.path.join(output_directory, output_filename)

    with open(output_path, 'w') as out_file:
        f = 0
        for filename in os.listdir(output_directory):
            if filename.endswith(".tsv") and filename != output_filename:
                path = os.path.join(output_directory, filename)
                sample = filename.split("_")[0]

                with open(path, 'r') as file:
                    for l,line in enumerate(file):

                        # Only write the header once
                        if f == 0 and l == 0:
                            out_file.write("sample\t" + line[1:])

                        # Always write the values
                        if l == 1:
                            out_file.write(sample + '\t' + line)

                os.remove(path)
                f += 1


def process_region(bam_paths, contig, start, stop, output_directory, token):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_subdirectory = region_string.replace(":","_")
    output_directory = os.path.join(output_directory, output_subdirectory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        sys.stderr.write("WARNING: duplicate region in BED file, skipping %s:%d-%d\n" % (contig, start, stop))
        return

    for path in bam_paths:
        local_bam_path = get_remote_region_as_bam(
            bam_path=path,
            contig=contig,
            start=start,
            stop=stop,
            output_directory=output_directory,
            token=token)

        if local_bam_path is None:
            return

        coverage_csv_path = get_region_coverage(
            bam_path=local_bam_path,
            contig=contig,
            start=start,
            stop=stop,
            output_directory=output_directory,
            token=token)

        if coverage_csv_path is None:
            return

        fasta_path = get_reads_from_bam(
            bam_path=local_bam_path,
            output_directory=output_directory,
            token=token)

        if fasta_path is None:
            return

        os.remove(local_bam_path)
        os.remove(local_bam_path + ".bai")

    merge_coverages(output_directory)

    with tarfile.open(output_directory + ".tar.gz", "w:gz") as tar:
        tar.add(output_directory, arcname=os.path.basename(output_directory))

    shutil.rmtree(output_directory)

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

    token = GoogleToken()

    with open(bed_path, 'r') as file:
        for l,line in enumerate(file):
            contig,start,stop = line.strip().split()
            start = int(start)
            stop = int(stop)

            args.append([bam_paths, contig, start, stop, output_directory, token])

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
