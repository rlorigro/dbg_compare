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


# Requires samtools installed!
def get_remote_region_as_bam(bam_path, output_path, contig, start, stop, token, index=True):
    region_string = "%s:%d-%d" % (contig, start, stop)

    samtools_view_args = ["samtools", "view", "-b", "-h", "-F", "4", "-o", output_path, bam_path, region_string]

    sys.stderr.write(" ".join(samtools_view_args)+'\n')

    # There is a small chance that this will fail if the token expires between updating and downloading...
    token.update_environment()
    try:
        p1 = subprocess.run(samtools_view_args, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return False
    except Exception as e:
        sys.stderr.write(str(e))
        return False

    if index:
        samtools_index_args = ["samtools", "index", output_path]
        sys.stderr.write(" ".join(samtools_index_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_index_args, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return True


def get_region_coverage(bam_path, output_path, contig, start, stop, token):
    region_string = "%s:%d-%d" % (contig, start, stop)
    samtools_args = ["samtools", "coverage", "-r", region_string, bam_path]

    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_args, stdout=file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return True


def get_reads_from_bam(bam_path, output_path, token):
    samtools_args = ["samtools", "fasta", bam_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return True


def merge_coverages(output_directory, expected_sample_count=None):
    output_filename = "coverage.tsv"
    output_path = os.path.join(output_directory, output_filename)

    with open(output_path, 'w') as out_file:
        f = 0
        for filename in os.listdir(output_directory):
            if filename.endswith(".tsv") and filename != output_filename:
                path = os.path.join(output_directory, filename)
                sample = filename.replace("_coverage.tsv","")

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

    if f != expected_sample_count:
        sys.stderr.write("ERROR: number of coverage TSV objects != number of samples: %d != %d\n" % (f, expected_sample_count))
        return False
    else:
        return True


def process_region(bam_paths, contig, start, stop, output_directory, token):
    region_string = "%s:%d-%d" % (contig, start, stop)
    output_subdirectory = region_string.replace(":","_")
    output_directory = os.path.join(output_directory, output_subdirectory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        sys.stderr.write("WARNING: duplicate region in BED file, skipping %s:%d-%d\n" % (contig, start, stop))
        return

    paths_to_validate = list()
    all_success = True

    for bam_path in bam_paths:
        sample_name = os.path.basename(bam_path).split('.')[0]

        local_bam_filename = sample_name + "_" + region_string.replace(":","_") + ".bam"
        local_bam_path = os.path.join(output_directory,local_bam_filename)

        fasta_path = os.path.join(output_directory, sample_name + ".fasta")
        coverage_path = os.path.join(output_directory, sample_name + "_coverage.tsv")

        paths_to_validate.append(fasta_path)
        paths_to_validate.append(coverage_path)

        success = get_remote_region_as_bam(
            bam_path=bam_path,
            output_path=local_bam_path,
            contig=contig,
            start=start,
            stop=stop,
            token=token)

        if not success:
            sys.stderr.write("ERROR: failed to fetch BAM: %s %s\n" % (bam_path, region_string))
            all_success = False

        get_region_coverage(
            bam_path=local_bam_path,
            output_path=coverage_path,
            contig=contig,
            start=start,
            stop=stop,
            token=token)

        if not success:
            sys.stderr.write("ERROR: failed to get region coverage: %s %s\n" % (coverage_path, region_string))
            all_success = False

        get_reads_from_bam(
            bam_path=local_bam_path,
            output_path=fasta_path,
            token=token)

        if not success:
            sys.stderr.write("ERROR: failed to get reads from BAM: %s %s\n" % (local_bam_path, region_string))
            all_success = False

        if os.path.exists(local_bam_path):
            os.remove(local_bam_path)
        else:
            sys.stderr.write("ERROR: local BAM does not exist: %s %s\n" % (local_bam_path, region_string))
            all_success = False

        # If possible, would want to keep the bai around to reduce re-downloading for each thread/region,
        # assuming there will be no duplicate BAM filenames. Logistically more annoying though.
        if os.path.exists(local_bam_path + ".bai"):
            os.remove(local_bam_path + ".bai")
        else:
            sys.stderr.write("ERROR: local BAM index does not exist: %s %s\n" % (local_bam_path, region_string))
            all_success = False

    # Do a sanity check to see that paths exist
    for path in paths_to_validate:
        if not os.path.exists(path):
            sys.stderr.write("ERROR: expected file path not found: %s\n" % path)
            all_success = False

    merge_coverages(output_directory=output_directory, expected_sample_count=len(bam_paths))

    if all_success:
        with tarfile.open(output_directory + ".tar.gz", "w:gz") as tar:
            tar.add(output_directory, arcname=os.path.basename(output_directory))
    else:
        sys.stderr.write("ERROR: region %s skipped because one or more samples resulted in error\n" % (region_string))

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
        sys.stderr.write("ERROR: failed to fetch: %s\n" % (" ".join([bam_path, contig, start, stop])))
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
        pool.close()
        pool.join()

    sys.stderr.write("Files prepared:\n")
    for filename in os.listdir(output_directory):
        sys.stderr.write(filename)
        sys.stderr.write('\n')


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
