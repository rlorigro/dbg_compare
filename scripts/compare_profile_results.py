from module.Authenticator import GoogleToken
from google.cloud import storage
import tarfile
import pandas
import argparse
import re
import os


def decode_gs_uri(uri):
    tokens = uri.split('/')
    bucket, file_path = tokens[2], '/'.join(tokens[3:])
    return bucket, file_path


def download_gs_uri(uri, output_directory):
    storage_client = storage.Client()
    bucket, file_path = decode_gs_uri(uri)

    bucket = storage_client.bucket(bucket)
    blob = bucket.blob(file_path)

    output_path = os.path.join(output_directory,os.path.basename(file_path))
    blob.download_to_filename(output_path)

    return output_path


def untar(tar_path):
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall()


def parse_coverage_file(file):
    total_coverage = 0
    depth_index = None

    for l,line in enumerate(file.readlines()):
        data = line.decode("utf8").strip().split('\t')

        if l == 0:
            depth_index = data.index("meandepth")
        else:
            total_coverage += float(data[depth_index])

    return total_coverage


def parse_time_as_minutes(time):
    tokens = time.split(":")

    minutes = None

    if len(tokens) == 3:
        minutes = 60*float(tokens[0]) + float(tokens[1]) + float(tokens[2])/60
    elif len(tokens) == 2:
        minutes = float(tokens[0]) + float(tokens[1])/60
    else:
        exit("ERROR: unparsable time string: " + time)

    return minutes


def parse_log_file(file):
    cpu_percent = None
    cpu_count = None
    elapsed_real_s = None
    ram_max_kbyte = None

    for l,line in enumerate(file.readlines()):
        data = line.decode("utf8").strip().split(',')

        if data[0] == "elapsed_real_s":
            elapsed_real_s = parse_time_as_minutes(data[1][:-1])
        if data[0] == "ram_max_kbyte":
            ram_max_kbyte = int(data[1][:-1])
        if data[0] == "cpu_percent":
            cpu_percent = int(data[1][:-1])
        if data[0] == "cpu_count":
            cpu_count = int(data[1])

    return elapsed_real_s, ram_max_kbyte, cpu_percent, cpu_count


def get_coverage_vs_time_for_each_tarball(tarball_uris, output_directory):
    for uri in tarball_uris[:4]:
        tar_path = download_gs_uri(uri, output_directory)

        total_coverage = None
        cpu_percent = None
        cpu_count = None
        elapsed_real_s = None
        ram_max_kbyte = None

        print(tar_path)
        with tarfile.open(tar_path, "r:gz") as tar:
            for item in tar.getmembers():
                name = os.path.basename(item.name)

                if name == "coverage.tsv":
                    f = tar.extractfile(item)
                    total_coverage = parse_coverage_file(f)

                if name == "log.csv":
                    f = tar.extractfile(item)
                    elapsed_real_s, ram_max_kbyte, cpu_percent, cpu_count = parse_log_file(f)

        print(total_coverage, elapsed_real_s, ram_max_kbyte, cpu_percent, cpu_count)


def main(tsv_paths, n_threads, output_directory):
    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for path in tsv_paths:
        df = pandas.read_table(path, sep='\t', header=0)

        n_rows, n_cols = df.shape

        for i in range(n_rows):
            bams = parse_comma_separated_string(df["bams"][i])
            bifrost_tarballs = parse_comma_separated_string(df["output_tarballs_bifrost"][i])

            n_samples = len(bams)
            n_tarballs = len(bifrost_tarballs)

            print(i, n_samples, n_tarballs)
            get_coverage_vs_time_for_each_tarball(bifrost_tarballs, output_directory)


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tsv",
        required=True,
        type=parse_comma_separated_string,
        help="Input tsvs containing fastas to be used for profiling"
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
        help="Output directory"
    )

    args = parser.parse_args()

    main(tsv_paths=args.tsv, n_threads=args.t, output_directory=args.o)
