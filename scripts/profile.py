from multiprocessing import Pool

import subprocess
import argparse
import tarfile
import sys
import re
import os


def run_bifrost(fasta_path, n_threads, output_directory):
    output_path = os.path.join(output_directory, fasta_path.replace("fasta", "gfa"))
    output_path = output_path.replace(".bam", "_coverage.tsv")

    print(output_path)

    time_args = ["/usr/bin/time","-f","elapsed_real,%E\\nelapsed_kernel,%S\\nram_max,%M\\nram_avg,%t\\ncpu_percent,%P\\n","-o","test.log"]
    args = time_args + ["Bifrost", "build", "-t", str(n_threads), "-r", fasta_path]

    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return None

    return output_path


def main(tar_paths, output_directory, n_cores):
    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for tar_path in tar_paths:
        output_prefix = os.path.basename(tar_path).split('.')[0]
        fasta_path = os.path.join(output_directory, output_prefix + ".fasta")

        with tarfile.open(tar_path, "r:gz") as tar, open(fasta_path, 'wb') as combined_fasta:
            for item in tar.getmembers():
                print(item.name)
                if item.name.endswith(".fasta"):
                    f = tar.extractfile(item)
                    for line in f.readlines():
                        combined_fasta.write(line)

                if item.name.endswith(".tsv"):
                    f = tar.extractfile(item)
                    for line in f.readlines():
                        print(line.decode('utf8').strip())

        run_bifrost(fasta_path, n_cores, output_directory)


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tars",
        required=True,
        type=parse_comma_separated_string,
        help="Input tarballs (.tar.gz) containing fastas to be used for profiling"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="Number of cores to use"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(tar_paths=args.tars, output_directory=args.o, n_cores=args.c)
