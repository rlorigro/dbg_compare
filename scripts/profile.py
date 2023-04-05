import subprocess
import argparse
import tarfile
import shutil
import sys
import re
import os


def run_cuttlefish(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    cuttlefish_prefix = os.path.join(output_directory, "cuttlefish")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P\\n","-o",log_path]

    # cuttlefish build -s refs1.fa -k 3 -t 4 -o cdbg -w temp/ --ref
    args = time_args + ["cuttlefish", "build", "-k", str(k), "-t", str(n_threads), "--ref", "-s", fasta_path, "-o", os.path.join(output_directory, cuttlefish_prefix)]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status : FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def run_ggcat(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    ggcat_prefix = os.path.join(output_directory, "ggcat")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P\\n","-o",log_path]

    # ggcat build -e -k <k_value> -j <threads_count> <input_files> -o <output_file>
    args = time_args + ["ggcat", "build", "-e", "-k", str(k), "-j", str(n_threads), fasta_path, "-o", os.path.join(output_directory, ggcat_prefix + ".fasta")]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status : FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def run_bifrost(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    bifrost_prefix = os.path.join(output_directory, "bifrost")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P\\n","-o",log_path]
    args = time_args + ["Bifrost", "build", "-n", "-k", str(k), "-t", str(n_threads), "-r", fasta_path, "-o", os.path.join(output_directory, bifrost_prefix)]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

        # Bifrost doesn't have a proper error signal smh
        if b'Error' in p1.stderr:
            exit(p1.stderr.decode('utf8'))

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status : FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def main(tar_paths, k, graph_builder, n_cores, timeout, output_directory):
    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for tar_path in tar_paths:
        output_prefix = os.path.basename(tar_path).split('.')[0]
        output_subdirectory = os.path.join(output_directory, output_prefix)

        if not os.path.exists(output_subdirectory):
            os.makedirs(output_subdirectory)

        combined_fasta_path = os.path.join(output_directory, output_prefix + ".fasta")
        output_tsv_path = os.path.join(output_subdirectory, "coverage.tsv")

        # FASTAs for this region need to be combined
        with tarfile.open(tar_path, "r:gz") as tar, open(combined_fasta_path, 'wb') as combined_fasta:
            for item in tar.getmembers():
                if item.name.endswith(".fasta"):
                    f = tar.extractfile(item)
                    for line in f.readlines():
                        combined_fasta.write(line)

                # Also untar the coverage file
                if item.name.endswith(".tsv"):
                    f = tar.extractfile(item)

                    # This is a dumb way to extract the file, but couldn't find anything simpler
                    with open(output_tsv_path, 'wb') as out_tsv:
                        for l in f.readlines():
                            out_tsv.write(l)

        # Add specified graph outputs to subdirectory
        if graph_builder == "bifrost":
            log_path = run_bifrost(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
        elif graph_builder == "ggcat":
            log_path = run_ggcat(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
        elif graph_builder == "cuttlefish":
            log_path = run_cuttlefish(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
        else:
            exit("ERROR: unrecognized choice for graph builder")

        if log_path is not None:
            # Update the log to contain the number of processors used
            with open(log_path, 'rb+') as file:
                # Move pointer to the last char of the file (for some reason there is extra newline)
                file.seek(-1, os.SEEK_END)
                file.write(("cpu_count,%d\n" % n_cores).encode())

            # Tar the outputs: coverage TSV, log CSV, and bifrost gfa/index
            with tarfile.open(output_subdirectory + ".tar.gz", "w:gz") as tar:
                tar.add(output_subdirectory, arcname=os.path.basename(output_subdirectory))

        # Remove intermediates
        os.remove(combined_fasta_path)
        shutil.rmtree(output_subdirectory)


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'{}"))


def parse_input_string(s):
    paths = None

    if os.path.isdir(os.path.abspath(s)):
        paths = [os.path.join(s,x) for x in os.listdir(s)]
        print(paths[:10])
    else:
        paths = parse_comma_separated_string(s)

    return paths


def parse_choice(s):
    s = s.lower()

    choices = {"bifrost", "cuttlefish", "ggcat"}
    if s not in choices:
        exit("ERROR: must select one of the following tools to profile: " + str(choices))

    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tars",
        required=True,
        type=parse_input_string,
        help="Input tarballs (.tar.gz) containing fastas to be used for profiling"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="Number of cores to use"
    )

    parser.add_argument(
        "-k",
        required=False,
        default=31,
        type=int,
        help="K value to use for de bruijn graph construction"
    )

    parser.add_argument(
        "--timeout",
        required=False,
        default=60*60*24,
        type=int,
        help="Don't let graph building step run for longer than this duration (skip without error if exceeded)"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "-g",
        required=True,
        type=parse_choice,
        help="Graph building tool to use. Must be one of the following: bifrost, cuttlefish, ggcat"
    )

    args = parser.parse_args()

    main(tar_paths=args.tars, k=args.k, graph_builder=args.g, output_directory=args.o, n_cores=args.c, timeout=args.timeout)
