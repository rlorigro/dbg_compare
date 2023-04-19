import subprocess
import argparse
import tarfile
import random
import shutil
import sys
import re
import os


def dry_run(output_directory):
    log_path = os.path.join(output_directory, "log.csv")

    with open(log_path,'w') as file:
        file.write("elapsed_real_s,0:0.0\nelapsed_kernel_s,0:0.0\nram_max_kbyte,0\nram_avg_kbyte,0\ncpu_percent,0%\n")

    return log_path


def run_cuttlefish(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    cuttlefish_prefix = os.path.join(output_directory, "cuttlefish")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P","-o",log_path]

    # cuttlefish build -s refs1.fa -k 3 -t 4 -o cdbg -w temp/ --ref
    args = time_args + ["cuttlefish", "build", "-k", str(k), "-t", str(n_threads), "--ref", "-s", fasta_path, "-o", os.path.join(output_directory, cuttlefish_prefix)]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def run_ggcat(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    ggcat_prefix = os.path.join(output_directory, "ggcat")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P","-o",log_path]

    # ggcat build -e --min-multiplicity 1 -k <k_value> -j <threads_count> <input_files> -o <output_file>
    args = time_args + ["ggcat", "build", "-e", "--min-multiplicity", "1", "-k", str(k), "-j", str(n_threads), fasta_path, "-o", os.path.join(output_directory, ggcat_prefix + ".fasta")]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def run_bifrost(fasta_path, k, output_directory, n_threads, timeout=60*60*24):
    log_path = os.path.join(output_directory, "log.csv")
    bifrost_prefix = os.path.join(output_directory, "bifrost")

    time_args = ["/usr/bin/time","-f","elapsed_real_s,%E\\nelapsed_kernel_s,%S\\nram_max_kbyte,%M\\nram_avg_kbyte,%t\\ncpu_percent,%P","-o",log_path]
    args = time_args + ["Bifrost", "build", "-n", "-k", str(k), "-t", str(n_threads), "-r", fasta_path, "-o", os.path.join(output_directory, bifrost_prefix)]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

        # Bifrost doesn't have a proper error signal smh
        if b'Error' in p1.stderr:
            exit(p1.stderr.decode('utf8'))

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return log_path


def main(tar_paths, k, graph_builder, n_cores, timeout, n_samples, output_directory):
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

        n = 0
        samples_visited = set()
        tsv_lines_per_sample = dict()
        all_empty = True

        # FASTAs for this region need to be combined
        with tarfile.open(tar_path, "r:gz") as tar, open(combined_fasta_path, 'wb') as combined_fasta:
            for item in tar.getmembers():
                if item.name.endswith(".fasta"):

                    # Arbitrarily sample the top n in the list
                    n += 1
                    if n > n_samples:
                        continue

                    samples_visited.add(os.path.basename(item.name).split('.')[0])

                    f = tar.extractfile(item)
                    for l,line in enumerate(f.readlines()):
                        combined_fasta.write(line)

                        if l > 0:
                            all_empty = False

                # Also untar the coverage file
                if item.name.endswith(".tsv"):
                    f = tar.extractfile(item)

                    # If we are subsampling, need to index the lines by sample so only relevant ones can be copied later
                    for l,line in enumerate(f.readlines()):
                        line = line.decode('utf8')
                        tokens = line.split('\t')

                        name = None
                        if l == 0:
                            name = "header"
                        else:
                            name = tokens[0]

                        tsv_lines_per_sample[name] = line

        # Write the coverage data for the samples that were visited
        with open(output_tsv_path, 'w') as out_tsv:
            out_tsv.write(tsv_lines_per_sample["header"])
            for s in samples_visited:
                out_tsv.write(tsv_lines_per_sample[s])

        # Add specified graph outputs to subdirectory
        if all_empty:
            # Don't bother trying to run the graph tool, it may crash on empty FASTA
            # Use dryrun function to generate 0 for all log stats
            log_path = dry_run(output_subdirectory)
            sys.stderr.write("WARNING: no coverage for region %s\n" % output_prefix)
        else:
            if graph_builder == "bifrost":
                log_path = run_bifrost(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
            elif graph_builder == "ggcat":
                log_path = run_ggcat(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
            elif graph_builder == "cuttlefish":
                log_path = run_cuttlefish(combined_fasta_path, k, output_subdirectory, n_cores, timeout=timeout)
            elif graph_builder == "test":
                log_path = dry_run(output_subdirectory)
            else:
                exit("ERROR: unrecognized choice for graph builder")

        if log_path is not None:
            # Update the log to contain the number of processors used
            with open(log_path, 'a') as file:
                file.write("cpu_count,%d\n" % n_cores)

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

    choices = {"bifrost", "cuttlefish", "ggcat", "test"}
    if s not in choices:
        exit("ERROR: must select one of the following tools to profile: " + str(choices))

    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tars",
        required=True,
        type=parse_input_string,
        help="Input tarballs (.tar.gz) containing fastas to be used for profiling. Can be a comma separated list OR a directory"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="Number of cores to use"
    )

    parser.add_argument(
        "-n",
        required=False,
        default=None,
        type=int,
        help="How many samples to use for profiling"
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

    main(tar_paths=args.tars, k=args.k, graph_builder=args.g, output_directory=args.o, n_cores=args.c, n_samples=args.n, timeout=args.timeout)
