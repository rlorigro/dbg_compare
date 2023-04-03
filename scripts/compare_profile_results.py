from multiprocessing import Pool
import argparse
import tarfile
import pandas
import math
import re
import os

from matplotlib.lines import Line2D
from matplotlib import pyplot
from google.cloud import storage


def decode_gs_uri(uri):
    tokens = uri.split('/')
    bucket, file_path = tokens[2], '/'.join(tokens[3:])
    return bucket, file_path


def download_gs_uri(uri, output_directory, cache=True):
    bucket, file_path = decode_gs_uri(uri)
    output_path = os.path.join(output_directory,os.path.basename(file_path))

    if (not os.path.exists(output_path)) or (not cache):
        storage_client = storage.Client()

        bucket = storage_client.bucket(bucket)
        blob = bucket.blob(file_path)

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


def get_resource_stats_for_each_tarball(tarball_paths, output_directory):
    for tar_path in tarball_paths:
        # tar_path = download_gs_uri(uri, output_directory)

        total_coverage = None
        cpu_percent = None
        cpu_count = None
        elapsed_real_s = None
        ram_max_kbyte = None
        ram_max_mbyte = None

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
                    ram_max_mbyte = float(ram_max_kbyte)/1000

        # Normalize CPU percent so it shows percent of total CPUs, instead of e.g. 233%
        adjusted_cpu_percent = cpu_percent / cpu_count
        print(cpu_percent, cpu_count, adjusted_cpu_percent)

        yield total_coverage, elapsed_real_s, ram_max_mbyte, adjusted_cpu_percent


def main(tsv_path, n_threads, output_directory):
    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    tool_names = ["bifrost", "ggcat", "cuttlefish"]
    fig,axes = pyplot.subplots(nrows=2,ncols=2)

    df = pandas.read_table(tsv_path, sep='\t', header=0)

    print(df.columns.values)

    n_rows, n_cols = df.shape

    colormap_name = {
        "bifrost":"Greens",
        "ggcat":"Blues",
        "cuttlefish":"Oranges"
    }

    for name in tool_names:
        colormap = pyplot.colormaps.get_cmap(colormap_name[name])

        for i in range(n_rows):
            total_coverage = list()
            elapsed_real_s = list()
            ram_max_mbyte = list()
            cpu_percent = list()

            bams = parse_comma_separated_string(df["bams"][i])

            try:
                tarballs = parse_comma_separated_string(df["output_tarballs_"+name][i])
            except Exception as e:
                print(e)
                continue

            n_samples = len(bams)

            color_index = float(math.log2(n_samples))/(n_rows+1)
            color = colormap(color_index)

            # Each tool downloads its regions to its own subdirectory to prevent overwriting (filenames are by region)
            output_subdirectory = os.path.join(output_directory, name)
            output_subdirectory = os.path.join(output_subdirectory, str(n_samples))
            if not os.path.exists(output_subdirectory):
                os.makedirs(output_subdirectory)

            args = [[t,output_subdirectory] for t in tarballs]
            print(len(tarballs))

            with Pool(n_threads) as pool:
                results = pool.starmap(download_gs_uri, args)

            for stats in get_resource_stats_for_each_tarball(results, output_subdirectory):
                total_coverage.append(stats[0])
                elapsed_real_s.append(stats[1])
                ram_max_mbyte.append(stats[2])
                cpu_percent.append(stats[3])

            axes[0][0].scatter(x=total_coverage, y=elapsed_real_s, s=0.5, color=color)
            axes[0][1].scatter(x=total_coverage, y=ram_max_mbyte, s=0.5, color=color)
            axes[1][0].scatter(x=total_coverage, y=cpu_percent, s=0.5, color=color)

    axes[0][0].set_xlabel("Average depth")
    axes[0][1].set_xlabel("Average depth")
    axes[1][0].set_xlabel("Average depth")
    axes[0][0].set_ylabel("Run time (min)")
    axes[0][1].set_ylabel("Peak RAM (MB)")
    axes[1][0].set_ylabel("CPU usage (% of available cores)")
    axes[1][1].axis("off")

    colors = ["green", "blue", "orange"]
    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    axes[0][1].legend(custom_lines, tool_names, bbox_to_anchor=(1.1, 1))

    pyplot.tight_layout()

    pyplot.show()
    pyplot.close()


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="Input tsv containing URIs, linking to tarballs to be parsed"
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

    main(tsv_path=args.tsv, n_threads=args.t, output_directory=args.o)
