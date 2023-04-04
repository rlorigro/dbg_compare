from multiprocessing import Pool
import argparse
import tarfile

import numpy
import pandas
import math
import sys
import re
import os

from matplotlib.lines import Line2D
from matplotlib import pyplot
from matplotlib import colors
from google.cloud import storage


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap


def decode_gs_uri(uri):
    tokens = uri.split('/')
    bucket, file_path = tokens[2], '/'.join(tokens[3:])
    return bucket, file_path


def download_gs_uri(uri, output_directory, cache=True):
    bucket, file_path = decode_gs_uri(uri)
    output_path = os.path.join(output_directory,os.path.basename(file_path))

    if (not os.path.exists(output_path)) or (not cache):
        sys.stderr.write("Downloading: %s\n" % output_path)
        sys.stderr.flush()
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


def get_resource_stats_for_each_tarball(tarball_paths):
    for tar_path in tarball_paths:
        yield get_resource_stats_for_tarball(tar_path)


def get_resource_stats_for_tarball(tar_path):
    total_coverage = None
    cpu_percent = None
    cpu_count = None
    elapsed_real_s = None
    ram_max_kbyte = None
    ram_max_mbyte = None

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

    return total_coverage, elapsed_real_s, ram_max_mbyte, adjusted_cpu_percent


def main(tsv_path, n_threads, output_directory):
    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    tool_names = ["bifrost", "ggcat", "cuttlefish"]
    fig,axes = pyplot.subplots(nrows=2,ncols=2)

    df = pandas.read_table(tsv_path, sep='\t', header=0)

    n_rows, n_cols = df.shape

    print(n_rows)

    colormaps = {
        "bifrost":pyplot.colormaps.get_cmap("Greens"),
        "ggcat":truncate_colormap(pyplot.colormaps.get_cmap("Blues"), maxval=0.8),
        "cuttlefish":pyplot.colormaps.get_cmap("Purples")
    }

    coverage_colormap = pyplot.colormaps.get_cmap("gist_heat")
    max_coverage = 0

    for n,name in enumerate(tool_names):
        print("---- %s ----" % name)
        colormap = colormaps[name]

        for i in range(n_rows):
            total_coverage = list()
            elapsed_real_s = list()
            ram_max_mbyte = list()
            cpu_percent = list()

            row_name = df.iloc[i][0]

            print(row_name)

            if "100Kbp" not in row_name:
                print("Skipping")
                continue

            bams = parse_comma_separated_string(df.iloc[i]["bams"])

            try:
                tarballs = parse_comma_separated_string(df.iloc[i]["output_tarballs_"+name])[:]
            except Exception as e:
                print(e)
                continue

            n_samples = len(bams)

            print(n_samples, len(tarballs))

            # Just use the same color for all dots within a dbg tool
            color = colormap(1.0)

            # Each tool downloads its regions to its own subdirectory to prevent overwriting (filenames are by region)
            output_subdirectory = os.path.join(output_directory, name)
            output_subdirectory = os.path.join(output_subdirectory, str(n_samples))
            if not os.path.exists(output_subdirectory):
                os.makedirs(output_subdirectory)

            # Multithread the downloading of files
            args = [[t,output_subdirectory] for t in tarballs]
            with Pool(n_threads) as pool:
                download_results = pool.starmap(download_gs_uri, args)

            # Multithread the parsing of results
            args = [[str(x)] for x in download_results]
            with Pool(n_threads) as pool:
                stats_results = pool.starmap(get_resource_stats_for_tarball, args)

            # Aggregate
            for stats in stats_results:
                total_coverage.append(stats[0])
                elapsed_real_s.append(stats[1])
                ram_max_mbyte.append(stats[2])
                cpu_percent.append(stats[3])

                if stats[0] > max_coverage:
                    max_coverage = stats[0]

            axes[0][0].scatter(x=total_coverage, y=elapsed_real_s, s=0.7, color=color, alpha=0.5)
            axes[0][1].scatter(x=total_coverage, y=ram_max_mbyte, s=0.7, color=color, alpha=0.5)
            axes[1][0].scatter(x=total_coverage, y=cpu_percent, s=0.7, color=color, alpha=0.5)

            # Only plot coverage histogram once
            if n == 0:
                n_bins = 400
                max = 1000
                step_size = float(max+1)/n_bins
                bins = numpy.arange(0,max,step_size)
                histogram,_ = numpy.histogram(total_coverage, bins=bins)

                bin_max = numpy.argmax(histogram)
                x_max = step_size*bin_max
                y_max = histogram[bin_max]

                v = float(math.log2(n_samples)-1)/5

                axes[1][1].plot(bins[:-1], histogram, color=coverage_colormap(v))
                axes[1][1].text(x_max, y_max, str(n_samples), horizontalalignment='left', verticalalignment='bottom')

    axes[0][0].set_xlabel("Average depth")
    axes[0][1].set_xlabel("Average depth")
    axes[1][0].set_xlabel("Average depth")
    axes[1][1].set_xlabel("Average depth")

    axes[0][0].set_ylabel("Run time (min)")
    axes[0][1].set_ylabel("Peak RAM (MB)")
    axes[1][0].set_ylabel("CPU usage (% of available cores)")
    axes[1][1].set_ylabel("Frequency")

    y_min, y_max = axes[1][1].get_ylim()
    axes[1][1].set_ylim([y_min, y_max*1.1])

    axes[0][0].set_xlim([0,600])
    axes[0][1].set_xlim([0,600])
    axes[1][0].set_xlim([0,600])
    axes[1][1].set_xlim([0,600])

    colors = ["green", "blue", "purple"]
    custom_lines = list()
    for i in range(len(colors)):
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))

    axes[0][1].legend(custom_lines, tool_names, bbox_to_anchor=(1.1, 1))

    fig.set_size_inches(12,9)
    pyplot.tight_layout()

    pyplot.savefig("resource_usage.png",dpi=200)

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
