from module.Edge import Edge
from collections import defaultdict
import argparse
import sys
import os


"""
>0 LN:i:33 L:+:35514:+ L:-:1:+ L:-:43315:+
>1 LN:i:31 L:+:5:+ L:+:35509:- L:+:48693:+ L:-:0:+
>2 LN:i:31 L:+:13:- L:-:10:- L:-:35586:-
>3 LN:i:31 L:+:16:+ L:+:61286:+ L:-:11:+ L:-:58289:-
>4 LN:i:35 L:+:35516:+ L:-:16:- L:-:58288:-
>5 LN:i:33 L:+:15:- L:-:1:-
"""
def main(fasta_path, output_path, no_sequence=False):
    output_directory = os.path.dirname(output_path)

    if not len(output_directory) == 0:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    if not output_path.endswith(".gfa"):
        exit("ERROR: output path does not have GFA suffix: " + output_path)

    print(fasta_path)
    print(output_path)

    nodes = defaultdict(str)
    edges = set()

    id = None
    sequence = None

    if not os.path.exists(fasta_path):
        sys.stderr.write("WARNING: fasta file not found, terminating early: %s" % fasta_path)
        return

    with open(fasta_path, 'r') as file:
        for l,line in enumerate(file):
            if len(line) <= 1:
                sys.stderr.write("WARNING: empty line detected at l=%d in file: %s" % (l, fasta_path))
                continue

            if line[0] == '>':
                if l > 0:
                    nodes[id] = sequence

                id, line_edges = Edge.parse_bcalm_string(line[1:])

                for e in line_edges:
                    e.canonicalize()
                    edges.add(e.to_gfa_line())

                sequence = ""

            else:
                sequence += line.strip()

    # Final line
    if sequence is not None and id is not None:
        nodes[id] = sequence

    with open(output_path, 'w') as file:
        file.write("H\tVN:Z:1.0\n")
        for id,seq in nodes.items():
            file.write("S\t")
            file.write(id)
            file.write('\t')

            if no_sequence:
                file.write("*\tLN:i:")
                file.write(str(len(seq)))
                file.write('\n')
            else:
                file.write(seq)
                file.write('\n')

        for item in edges:
            file.write(item)
            file.write('\n')


def str_as_bool(s):
    if s in {'Y','y','1','true','True','on','yes'}:
        return True
    elif s in {'N','n','0','false','False','off','no'}:
        return False
    else:
        exit("ERROR: unparsable boolean string: %s" % s)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input fasta containing ggcat formatted sequences, with BCALM Link strings as annotation"
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output path, any non-existent directories will be created"
    )

    parser.add_argument(
        "--no_sequence",
        required=True,
        type=str_as_bool,
        help="Don't write any sequence to the GFA"
    )

    args = parser.parse_args()

    main(fasta_path=args.i, output_path=args.o, no_sequence=args.no_sequence)
