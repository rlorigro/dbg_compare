from module.Edge import Edge
from collections import defaultdict
import argparse


"""
>0 LN:i:33 L:+:35514:+ L:-:1:+ L:-:43315:+
>1 LN:i:31 L:+:5:+ L:+:35509:- L:+:48693:+ L:-:0:+
>2 LN:i:31 L:+:13:- L:-:10:- L:-:35586:-
>3 LN:i:31 L:+:16:+ L:+:61286:+ L:-:11:+ L:-:58289:-
>4 LN:i:35 L:+:35516:+ L:-:16:- L:-:58288:-
>5 LN:i:33 L:+:15:- L:-:1:-
"""
def main(fasta_path, no_sequence=False):
    output_path = '.'.join(fasta_path.split('.')[:-1]+["gfa"])
    print(fasta_path)
    print(output_path)

    nodes = defaultdict(str)
    edges = set()

    id = None
    sequence = None

    with open(fasta_path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == '>':
                id, line_edges = Edge.parse_bcalm_string(line[1:])

                for e in line_edges:
                    e.canonicalize()
                    edges.add(e.to_gfa_line())

                if l > 0:
                    nodes[id] = sequence
                sequence = ""

            else:
                sequence += line.strip()

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
        help="Input fasta containing ggcat formated sequences, with BCALM Link strings as annotation"
    )

    parser.add_argument(
        "--no_sequence",
        required=True,
        type=str_as_bool,
        help="Don't write any sequence to the GFA"
    )

    args = parser.parse_args()

    main(fasta_path=args.i, no_sequence=args.no_sequence)
