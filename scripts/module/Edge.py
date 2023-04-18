
class Edge:
    def __init__(self,id_a,reversal_a,id_b,reversal_b):
        self.id_a = id_a
        self.reversal_a = reversal_a
        self.id_b = id_b
        self.reversal_b = reversal_b

    def to_gfa_line(self):
        return "L\t%s\t%s\t%s\t%s\t*" % (self.id_a, self.reversal_as_char(self.reversal_a), self.id_b, self.reversal_as_char(self.reversal_b))

    """
    Flip the edge to its equivalent but reversed representation
    """
    def flip(self):
        id_a_flipped = self.id_b
        id_b_flipped = self.id_a
        reversal_a_flipped = not self.reversal_b
        reversal_b_flipped = not self.reversal_a

        self.id_a = id_a_flipped
        self.id_b = id_b_flipped
        self.reversal_a = reversal_a_flipped
        self.reversal_b = reversal_b_flipped

        return

    """
    Order the edge so that the lexicographically lower ID is first
    """
    def canonicalize(self):
        if not self.is_canonical():
            self.flip()

        return

    def is_canonical(self):
        return self.id_a < self.id_b

    def __str__(self):
        return self.to_gfa_line()

    @staticmethod
    def parse_bcalm_string(s):
        tokens = s.strip().split(' ')

        id_a = tokens[0]
        edges = list()

        edge_tokens = [t[2:] for t in tokens if t.startswith("L:")]

        for t in edge_tokens:
            # Parse what remains of the token, now appearing like this: "+:35514:+"
            reversal_a = Edge.parse_char_as_reversal(t[0])
            reversal_b = Edge.parse_char_as_reversal(t[-1])
            id_b = t[2:-2]

            e = Edge(id_a,reversal_a,id_b,reversal_b)
            edges.append(e)

        return id_a, edges

    @staticmethod
    def parse_char_as_reversal(c):
        if c == '+':
            return False
        elif c == '-':
            return True
        else:
            exit("ERROR: unparsable reversal char: %s" % c)

    @staticmethod
    def reversal_as_char(reversal):
        if reversal:
            return '-'
        else:
            return '+'


def test_edge():
    bcalm_lines = [
        "0 LN:i:33 L:+:35514:+ L:-:1:+ L:-:43315:+",
        "1 LN:i:31 L:+:5:+ L:+:35509:- L:+:48693:+ L:-:0:+",
        "2 LN:i:31 L:+:13:- L:-:10:- L:-:35586:-",
        "3 LN:i:31 L:+:16:+ L:+:61286:+ L:-:11:+ L:-:58289:-",
        "4 LN:i:35 L:+:35516:+ L:-:16:- L:-:58288:-",
        "5 LN:i:33 L:+:15:- L:-:1:-",
    ]

    for l in bcalm_lines:
        node, edges = Edge.parse_bcalm_string(l)

        print(l)
        for e in edges:
            print(e)
            if not e.is_canonical():
                e.canonicalize()
                print(" --> %s" % str(e))


if __name__ == "__main__":
    test_edge()
