"""
schematic.py -- draw the schematic of a series/parallel RLC one-port with
schemdraw, straight from its m-expression.

    >>> from mexpr_re import MExpr
    >>> from schematic import draw
    >>> draw(MExpr('B s (C p (B s A))'))         # displayed inline in Jupyter

Atoms are drawn as   B -> resistor R,  C -> inductor L,  A -> capacitor C.
Every element carries two labels: its circuit name above (R1, L1, C1) and,
in grey below, the label used by MExpr in the fault tables (B1, C1, A1).
A collapsed network ('K' or 'P', e.g. a residual network after a fault) is
drawn as a plain wire (short) or as two open terminals (open).
"""

import schemdraw
import schemdraw.elements as elm

SYMBOL = {'B': elm.ResistorIEC, 'C': elm.Inductor2, 'A': elm.Capacitor}
NAME = {'B': 'R', 'C': 'L', 'A': 'C'}

ATOM_W, ATOM_H = 3.0, 2.4          # one element with its labels (schemdraw units)
LEAD = 0.8                         # wire on each side of a parallel block
PORT = 1.0                         # wire from each port terminal


# --------------------------------------------------------------------
#  'B1 s (C1 p (B2 s A1))'  ->  tree
# --------------------------------------------------------------------
def parse(text):
    """-> 'B1' (a labelled atom) | ('s' | 'p', [subtrees])"""
    tokens = text.replace('(', ' ( ').replace(')', ' ) ').split()
    pos = 0

    def expr():
        nonlocal pos
        items, ops = [term()], set()
        while pos < len(tokens) and tokens[pos] in ('s', 'p'):
            ops.add(tokens[pos])
            pos += 1
            items.append(term())
        if len(ops) > 1:
            raise ValueError(f"mixed s and p without parentheses in '{text}'")
        return items[0] if not ops else (ops.pop(), items)

    def term():
        nonlocal pos
        tok = tokens[pos]
        pos += 1
        if tok == '(':
            node = expr()
            pos += 1                                # ')'
            return node
        return tok

    return expr()


def size(node):
    """(width, height) of a subtree."""
    if isinstance(node, str):
        return ATOM_W, ATOM_H
    sizes = [size(c) for c in node[1]]
    if node[0] == 's':
        return sum(w for w, _ in sizes), max(h for _, h in sizes)
    return max(w for w, _ in sizes) + 2 * LEAD, sum(h for _, h in sizes)


def _place(d, node, x, y):
    """Draw the subtree between (x, y) and (x + width, y)."""
    w, h = size(node)
    if isinstance(node, str):                       # a labelled atom, e.g. 'B2'
        kind, idx = node[0], node[1:]
        d.add(SYMBOL[kind]().endpoints((x, y), (x + w, y))
              .label(f'${NAME[kind]}_{{{idx}}}$', loc='top')
              .label(f'${kind}_{{{idx}}}$', loc='bottom', color='gray', fontsize=10))
        return
    if node[0] == 's':
        for child in node[1]:
            _place(d, child, x, y)
            x += size(child)[0]
        return
    xl, xr = x + LEAD, x + w - LEAD                 # parallel: branches top to bottom
    top = y + h / 2
    ys = []
    for child in node[1]:
        cw, ch = size(child)
        yc = top - ch / 2
        ys.append(yc)
        _place(d, child, xl, yc)
        if xl + cw < xr - 1e-9:
            d.add(elm.Line().endpoints((xl + cw, yc), (xr, yc)))
        top -= ch
    d.add(elm.Line().endpoints((xl, ys[0]), (xl, ys[-1])))
    d.add(elm.Line().endpoints((xr, ys[0]), (xr, ys[-1])))
    d.add(elm.Line().endpoints((x, y), (xl, y)))
    d.add(elm.Line().endpoints((xr, y), (x + w, y)))


def draw(network):
    """
    A schemdraw Drawing of the network; `network` is an MExpr or its text
    with labelled atoms ('B1 s (C1 p A1)').  In Jupyter the drawing is
    shown when it is the last value of a cell (or pass it to display()).
    """
    text = str(network)
    d = schemdraw.Drawing(show=False)
    d.config(fontsize=12)
    if text in ('K', 'P'):                          # the whole network collapsed
        d.add(elm.Dot(open=True).at((0, 0)))
        if text == 'K':
            d.add(elm.Line().endpoints((0, 0), (ATOM_W + 2 * PORT, 0)).label('short', loc='top'))
        else:
            d.add(elm.Gap().endpoints((0, 0), (ATOM_W + 2 * PORT, 0)).label('open', loc='top'))
        d.add(elm.Dot(open=True).at((ATOM_W + 2 * PORT, 0)))
        return d
    tree = parse(text)
    w, _ = size(tree)
    d.add(elm.Dot(open=True).at((0, 0)))
    d.add(elm.Line().endpoints((0, 0), (PORT, 0)))
    _place(d, tree, PORT, 0)
    d.add(elm.Line().endpoints((PORT + w, 0), (2 * PORT + w, 0)))
    d.add(elm.Dot(open=True).at((2 * PORT + w, 0)))
    return d
