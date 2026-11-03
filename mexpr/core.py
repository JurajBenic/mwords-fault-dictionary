"""
mexpr.py -- the structural fault dictionary of a series/parallel RLC one-port,
computed from its m-expression (no element values, no simulation).

    python mexpr.py "B p (B s C) p C"

An m-expression is a text such as 'B s (C p (B s A))': atoms A (capacitor),
B (resistor), C (inductor) joined by
    s ... series connection      (MNumber &)
    p ... parallel connection    (MNumber |)
Mixed s and p on one level need parentheses: 'B s (A p C)', not 'B s A p C'.
Every atom gets a label (B1, B2, C1, ...) so that a reduced network can say
WHICH element survived.  Two ideal atoms may appear in a reduced network,
    K ... short circuit (Z = 0)        P ... open circuit (Z = inf)
with the absorption rules   x s K = x,  x s P = P,  x p P = x,  x p K = K.

Everything is computed on the text itself, from the innermost parentheses
outwards (_levels, _fold).  The same walk gives the m-number of the network
(atoms = MNumber) and any reduced network (atoms = text, K, P).

What the dictionary contains (Sections 2-4 of the paper):
    word        canonical form of the m-word (capital irreducible reduct)
    quadruple   Phi = (v(l), lambda_s, lambda_p, v(r))  -- the signature
    corners     delta_0 = x[A->P, C->K], delta_inf = x[A->K, C->P],
                delta_K = x[A,C->K], delta_P = x[A,C->P]
    faults      every element shorted (K) and opened (P): residual network,
                its word and signature, which test detects it
                (DC: v(l) changes, HF: v(r) changes, sweep: lambda changes,
                port: the port itself is shorted/opened); a fault whose
                signature equals the nominal one is structurally undetectable
                (none of the three tests detects it for any element values;
                value-dependent features of the response may still reveal it)
    groups      structural ambiguity group = same residual network after
                merging of same-type neighbours (same topology and element
                types, possibly different values);
                absolute ambiguity group = the sub-group of a structural
                group formed by all opens of one series node or all shorts
                of one parallel node (identical residual networks, element
                values included);
                signature ambiguity group = same quadruple

    >>> e = MExpr('B s (C p (B s A))')                  # network #63
    >>> print(e, e.word, e.quadruple)
    B1 s (C1 p (B2 s A1)) BCAB (0, 2, 2, 0)
    >>> e.corners()
    {'delta_0': 'B1', 'delta_inf': 'B1 s B2', 'delta_K': 'B1', 'delta_P': 'P'}
    >>> for f in e.faults()[:3]: print(f)
    Fault(element='B1', kind='K', residual='C1 p (B2 s A1)', word='CAB', quadruple=(1, 2, 1, 0), test='DC+sweep', undetectable=False)
    Fault(element='B1', kind='P', residual='P', word='P', quadruple=None, test='port', undetectable=False)
    Fault(element='C1', kind='K', residual='B1', word='B', quadruple=(0, 0, 0, 0), test='sweep', undetectable=False)
    >>> e.ambiguity_groups()
    {'structural': [['B2 P', 'A1 P']], 'absolute': [['B2 P', 'A1 P']], 'signature': [['nominal', 'B2 K'], ['B2 P', 'A1 P']]}
"""

import re
import sys
from collections import Counter, namedtuple
from functools import lru_cache
from mnumbers import MNumber

MNumber.set_alphabet('abc', m_continuum=True)       # v(a) = -1, v(b) = 0, v(c) = 1

GROUP = re.compile(r'\(([^()]*)\)')                 # innermost parentheses
UNIT = {'s': 'K', 'p': 'P'}                         # x s K = x,  x p P = x
KILLER = {'s': 'P', 'p': 'K'}                       # x s P = P,  x p K = K
REACTIVE = {'A', 'C'}

Fault = namedtuple('Fault', 'element kind residual word quadruple test undetectable')
# undetectable = structurally undetectable: the signature (quadruple) is unchanged by the fault


class MExprError(Exception):
    pass


# --------------------------------------------------------------------
#  Text helpers
# --------------------------------------------------------------------
def _words(text):
    """'B s (C p A)' -> ['B', 's', '(', 'C', 'p', 'A', ')']"""
    return text.replace('(', ' ( ').replace(')', ' ) ').split()


def _text(words):
    return ' '.join(words).replace('( ', '(').replace(' )', ')')


def _is_atom(word):
    return word[0] in 'ABC'


def _levels(expr):
    """
    The levels of an expression, innermost first; every group of
    parentheses becomes a placeholder #0, #1, ...

    >>> _levels('B p (C s (B p A))')
    ['B p A', 'C s #0', 'B p #1']
    """
    out = []
    while (m := GROUP.search(expr)):
        expr = expr[:m.start()] + f'#{len(out)}' + expr[m.end():]
        out.append(m.group(1))
    return out + [expr]


def _fold(expr, leaf, combine):
    """
    Compute an expression from the inside out:
        leaf(word)          ... the value of an atom (or K, P),
        combine(op, values) ... the value of one level 'x op y op z'.
    """
    values = {}
    for n, level in enumerate(_levels(expr)):
        words = level.split()
        op = 'p' if 'p' in words else 's'
        values[f'#{n}'] = combine(op, [values.get(w) or leaf(w) for w in words[::2]])
    return values[f'#{n}']


def _top_op(text):
    """Operator of the outermost level, or None for a single atom."""
    top = _levels(text)[-1].split()
    return 'p' if 'p' in top else 's' if 's' in top else None


# ----- the three meanings of one level 'x op y op z'
def _connect(op, atoms):
    """m-numbers: & in series, | in parallel."""
    result = atoms[0]
    for x in atoms[1:]:
        result = (result & x) if op == 's' else (result | x)
    return result


def _join(op, items):
    """
    Texts, with the absorption rules of K and P.

    >>> _join('s', ['B1', 'K', 'B2 p A1']), _join('s', ['B1', 'P']), _join('p', ['P', 'P'])
    ('B1 s (B2 p A1)', 'P', 'P')
    """
    if KILLER[op] in items:
        return KILLER[op]
    items = [t for t in items if t != UNIT[op]]
    if not items:
        return UNIT[op]
    if len(items) == 1:
        return items[0]
    return f' {op} '.join(t if _top_op(t) in (None, op) else f'({t})' for t in items)


def _shape(op, items):
    """
    The network up to renaming and merging of same-type neighbours
    (R s R = R): nested sorted tuples, the key of an absolute group.

    >>> _fold('B1 s B2 s (A1 p B3)', lambda w: w[0], _shape)
    ('s', ('B', ('p', ('A', 'B'))))
    """
    flat = []
    for it in items:
        flat.extend(it[1] if isinstance(it, tuple) and it[0] == op else [it])
    atoms = sorted({x for x in flat if isinstance(x, str)})
    nodes = sorted((x for x in flat if isinstance(x, tuple)), key=repr)
    out = atoms + nodes
    return out[0] if len(out) == 1 else (op, tuple(out))


# --------------------------------------------------------------------
#  Words: canonical form and equality as m-numbers
# --------------------------------------------------------------------
@lru_cache(maxsize=None)
def _reducts(w):
    """Irreducible words reachable from w by prst -> pt with equal lambda."""
    out, found = set(), False
    for i in range(len(w) - 3):
        s = w[i:i + 4]
        if MNumber(s).lmb == MNumber(s[0] + s[3]).lmb:
            found = True
            out |= _reducts(w[:i] + s[0] + s[3] + w[i + 4:])
    return frozenset(out) if found else frozenset([w])


def canonical(m):
    """
    The canonical form of an m-word: its unique irreducible reduct written in
    capital letters only (Section 2.3).

    >>> canonical(MNumber('BBcbCCacABbAaBBB'))
    'BCAB'
    """
    caps = [w for w in _reducts(m.m_word) if all(w[i] == w[i + 1] for i in range(0, len(w), 2))]
    if len(caps) != 1:
        raise MExprError(f'no unique capital reduct of {m.l2U}')
    return MNumber(caps[0]).l2U


def word_class(w):
    """
    All capital irreducible words equal to w as m-numbers: Sare's reversal
    rule p u p = p u<- p applied to any pair of equal letters, then
    compression.

    >>> sorted(word_class('BACB')), sorted(word_class('BAB'))
    (['BACB', 'BCAB'], ['BAB'])
    """
    seen, todo = {w}, [w]
    while todo:
        u = todo.pop()
        for i in range(len(u)):
            for j in range(i + 2, len(u)):
                if u[i] == u[j]:
                    for red in _reducts(MNumber(u[:i + 1] + u[i + 1:j][::-1] + u[j:]).m_word):
                        if all(red[k] == red[k + 1] for k in range(0, len(red), 2)):
                            cap = MNumber(red).l2U
                            if cap not in seen:
                                seen.add(cap); todo.append(cap)
    return frozenset(seen)


# --------------------------------------------------------------------
#  MExpr: the network and its dictionary
# --------------------------------------------------------------------
class MExpr:
    """A series/parallel network of atoms A, B, C given as a text."""

    def __init__(self, text):
        counts = Counter()
        words = []
        for w in _words(text):
            if _is_atom(w) and w.isalpha():                 # unlabelled atom
                counts[w] += 1
                w = f'{w}{counts[w]}'
            words.append(w)
        self.expr = _text(words)                            # 'B1 s (C1 p (B2 s A1))'
        for level in _levels(self.expr):
            if ' s ' in level and ' p ' in level:
                raise MExprError(f"Mixed s and p need parentheses: '{level}'.")
            if not re.fullmatch(r'(#\d+|[ABC]\d+|[KP])( [sp] (#\d+|[ABC]\d+|[KP]))*', level):
                raise MExprError(f"Bad expression '{self.expr}'.")

    def __str__(self):
        return self.expr

    @property
    def labels(self):
        return [w for w in _words(self.expr) if _is_atom(w)]

    @property
    def ideal(self):
        """'K' or 'P' if the whole network collapsed, else None."""
        return self.expr if self.expr in ('K', 'P') else None

    @property
    def mnumber(self):
        """The raw m-word of the network (8n-6 letters), by MNumber's & and |."""
        return _fold(self.expr, lambda w: MNumber(w[0]), _connect)

    @property
    def word(self):
        """Canonical form, e.g. 'BCAB'."""
        return canonical(self.mnumber)

    @property
    def quadruple(self):
        """(v(l), lambda_s, lambda_p, v(r)) -- the signature."""
        return self.mnumber.quadruple

    def substitute(self, choice):
        """
        Replace atoms by K / P and simplify; choice(label) returns 'K', 'P'
        or None (the atom stays).

        >>> print(MExpr('C p (B s (B p A))').substitute({'A1': 'P', 'C1': 'P'}.get))
        B1 s B2
        """
        return MExpr(_fold(self.expr, lambda w: choice(w) or w, _join))

    def corners(self):
        """The four corners of Definition 4.2 (Z(0) = delta_0, Z(inf) = delta_inf)."""
        at = lambda a, c: str(self.substitute(lambda w: {'A': a, 'C': c}.get(w[0])))
        return {'delta_0': at('P', 'K'), 'delta_inf': at('K', 'P'),
                'delta_K': at('K', 'K'), 'delta_P': at('P', 'P')}

    def faults(self):
        """Every single hard fault: element shorted (K) and opened (P)."""
        v0, ls0, lp0, r0 = self.quadruple
        out = []
        for label in self.labels:
            for kind in 'KP':
                res = self.substitute({label: kind}.get)
                if res.ideal:
                    out.append(Fault(label, kind, res.ideal, res.ideal, None, 'port', False))
                    continue
                v, ls, lp, r = q = res.quadruple
                tests = [t for t, hit in (('DC', v != v0), ('HF', r != r0), ('sweep', ls + lp != ls0 + lp0)) if hit]
                out.append(Fault(label, kind, str(res), res.word, q, '+'.join(tests) or 'none',
                                 q == self.quadruple))
        return out

    def _node_of(self):
        """label -> (node index, 's' | 'p') of the node of the expression that contains the atom."""
        out = {}
        for i, level in enumerate(_levels(self.expr)):
            words = level.split()
            op = 'p' if 'p' in words else 's' if 's' in words else None
            for w in words[::2]:
                if _is_atom(w):
                    out[w] = (i, op)
        return out

    def ambiguity_groups(self):
        """{'structural': [[...], ...], 'absolute': [[...], ...], 'signature': [[...], ...]}
        over the faults, as in the paper.  Structural: the same residual network after
        merging of same-type neighbours.  Absolute: the sub-groups of a structural group
        consisting of all opens of one series node or all shorts of one parallel node
        (identical residual networks, element values included).  Signature: the same
        quadruple.  A collapsed port (K or P) counts as a network and as a signature of
        its own, and 'nominal' is listed in the signature group of the structurally
        undetectable faults."""
        node = self._node_of()
        shape, signature = {}, {self.quadruple: ['nominal']}
        for f in self.faults():
            key = f.residual if f.quadruple is None else _fold(f.residual, lambda w: w[0], _shape)
            shape.setdefault(key, []).append(f)
            signature.setdefault(f.quadruple or f.residual, []).append(f'{f.element} {f.kind}')
        structural, absolute = [], []
        for g in shape.values():
            if len(g) < 2:
                continue
            structural.append([f'{f.element} {f.kind}' for f in g])
            sub = {}
            for f in g:
                level, op = node[f.element]
                if (f.kind == 'P' and op == 's') or (f.kind == 'K' and op == 'p'):
                    sub.setdefault((level, f.kind), []).append(f'{f.element} {f.kind}')
            absolute.extend(v for v in sub.values() if len(v) >= 2)
        return {'structural': structural, 'absolute': absolute,
                'signature': [g for g in signature.values() if len(g) > 1]}

    def dictionary(self, file=sys.stdout):
        """Print the whole structural fault dictionary."""
        p = lambda *a: print(*a, file=file)
        v, ls, lp, r = self.quadruple
        others = sorted(word_class(self.word) - {self.word})
        p(f'network     {self}')
        p(f'm-word      {self.mnumber.l2U}  ->  {self.word}'
          + (f'  (= {", ".join(others)} by the reversal rule)' if others else ''))
        p(f'signature   Phi = {self.quadruple}   Z(s) = s^{v} P/Q, deg P = {lp}, deg Q = {ls}, lambda = {ls + lp}')
        p('corners     ' + '   '.join(f'{k} = {val}' for k, val in self.corners().items()))
        p('')
        p(f"{'fault':<9}{'faulty network':<28}{'m-word':<8}{'Phi':<16}{'detected by':<12}")
        for f in self.faults():
            note = 'STRUCTURALLY UNDETECTABLE' if f.undetectable else ''
            p(f"{f.element + ' ' + f.kind:<9}{f.residual:<28}{f.word:<8}{str(f.quadruple or '--'):<16}{f.test:<12}{note}")
        g = self.ambiguity_groups()
        p('')
        p(f"structural ambiguity groups (same topology after merging): {g['structural'] or 'none'}")
        p(f"absolute ambiguity groups (identical residual networks):    {g['absolute'] or 'none'}")
        p(f"signature ambiguity groups (same quadruple):                {g['signature'] or 'none'}")


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] != '--test':
        MExpr(sys.argv[1]).dictionary()
    else:
        import doctest
        print(doctest.testmod())
