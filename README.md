# mNumbers fault dictionary

Structural fault dictionaries of series-parallel RLC one-ports, computed from the schematic alone
with Šare's m-theory: no element values, no simulation, no symbolic network functions.

This repository accompanies the paper

> *From the Schematic to the Fault Dictionary: Structural Signatures of Hard Faults in RLC One-Ports*

and contains the code that produces its examples and benchmark, together with the supplementary
document `ladenheim_dictionaries.pdf` (one page per network for the 99 series-parallel networks
of Ladenheim's catalogue, 798 single hard faults).

## What it computes

A series-parallel one-port is written as an **m-expression**, a text such as `B s (C p (B s A))`:
the atoms `A`, `B`, `C` stand for a capacitor, a resistor and an inductor, `s` for the series and
`p` for the parallel connection. From this expression `MExpr` derives

| quantity | meaning |
|---|---|
| `word` | the m-word of the network in canonical form, e.g. `BCAB`; its first and last letter are the Bode asymptotes at $s\to0$ and $s\to\infty$, the letters in between count the turns of the magnitude plot |
| `quadruple` | the **signature** $\Phi=(v(l),\lambda_s,\lambda_p,v(r))$: the impedance is $Z(s)=s^{v(l)}P(s)/Q(s)$ with $\deg P=\lambda_p$ and $\deg Q=\lambda_s$ for every choice of element values |
| `corners()` | the resistive networks $\delta_0=Z(0)$, $\delta_\infty=Z(\infty)$ and the floor and ceiling $\delta_K$, $\delta_P$, obtained by shorting or opening all reactive elements |
| `faults()` | every single hard fault (each element shorted, `K`, and opened, `P`): the residual network, its word and signature, and the simplest test that detects the fault (`DC` if $v(l)$ changes, `HF` if $v(r)$ changes, `sweep` if the corner count changes, `port` if the port itself is shorted or opened) |
| `ambiguity_groups()` | faults that a port measurement cannot tell apart, in the three kinds of the paper: **structural** groups (same residual network after merging of same-type neighbours), **absolute** groups (sub-groups of a structural group formed by all opens of one series node or all shorts of one parallel node, hence identical residual networks, values included) and **signature** groups (same quadruple) |

A fault whose signature equals the nominal one is **structurally undetectable** (field `undetectable`
in `Fault`): none of the three signature tests detects it for any choice of element values. Such a
fault may still be detectable through value-dependent features of the response, such as a change
of level or of damping. In the benchmark every such fault is a resistor fault.

## Installation

Python 3.10 or newer. The only hard dependency is the package
[`mnumbers`](https://github.com/JurajBenic/mNumbers), which implements the m-word arithmetic;
`schemdraw` and `pandas` are needed for the notebook only.

```bash
git clone <this repository>
cd mnumbers-fault-dictionary
pip install -r requirements.txt
```

## Quick start

```python
from mexpr import MExpr

e = MExpr('B s (C p (B s A))')          # Ladenheim network #63: R2 + (L1 || (R1 + C1))
print(e)                                 # B1 s (C1 p (B2 s A1))   (atoms are labelled)
print(e.word, e.quadruple)               # BCAB (0, 2, 2, 0)
print(e.corners())                       # {'delta_0': 'B1', 'delta_inf': 'B1 s B2', 'delta_K': 'B1', 'delta_P': 'P'}

for f in e.faults():
    print(f.element, f.kind, f.residual, f.word, f.quadruple, f.test, f.undetectable)

print(e.ambiguity_groups())
```

The whole dictionary of a network, formatted as in the paper, is printed from the command line:

```bash
python -m mexpr.core "B s (C p (B s A))"
```

```
network     B1 s (C1 p (B2 s A1))
m-word      BBcbCCacABbAaBBB  ->  BCAB  (= BACB by the reversal rule)
signature   Phi = (0, 2, 2, 0)   Z(s) = s^0 P/Q, deg P = 2, deg Q = 2, lambda = 4
corners     delta_0 = B1   delta_inf = B1 s B2   delta_K = B1   delta_P = P

fault    faulty network              m-word  Phi             detected by
B1 K     C1 p (B2 s A1)              CAB     (1, 2, 1, 0)    DC+sweep
B1 P     P                           P       --              port
C1 K     B1                          B       (0, 0, 0, 0)    sweep
C1 P     B1 s B2 s A1                AB      (-1, 0, 1, 0)   DC+sweep
B2 K     B1 s (C1 p A1)              BCAB    (0, 2, 2, 0)    none        STRUCTURALLY UNDETECTABLE
B2 P     B1 s C1                     BC      (0, 0, 1, 1)    HF+sweep
A1 K     B1 s (C1 p B2)              BCB     (0, 1, 1, 0)    sweep
A1 P     B1 s C1                     BC      (0, 0, 1, 1)    HF+sweep

structural ambiguity groups (same topology after merging): [['B2 P', 'A1 P']]
absolute ambiguity groups (identical residual networks):    [['B2 P', 'A1 P']]
signature ambiguity groups (same quadruple):                [['nominal', 'B2 K'], ['B2 P', 'A1 P']]
```

Any residual network can be built by hand with `substitute`, which replaces chosen atoms by `K`
or `P` and simplifies with the absorption rules `x s K = x`, `x s P = P`, `x p P = x`, `x p K = K`:

```python
print(e.substitute({'B2': 'K'}.get))     # B1 s (C1 p A1)
```

The doctests of the module are run with `python -m mexpr.core --test`.

## Writing m-expressions

* Atoms: `A` capacitor, `B` resistor, `C` inductor. Labels `A1`, `B2`, ... are assigned
  automatically in reading order and are used in the fault tables.
* `x s y` is the series and `x p y` the parallel connection. Several operands of the same kind
  may be chained without parentheses (`B s C s A`), mixed `s` and `p` on one level need them
  (`B s (A p C)`, not `B s A p C`).
* Only series-parallel networks are covered. Bridge networks cannot be written as m-expressions.

## The notebook

`example.ipynb` walks through the catalogue: it loads `data/ladenheim_networks.csv` (network
number, realisation with circuit elements, m-expression), draws the schematic of any network with
`schematic.draw` (schemdraw; the grey labels under the elements are the labels of the fault
tables), prints its fault dictionary as a `pandas` table, draws residual networks, and shows how
to analyse a network of your own.

## The supplementary document

`ladenheim_dictionaries.pdf` contains, for each of the 99 series-parallel networks of
Ladenheim's catalogue, the schematic, the m-word and signature, the corners, the table of all
single hard faults and the ambiguity groups, preceded by the summary tables of the paper: of the
798 faults, 744 (93.2 %) are structurally detectable and 54 (6.8 %) are structurally undetectable;
410 faults (51.4 %) lie in 146 structural ambiguity groups, of which 272 faults (34.1 %) form 134
absolute groups with identical residual networks.

## Repository layout

```
mexpr/core.py              MExpr, Fault, canonical, word_class  (the dictionary engine)
schematic.py               schemdraw drawing of a network from its m-expression
example.ipynb              the catalogue, one network in detail, your own network
data/ladenheim_networks.csv  the 99 series-parallel networks of Ladenheim's catalogue
ladenheim_dictionaries.pdf   supplementary document, one page per network
requirements.txt
```

## References

* M. Šare, *Jorbologija* (Jorbology), Element, Zagreb, 2000 (in Croatian).
* M. Orešković, I. Kuzmanović Ivičić, J. Benić, M. Essert, A perspective on the algebra,
  topology, and logic of electrical networks, *Int. J. Circuit Theory Appl.* 54 (2026) 4588–4604,
  doi:10.1002/cta.70359.
* E. L. Ladenheim, *A Synthesis of Biquadratic Impedances*, M.Sc. thesis, Polytechnic Institute of
  Brooklyn, 1948; A. Morelli, M. C. Smith, *Passive Network Synthesis: An Approach to
  Classification*, SIAM, 2019.
* The `mnumbers` package: https://github.com/JurajBenic/mNumbers

## License

MIT, see `LICENSE`.
