# rootphi

`rootphi` is a symbolic arithmetic library for representing and manipulating numbers of the form:

$$
\frac{a + b \sqrt{\varphi} + c \varphi + d \varphi \sqrt{\varphi}}{e}
$$

Where $\varphi = \frac{1 + \sqrt{5}}{2}$ is the **golden ratio**, and $a, b, c, d, e \in \mathbb{Z}$.

---

## ✨ Features

* Exact arithmetic with irrational algebraic numbers ($\varphi$, $\sqrt{\varphi}$)
* Fully symbolic operations: `+`, `-`, `*`, `/`, power, negation
* String and decimal representations
* Integration with SymPy for symbolic math
* Precision-preserving `.to_decimal_string()` and `.digits()`
* Hashing and equality for use in sets and dicts

---

## 📦 Installation

```bash
pip install rootphi  # (once published)
```

Or clone and use directly:

```bash
git clone https://github.com/yourname/rootphi.git
cd rootphi
```

---

## 🧪 Example

```python
from rootphi import rootphi

phi = rootphi.phi()
x = rootphi([1, 0, 1, 0, 1])  # 1 + φ

print(x.to_string())  # "(1 + φ)"
print((x / phi).to_decimal_string())  # ~1.618...
```

---

## 🧠 Applications

* Symbolic geometry
* Quasicrystals and Penrose tilings
* Algebraic number theory
* Exact irrational vectors and golden-angle rotations
* Educational tools

---

## 🧰 Development

Run tests:

```bash
python -m unittest discover tests
```

---

## 📜 License

MIT License

---

## 🙌 Acknowledgments

* Inspired by algebraic number theory and symbolic computation
* Golden ratio $\varphi$ and its powers appear in many mathematical and physical systems

---

## 🔗 See Also

* [SymPy](https://www.sympy.org/)
* [NumPy](https://numpy.org/)
* [Golden Ratio - Wikipedia](https://en.wikipedia.org/wiki/Golden_ratio)

---

> 📣 Want to contribute? Fork this project, add new features (like `from_decimal`, LaTeX rendering, or vector support), and submit a pull request!
