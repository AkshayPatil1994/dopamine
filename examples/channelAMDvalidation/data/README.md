# DNS reference data

`postprocess.py --dns` expects the Re_tau = 2000 channel DNS of Hoyas & Jimenez (2006) as a two-column
file, `data/hoyas_jimenez_retau2000.csv`:

```
# y/delta, U+
0.0, 0.0
...
```

The data is not shipped with this repository (it is not available offline here); obtain the mean profile from
the Hoyas-Jimenez / Jimenez-group channel database and convert `y/delta` and `U+` to this format. Without the
file, `postprocess.py` still writes the numerical profiles but makes no comparison.

For a Re_tau = 395 sanity check of the machinery, the Moser-Kim-Mansour means ship with OpenFOAM
(`.../turbulentInflow/oneCellThickPlaneChannel/resources/dataset/chan395.means`).
