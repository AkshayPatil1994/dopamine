---
title: Pre-processing
nav_order: 5
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Pre-processing

Everything that happens between a geometry and a running solver: building the mesh, naming boundary
patches, and optionally decomposing the mesh for large runs. Case parameters are in the
[Input Reference](input-reference.md).

## Boundary patch names

The names in `patch_name(i)` of `&dopamine_bc` must match the boundary patches of the mesh exactly:

- **Gmsh:** the names of the `Physical Surface` groups (`$PhysicalNames`). A boundary face that belongs to no
  physical surface stops the run.
- **OpenFOAM:** the patch names in `constant/polyMesh/boundary`.
- **Periodic pairs:** periodic patches must come in pairs of matching faces; the solver stops if a patch has no partner.

A mesh patch that is not listed in `&dopamine_bc` takes the default type with a warning; a name listed in
`&dopamine_bc` that is absent from the mesh also only prints a warning, so check the startup log.

## Choosing a Mesh Path

`mesh_type` selects between two quite different startup paths. The choice is
about where the global mesh has to fit, not about accuracy — both produce the
same partitions.

| | Direct read (`gmsh`, `openfoam`) | Preprocessed |
|---|---|---|
| Preparation | none | one offline `dopamine-decompose` run |
| Global mesh held by | rank 0, entirely | nobody at solver time |
| Startup | re-partitioned every run | each rank reads its own file |
| Rank count | change freely | fixed at decomposition time |
| Practical range | up to ~10M cells | above that |

**Direct read** is the default and the right choice while the mesh still fits
on one rank. Rank 0 reads the whole mesh, partitions it and scatters; nothing to
prepare and you can vary `-np` between runs.

```fortran
&dopamine_domain
  mesh_file = 'channel.msh'
  mesh_type = 'gmsh'
/
```

**Preprocessed** removes the rank-0 ceiling. Decompose once, then every rank
reads only its own `processorN/mesh.bin`:

```bash
mpirun -np 1 ./build/dopamine-decompose channel.input --nranks 128 \
  --output-dir ./channel_128 --threads 16
```

```fortran
&dopamine_domain
  mesh_file = './channel_128'
  mesh_type = 'preprocessed'
/
```

Run it with exactly the rank count it was decomposed for
(`mpirun -np 128`); changing `-np` means decomposing again.

Switch to preprocessed when rank 0 runs out of memory during startup, when the
mesh is beyond roughly 10M cells, or when you will launch the same mesh many
times at a fixed rank count and want the startup cost paid once.

### Decomposer options

| Flag | Meaning |
|---|---|
| `--nranks N` | partitions to write (required) |
| `--output-dir DIR` | destination for `processorN/mesh.bin` |
| `--threads N` | writer threads; defaults to all cores |
| `--max-memory-gb G` | cap peak memory by staging face geometry to disk |
| `--chunk-size M` | face scan chunk, default 2M |
| `--weight-by-volume` | balance total volume per rank instead of cell count |

Partitions carry an equal number of cells by default, because solver cost scales
with cell count. `--weight-by-volume` restores the older behaviour of equalising
total volume, which on a wall-graded channel or ABL mesh leaves the near-wall
ranks holding several times more cells than average — roughly 3.5x at 1024
ranks on a tanh-graded mesh. Leave it off unless you specifically want it.

The decomposer itself runs on a single core (`-np 1`); `--threads` parallelises
the write phase only. It prints the resulting balance, which is the number to
check:

```
Cells/rank  min=130918  max=131204  imbalance=1.001
Load balance: equal cell count
```

Note that the direct-read path has always balanced by cell count, so this
setting only affects meshes that go through the decomposer.

## Supported Mesh Formats

- **Gmsh v4 (.msh):** tet4, hex8, prism6, pyr5, tet10, hex20. Use
  `$PhysicalNames` for boundaries. Generate with
  `gmsh -3 -format msh4 geometry.geo`. Every boundary surface must belong to a
  physical group: the reader stops with an error if a cell face has neither a
  neighbour nor a labelled surface element.
- **OpenFOAM polyMesh:** Reads `constant/polyMesh/{points,faces,owner,
  neighbour,boundary}`. All OpenFOAM-supported cell types.

## Next steps

- Run on the decomposed mesh: [Usage](usage.md#run).
- All flags of the decomposition and checkpoint tools: [Tools](tools.md).
