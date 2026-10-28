<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Contributing

Bug reports, documentation fixes and pull requests are welcome through GitHub issues and pull requests.

## Before you open a pull request

1. Build and run the tests:

   ```sh
   cmake -B build -S .
   cmake --build build -j
   cmake --build build --target test_amd test_inflow_db test_inflow_gen test_inflow_planes test_sem_dfsem
   ctest --test-dir build --output-on-failure
   ```

   Regression cases need `gmsh` on `PATH`. A change to numerics or physics should come with a test that
   exercises the invariant it affects.
2. Follow the style rules in [docs/development.md](docs/development.md#style): `implicit none`, a public
   interface block at the top of each `mod_<area>.f90`, short comments, `error stop` messages under 80
   characters, lines under 133 characters.
3. Kernels that have an OpenACC port keep both the AoS and the SoA code paths in step.
4. New behaviour goes behind an input flag whose default keeps existing input files working, and is
   documented in [docs/input-reference.md](docs/input-reference.md).
5. Keep commits focused, with a short one-line message such as `numerics: Fix ...` or `docs: Add ...`.

## Reporting a problem

Use the issue templates. Include the input file, the mesh type and size, the rank count, and the compiler
and MPI versions.

## License

By contributing you agree that your contribution is licensed under AGPL-3.0-only, the license of the project.
