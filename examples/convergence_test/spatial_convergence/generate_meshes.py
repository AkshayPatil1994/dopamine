#!/usr/bin/env python3
"""
Generate structured gmsh meshes for 3D Taylor-Green vortex convergence tests.
Creates periodic cubic domains [0, 2pi]^3 at multiple resolutions.
Based on taylor_green.geo from examples/taylor_green/
"""

import sys
import os
import subprocess
from pathlib import Path

def generate_geo_file(resolution, output_file):
    """Generate gmsh .geo file from taylor_green template."""

    L = 6.283185307179586  # 2pi

    geo_content = f"""// Taylor-Green vortex periodic box mesh (convergence test)
// Generates a structured hex mesh over [0, 2pi]^3 at {resolution}^3 resolution
// Usage:
//   gmsh -3 -format msh4 -o mesh_{resolution}x{resolution}.msh mesh_{resolution}x{resolution}.geo

L = {L};

// Resolution
n_points = {resolution};
Nx = n_points;
Ny = n_points;
Nz = n_points;

Point(1) = {{0, 0, 0, 1.0}};
Point(2) = {{L, 0, 0, 1.0}};
Point(3) = {{L, L, 0, 1.0}};
Point(4) = {{0, L, 0, 1.0}};
Point(5) = {{0, 0, L, 1.0}};
Point(6) = {{L, 0, L, 1.0}};
Point(7) = {{L, L, L, 1.0}};
Point(8) = {{0, L, L, 1.0}};

Line(1) = {{1, 2}};
Line(2) = {{2, 3}};
Line(3) = {{3, 4}};
Line(4) = {{4, 1}};
Line(5) = {{5, 6}};
Line(6) = {{6, 7}};
Line(7) = {{7, 8}};
Line(8) = {{8, 5}};
Line(9) = {{1, 5}};
Line(10) = {{2, 6}};
Line(11) = {{3, 7}};
Line(12) = {{4, 8}};

Line Loop(1) = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}}; // z=0
Line Loop(2) = {{5, 6, 7, 8}};
Plane Surface(2) = {{2}}; // z=L
Line Loop(3) = {{-1, 9, 5, -10}};
Plane Surface(3) = {{3}}; // y=0
Line Loop(4) = {{-3, 11, 7, -12}};
Plane Surface(4) = {{4}}; // y=L
Line Loop(5) = {{-4, 12, 8, -9}};
Plane Surface(5) = {{5}}; // x=0
Line Loop(6) = {{-2, 10, 6, -11}};
Plane Surface(6) = {{6}}; // x=L

Surface Loop(1) = {{1, 2, 3, 4, 5, 6}};
Volume(1) = {{1}};

// Structured mesh using Transfinite
Transfinite Line{{1, 3, 5, 7}} = Nx + 1;
Transfinite Line{{2, 4, 6, 8}} = Ny + 1;
Transfinite Line{{9, 10, 11, 12}} = Nz + 1;
Transfinite Surface{{1, 2, 3, 4, 5, 6}} Alternated;
Transfinite Volume{{1}};
Recombine Surface{{1, 2, 3, 4, 5, 6}};
Recombine Volume{{1}};

Physical Surface("periodic_x0", 1) = {{5}};
Physical Surface("periodic_xL", 2) = {{6}};
Physical Surface("periodic_y0", 3) = {{3}};
Physical Surface("periodic_yL", 4) = {{4}};
Physical Surface("periodic_z0", 5) = {{1}};
Physical Surface("periodic_zL", 6) = {{2}};
Physical Volume("fluid", 7) = {{1}};

Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
"""

    with open(output_file, 'w') as f:
        f.write(geo_content)

    print(f"  [ok] Generated {output_file}")


def run_gmsh(geo_file, msh_file):
    """Run gmsh to generate .msh file from .geo file."""
    try:
        result = subprocess.run(
            ['gmsh', '-3', '-format', 'msh4', '-o', msh_file, geo_file],
            capture_output=True,
            timeout=300
        )
        if result.returncode != 0:
            print(f"  [x] gmsh failed: {result.stderr.decode()}")
            return False
        print(f"  [ok] Generated {msh_file}")
        return True
    except FileNotFoundError:
        print("  [x] gmsh not found. Install with: apt-get install gmsh")
        return False
    except subprocess.TimeoutExpired:
        print(f"  [x] gmsh timeout (increase timeout for large meshes)")
        return False


def main():
    script_dir = Path(__file__).parent
    os.chdir(script_dir)

    resolutions = [16, 32, 64, 128]

    print("=" * 70)
    print("Generating 3D structured meshes for Taylor-Green vortex convergence")
    print("=" * 70)
    print("Template: Based on taylor_green.geo (transfinite structured hex mesh)")
    print("")

    all_success = True
    for res in resolutions:
        print(f"\nResolution: {res}^3 ({res**3:,d} cells)")

        geo_file = f"mesh_{res}x{res}.geo"
        msh_file = f"mesh_{res}x{res}.msh"

        # Generate .geo file
        generate_geo_file(res, geo_file)

        # Run gmsh
        if not run_gmsh(geo_file, msh_file):
            all_success = False
            print(f"  WARNING: Mesh generation failed for {res}^3")
            continue

        print(f"  Mesh size: {msh_file}")

    print("\n" + "=" * 70)
    if all_success:
        print("[ok] All meshes generated successfully")
        print("\nNext: Run ./run_convergence.sh to start convergence tests")
    else:
        print("[!] Some meshes failed to generate")
        print("  Install gmsh: apt-get install gmsh")
    print("=" * 70)

    return 0 if all_success else 1


if __name__ == '__main__':
    sys.exit(main())
