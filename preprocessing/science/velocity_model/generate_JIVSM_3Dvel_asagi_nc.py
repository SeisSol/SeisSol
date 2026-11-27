#!/usr/bin/env python3
"""Build an ASAGI NetCDF material model from the JIVSM layer model.
The model can be obtained via https://www.j-shis.bosai.go.jp/map/JSHIS2/download.html

Please cite JIVSM references: 
- Koketsu, K., Miyake, H., & Suzuki, H. (2012). Japan integrated velocity structure model version 1. In Proceedings of the 15th world conference on earthquake engineering (Vol. 1, p. 4). Lisbon.
https://www.iitk.ac.in/nicee/wcee/article/WCEE2012_1773.pdf
- Koketsu, K., Miyake, H., Fujiwara, H., & Hashimoto, T. (2008). Progress towards a Japan integrated velocity structure model and long-period ground motion hazard map. In Proceedings of the 14th World conference on earthquake engineering (pp. S10-038). Beijing: China Seismological Society.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Interpolate JIVSM onto a projected Cartesian grid and write an "
            "ASAGI-compatible NetCDF file plus a YAML reader snippet."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Example Tohoku run:\n"
            "  python %(prog)s --proj EPSG:6677 \\\n"
            "    --lon0 140 --lon1 140.02 --lat0 38 --lat1 38.02 \\\n"
            "    --dx 500 --dy 500 --dz 100 --min-vs 1700 --output Tohoku_JIVSM.nc"
        ),
    )
    parser.add_argument("--proj", required=True, help="Target CRS understood by PROJ (for example EPSG:6677 or a PROJ string)")
    parser.add_argument("--lon0", required=True, type=float, help="Minimum longitude in degrees")
    parser.add_argument("--lon1", required=True, type=float, help="Maximum longitude in degrees")
    parser.add_argument("--lat0", required=True, type=float, help="Minimum latitude in degrees")
    parser.add_argument("--lat1", required=True, type=float, help="Maximum latitude in degrees")
    parser.add_argument("--dx", required=True, type=float, help="Grid spacing in x, metres")
    parser.add_argument("--dy", required=True, type=float, help="Grid spacing in y, metres")
    parser.add_argument("--dz", required=True, type=float, help="Grid spacing in z, metres")
    parser.add_argument("-o", "--output", type=Path, default=Path("JIVSM_material.nc"), help="Output NetCDF path")
    parser.add_argument("--yaml", type=Path, help="Output YAML path (default: same stem as --output)")
    parser.add_argument("--jivsm", type=Path, default=here / "D-V4-STRUCT_DEEP-LYRD.csv", help="J-SHIS V4 layer-depth CSV")
    parser.add_argument("--velocity-csv", type=Path, default=here / "D-V4-STRUCT_DEEP-PYS.csv", help="J-SHIS V4 material-property CSV")
    parser.add_argument("--min-vs", type=float, default=0.0, help="Minimum shear velocity in m/s. Vp is raised where needed to keep elastic moduli stable.")
    return parser.parse_args()


def inclusive_axis(vmin: float, vmax: float, spacing: float) -> np.ndarray:
    """Return a regular axis which covers [vmin, vmax]."""
    start = math.floor(vmin / spacing) * spacing
    stop = math.ceil(vmax / spacing) * spacing
    count = int(round((stop - start) / spacing)) + 1
    return start + np.arange(count, dtype=np.float64) * spacing


def load_materials(csv_path: Path, min_vs: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(csv_path, delimiter=",", usecols=(0, 1, 2, 3), ndmin=2)
    table = table[np.argsort(table[:, 0])]
    if not np.isfinite(table).all() or not np.array_equal(table[:, 0], np.arange(1, len(table) + 1)):
        raise ValueError(f"Could not read Vp, Vs, and density from {csv_path}")
    _, vp, vs, rho = table.T  # PYS values are already in SI units.
    floored = vs < min_vs
    vs = np.maximum(vs, min_vs)
    vp[floored] = np.maximum(vp[floored], np.sqrt(2.0) * vs[floored])
    mu = rho * vs**2
    lam = rho * vp**2 - 2.0 * mu
    return rho, mu, lam


def mesh_code(row: np.ndarray, col: np.ndarray) -> np.ndarray:
    """Encode indices on the Japanese third-level (30 by 45 arcsec) mesh."""
    return (
        (row // 80) * 1_000_000 + (col // 80) * 10_000
        + (row % 80 // 10) * 1_000 + (col % 80 // 10) * 100
        + (row % 10) * 10 + col % 10
    )


def interpolate_interfaces(
    path: Path,
    layer_count: int,
    query_lon: np.ndarray,
    query_lat: np.ndarray,
) -> np.ndarray:
    # Mesh values are sampled at cell centers. Work in fractional mesh indices.
    row = query_lat * 120 - 0.5
    col = (query_lon - 100) * 80 - 0.5
    row0, col0 = np.floor(row).astype(int), np.floor(col).astype(int)
    fy, fx = row - row0, col - col0
    corners = [
        (mesh_code(row0, col0), (1 - fy) * (1 - fx)),
        (mesh_code(row0, col0 + 1), (1 - fy) * fx),
        (mesh_code(row0 + 1, col0), fy * (1 - fx)),
        (mesh_code(row0 + 1, col0 + 1), fy * fx),
    ]
    needed = {int(code) for codes, weights in corners for code in codes[weights > 0]}
    profiles = {}
    # Stream the CSV and parse depth columns only for cells needed by this grid.
    with path.open("rb") as source:
        for line in source:
            if not line.strip() or line.startswith(b"#"):
                continue
            code_text, _, values = line.partition(b",")
            code = int(code_text)
            if code in needed:
                depths = np.fromstring(values.decode("ascii"), sep=",")
                if depths.size != layer_count + 1 or not np.isfinite(depths).all():
                    raise ValueError(f"Invalid LYRD profile for mesh {code} in {path}")
                # D0 is the bottom of layer 0 and the top of PYS layer 1.
                profiles[code] = -np.maximum.accumulate(depths[1:])
    missing = needed.difference(profiles)
    if missing:
        raise ValueError(f"Requested grid extends outside LYRD coverage ({len(missing)} missing mesh cells, e.g. {min(missing):08d})")
    interfaces = np.zeros((*query_lon.shape, layer_count), dtype=np.float64)
    for codes, weights in corners:
        active = weights > 0
        interfaces[active] += np.array([profiles[code] for code in codes[active]]).reshape(-1, layer_count) * weights[active, None]
    return interfaces


def build_material_grid(
    z: np.ndarray,
    interfaces: np.ndarray,
    rho: np.ndarray,
    mu: np.ndarray,
    lam: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = (z.size, *interfaces.shape[:2])
    rho_grid = np.empty(shape, dtype=np.float32)
    mu_grid = np.empty(shape, dtype=np.float32)
    lam_grid = np.empty(shape, dtype=np.float32)

    # searchsorted requires increasing values; negate depths because JIVSM
    # interfaces normally run from the surface downward (decreasing z).
    for iy in range(shape[1]):
        for ix in range(shape[2]):
            boundaries = interfaces[iy, ix, :]
            layer = np.searchsorted(-boundaries, -z, side="right") - 1
            layer = np.clip(layer, 0, rho.size - 1)
            rho_grid[:, iy, ix] = rho[layer]
            mu_grid[:, iy, ix] = mu[layer]
            lam_grid[:, iy, ix] = lam[layer]
    return rho_grid, mu_grid, lam_grid


def write_netcdf(
    path: Path,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    arrays: tuple[np.ndarray, np.ndarray, np.ndarray],
    proj: str,
) -> None:
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise SystemExit("Missing dependency: install netCDF4 (for example: pip install netCDF4)") from exc

    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w", format="NETCDF4") as nc:
        for name, values in (("x", x), ("y", y), ("z", z)):
            nc.createDimension(name, values.size)
            coordinate = nc.createVariable(name, "f8", (name,))
            coordinate[:] = values
            coordinate.units = "m"
        dtype = np.dtype([("rho", "<f4"), ("mu", "<f4"), ("lambda", "<f4")])
        material_type = nc.createCompoundType(dtype, "material")
        data = nc.createVariable("data", material_type, ("z", "y", "x"))
        packed = np.empty(arrays[0].shape, dtype=dtype)
        for name, values in zip(("rho", "mu", "lambda"), arrays):
            packed[name] = values
        data[:] = packed
        nc.title = "JIVSM material model interpolated onto a Cartesian grid"
        nc.projection = proj


def write_yaml(path: Path, nc_path: Path, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"""# ASAGI material model generated with {Path(__file__).name}
!AxisAlignedCuboidalDomainFilter
limits:
  x: [{x[0]:.10g}, {x[-1]:.10g}]
  y: [{y[0]:.10g}, {y[-1]:.10g}]
  z: [{z[0]:.10g}, {z[-1]:.10g}]
components: !ASAGI
  file: {nc_path.resolve()}
  parameters: [rho, mu, lambda]
  var: data
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if not (args.lon0 < args.lon1 and args.lat0 < args.lat1):
        raise SystemExit("Require lon0 < lon1 and lat0 < lat1")
    if min(args.dx, args.dy, args.dz) <= 0:
        raise SystemExit("dx, dy, and dz must be positive")
    if not math.isfinite(args.min_vs) or args.min_vs < 0:
        raise SystemExit("--min-vs must be finite and nonnegative")
    try:
        from pyproj import Transformer
    except ImportError as exc:
        raise SystemExit("Missing dependency: install pyproj (for example: pip install pyproj)") from exc

    rho, mu, lam = load_materials(args.velocity_csv, args.min_vs)

    forward = Transformer.from_crs("EPSG:4326", args.proj, always_xy=True)
    inverse = Transformer.from_crs(args.proj, "EPSG:4612", always_xy=True)
    corner_lon = [args.lon0, args.lon0, args.lon1, args.lon1]
    corner_lat = [args.lat0, args.lat1, args.lat0, args.lat1]
    corner_x, corner_y = forward.transform(corner_lon, corner_lat)
    x = inclusive_axis(float(np.min(corner_x)), float(np.max(corner_x)), args.dx)
    y = inclusive_axis(float(np.min(corner_y)), float(np.max(corner_y)), args.dy)
    x_grid, y_grid = np.meshgrid(x, y)
    query_lon, query_lat = inverse.transform(x_grid, y_grid)
    try:
        interfaces = interpolate_interfaces(args.jivsm, rho.size, query_lon, query_lat)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    z = inclusive_axis(float(np.nanmin(interfaces)), float(np.nanmax(interfaces[:, :, 0])), args.dz)
    cell_count = x.size * y.size * z.size
    print(f"Grid: {x.size} x {y.size} x {z.size} = {cell_count:,} cells")
    arrays = build_material_grid(z, interfaces, rho, mu, lam)
    write_netcdf(args.output, x, y, z, arrays, args.proj)
    yaml_path = args.yaml or args.output.with_suffix(".yaml")
    write_yaml(yaml_path, args.output, x, y, z)
    print(f"Wrote {args.output}")
    print(f"Wrote {yaml_path}")


if __name__ == "__main__":
    main()




"""
Example of nested 3D velocity model with varying resolution. 
Each asagi nc file cannot have varying spacing.

!IdentityMap
  components:
    # fine grid
    - !AxisAlignedCuboidalDomainFilter
      limits:
        x: [-330000.0, 149500.0]
        y: [-460000.0, 239500.0]
        z: [-80000.0, -30000.0]
      components: !ASAGI
        file: /leonardo/home/userexternal/jwong001/work/Tohoku_FFM/EKI_FL103_Pf_0.90_plas/Tohoku_material/Tohoku_3D_velocity_model_fine_z80_25km.nc
        parameters: [rho, mu, lambda]
        var: data
    # fine grid
    - !AxisAlignedCuboidalDomainFilter
      limits:
        x: [-330000.0, 149500.0]
        y: [-460000.0, 239500.0]
        z: [-30000.0, 2000.0]
      components: !ASAGI
        file: /leonardo/home/userexternal/jwong001/work/Tohoku_FFM/EKI_FL103_Pf_0.90_plas/Tohoku_material/Tohoku_3D_velocity_model_very_fine_z30_2km_minVs1700.nc #Tohoku_3D_velocity_model_very_fine_z30_2km.nc
        parameters: [rho, mu, lambda]
        var: data
    # elsewhere
    - !AxisAlignedCuboidalDomainFilter
      limits:
        x: [-851638.0, 452361.0]
        y: [-740715.0, 514284.0]
        z: [-15000.0, 10000.0]
      components: !ASAGI
        file: /leonardo/home/userexternal/jwong001/work/Tohoku_FFM/EKI_FL103_Pf_0.90_plas/Tohoku_material/Tohoku_3D_velocity_model_large_shallow_15km.nc
        parameters: [rho, mu, lambda]
        var: data
    - !ASAGI
      file: /leonardo/home/userexternal/jwong001/work/Tohoku_FFM/EKI_FL103_Pf_0.90_plas/Tohoku_material/Tohoku_3D_velocity_model_large.nc
      parameters: [rho, mu, lambda]
      var: data

    - !ConstantMap
      map:
         rho: 3400.0
         mu: 71944003584.0
         lambda: 127148662784.0
"""
