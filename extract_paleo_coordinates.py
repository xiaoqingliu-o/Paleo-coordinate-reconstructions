"""
extract_paleo_coordinates.py
============================
Reconstruct paleo-latitude and paleo-longitude for site locations using
pyGPlates plate-motion models.

Supported plate models
----------------------
Müller et al. (2008)
    Global ocean-crust age and spreading-rate model.
    anchor_plate_id = 0  (optimized mantle reference frame)
    Recommended for reconstructing paleo-coordinates when performing point-wise model-data 
    comparisons with published MioMIP1 and DeepMIP-Eocene simulations, as most of those 
    simulations use paleogeography based on the Müller et al. (2008) model.

Zahirovic et al. (2022)
    anchor_plate_id = 701701: paleomagnetic reference frame
    Recommended for reconstructing paleo-coordinates when performing point-wise model-data 
    comparisons with future DeepMIP simulations, as those simulations are planned to use 
    paleogeography based on the Zahirovic et al. (2022) model.

Plate model files
-----------------
Rotation (.rot) and static-polygon (.shp / .gpml) files for
``muller2008`` and ``zahirovic2022`` are included in this repository
under ``Plate_model_rotation_and_StaticPolygons_files/``.
Additional plate models can be downloaded with plate-model-manager:
    https://gplates.github.io/plate-model-manager/latest/index.html

Input CSV columns (order matters)
----------------------------------
  Site        : site identifier (string or integer)
  Longitude   : modern longitude (decimal degrees)
  Latitude    : modern latitude  (decimal degrees)
  Age (Ma)    : reconstruction age (millions of years before present)

Output
------
A CSV file written alongside the input file, with the suffix
``_with_paleo_coordinates.csv``. Two columns are appended:
  paleolatitude   : reconstructed latitude  (decimal degrees, rounded to 3 d.p.)
  paleolongitude  : reconstructed longitude (decimal degrees, rounded to 3 d.p.)

References
----------
Müller, R. D., et al. (2008). Age, spreading rates, and spreading
  asymmetry of the world's ocean crust. Geochem. Geophys. Geosyst., 9,
  Q04006. https://doi.org/10.1029/2007GC001743

Zahirovic, S., et al. (2022). Subduction kinematics and carbonate
  platform interactions. Geosci. Data J., 9, 1–21.
  https://doi.org/10.1002/gdj3.146
"""

import sys
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pygplates

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Required CSV columns
# ---------------------------------------------------------------------------
REQUIRED_COLUMNS = {"Site", "Longitude", "Latitude", "Age (Ma)"}


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def reconstruct_paleo_coordinates(
    input_path,
    rotation_path,
    polygon_path,
    anchor_plate_id=701701,
):
    """
    Reconstruct paleo-coordinates for all rows in *input_path* and write
    results to a new CSV file.

    Parameters
    ----------
    input_path : str or Path
        Path to the input CSV file.
    rotation_path : str or Path
        Path to the plate-model rotation file (.rot).
    polygon_path : str or Path
        Path to the static-polygon file (.gpml, .gpmlz, or .shp).
    anchor_plate_id : int, optional
        Anchor plate ID that defines the absolute reference frame.
        Default is 701701 (paleomagnetic frame, Merdith et al., 2021).

    Returns
    -------
    Path
        Path to the output CSV file.

    Raises
    ------
    FileNotFoundError
        If any of the three input paths do not exist.
    ValueError
        If the input CSV is missing required columns.
    """
    input_path    = Path(input_path)
    rotation_path = Path(rotation_path)
    polygon_path  = Path(polygon_path)

    # --- Validate paths ---------------------------------------------------
    for path in (input_path, rotation_path, polygon_path):
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    # --- Load and validate the input table --------------------------------
    log.info("Reading input file: %s", input_path)
    df = pd.read_csv(input_path)

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Input CSV is missing required column(s): {sorted(missing)}\n"
            f"Expected: {sorted(REQUIRED_COLUMNS)}"
        )

    n_rows = len(df)
    log.info("  %d rows loaded.", n_rows)

    # --- Build rotation model ---------------------------------------------
    log.info("Loading rotation model (anchor plate ID = %d) ...", anchor_plate_id)
    rotation_model = pygplates.RotationModel(
        str(rotation_path),
        default_anchor_plate_id=anchor_plate_id,
    )

    # --- Create point features (one per row) ------------------------------
    log.info("Creating point features ...")
    point_features = []
    for idx, row in df.iterrows():
        point   = pygplates.PointOnSphere(float(row["Latitude"]), float(row["Longitude"]))
        feature = pygplates.Feature()
        feature.set_geometry(point)
        feature.set_name(str(idx))   # row index used as a unique key
        point_features.append(feature)

    # --- Assign plate IDs via static polygons -----------------------------
    log.info("Partitioning points into tectonic plates ...")
    partitioned_features = pygplates.partition_into_plates(
        str(polygon_path),
        rotation_model,
        point_features,
    )

    # --- Reconstruct, grouped by age to minimise pygplates calls ----------
    reconstructed_lats = [np.nan] * n_rows
    reconstructed_lons = [np.nan] * n_rows

    ages = df["Age (Ma)"].unique()
    log.info("Reconstructing across %d unique age(s) ...", len(ages))

    for age, group in df.groupby("Age (Ma)"):
        indices  = set(group.index.tolist())
        features = [f for f in partitioned_features if int(f.get_name()) in indices]

        reconstructed_geometries = []
        pygplates.reconstruct(
            features,
            rotation_model,
            reconstructed_geometries,
            age,
            anchor_plate_id=anchor_plate_id,
        )

        n_reconstructed = 0
        for geom in reconstructed_geometries:
            idx = int(geom.get_feature().get_name())
            lat, lon = geom.get_reconstructed_geometry().to_lat_lon()
            reconstructed_lats[idx] = round(lat, 3)
            reconstructed_lons[idx] = round(lon, 3)
            n_reconstructed += 1

        n_expected = len(group)
        if n_reconstructed < n_expected:
            log.warning(
                "  Age %.4f Ma: %d/%d points reconstructed "
                "(missing points may lie outside polygon coverage).",
                age, n_reconstructed, n_expected,
            )

    # --- Append results and write output ----------------------------------
    df["paleolatitude"]  = reconstructed_lats
    df["paleolongitude"] = reconstructed_lons

    n_missing = df["paleolatitude"].isna().sum()
    if n_missing:
        log.warning(
            "%d row(s) could not be reconstructed and will contain NaN.", n_missing
        )

    output_path = input_path.parent / f"{input_path.stem}_with_paleo_coordinates.csv"
    df.to_csv(output_path, index=False)
    log.info("Results written to: %s", output_path)

    return output_path


# ---------------------------------------------------------------------------
# Command-line entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    BASE_DIR = Path(__file__).parent

    # Input data
    INPUT_CSV = BASE_DIR / "test_dataset" / "Template_Site982.csv"

    # --- Choose one plate model by uncommenting the relevant block --------

    # Müller et al. (2008)  —  use anchor_plate_id=0
    # ROT_FILE  = (BASE_DIR / "Plate_model_rotation_and_StaticPolygons_files"
    #              / "muller2008" / "Rotations"
    #              / "Global_EarthByte_GPlates_Rotation_20100927.rot")
    # POLY_FILE = (BASE_DIR / "Plate_model_rotation_and_StaticPolygons_files"
    #              / "muller2008" / "StaticPolygons"
    #              / "Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_20100927.gpml")
    # ANCHOR_ID = 0

    # Zahirovic et al. (2022) — paleomagnetic frame (default)
    ROT_FILE  = (BASE_DIR / "Plate_model_rotation_and_StaticPolygons_files"
                 / "zahirovic2022" / "Rotations" / "CombinedRotations.rot")
    POLY_FILE = (BASE_DIR / "Plate_model_rotation_and_StaticPolygons_files"
                 / "zahirovic2022" / "StaticPolygons"
                 / "Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.shp")
    ANCHOR_ID = 701701

    try:
        reconstruct_paleo_coordinates(INPUT_CSV, ROT_FILE, POLY_FILE, anchor_plate_id=ANCHOR_ID)
    except (FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        sys.exit(1)
