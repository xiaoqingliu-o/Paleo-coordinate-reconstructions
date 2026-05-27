# Paleo-Coordinate Reconstructions

A Python tool for reconstructing **paleo-latitude** and **paleo-longitude** from
modern site coordinates using [pyGPlates](https://www.gplates.org/docs/pygplates/)
plate-motion models.

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data Format](#input-data-format)
- [Usage](#usage)
- [Plate Models](#plate-models)
- [Output](#output)
- [References](#references)

---

## Repository Structure

```
.
├── extract_paleo_coordinates.py          # Main processing script
├── test_dataset/
│   └── Template_Site982.csv             # Example input (ODP Site 982)
└── Plate_model_rotation_and_StaticPolygons_files/
    ├── muller2008/                       # Müller et al. (2008) plate model
    │   ├── Rotations/
    │   └── StaticPolygons/
    └── zahirovic2022/                    # Zahirovic et al. (2022) plate model
        ├── Rotations/
        └── StaticPolygons/
```

Plate-model files for `muller2008` and `zahirovic2022` are included in this
repository and ready to use. If you need additional plate models, they can be
downloaded with [**plate-model-manager**](https://gplates.github.io/plate-model-manager/latest/index.html)
(see [Installation](#installation)).

---

## Requirements

| Package            | Purpose                           | Included with pygplates? |
|--------------------|-----------------------------------|--------------------------|
| `pygplates`        | Plate-tectonic reconstruction     | —                        |
| `numpy`            | Numerical operations              | Yes (conda dependency)   |
| `pandas`           | CSV I/O and tabular data handling | **No — install separately** |

> **Note:** When installing `pygplates` via conda, `numpy` is pulled in
> automatically as a required dependency. `pandas` must be installed
> explicitly (see [Installation](#installation)).

---

## Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/<your-username>/Paleo-coordinate-reconstructions.git
   cd Paleo-coordinate-reconstructions
   ```

2. **Create and activate a conda environment with pyGPlates**

   Following the [official pyGPlates documentation](https://www.gplates.org/docs/pygplates/pygplates_getting_started),
   create a new conda environment. `numpy` is installed automatically as
   a dependency of `pygplates`.

   ```bash
   conda create -n pygplates -c conda-forge python=3.13 pygplates
   conda activate pygplates
   ```

3. **Install `pandas`**

   ```bash
   conda install -c conda-forge pandas
   ```

4. **(Optional) Download additional plate models**

   The `muller2008` and `zahirovic2022` model files are already included in
   this repository. If you need other plate models, first install
   `plate-model-manager` into the same conda environment:

   ```bash
   conda install conda-forge::plate-model-manager
   ```
   To download the rotation file and the StaticPolygons layer from a plate model, follow the [plate-model-manager instructions] (https://gplates.github.io/plate-model-manager/latest/basic_usages.html#download-rotation-files)
  
   For available model names see the
   [plate-model-manager documentation](https://gplates.github.io/plate-model-manager/latest/plate_models.html).

---

## Input Data Format

Place your site-data CSV files in `test_dataset/`. The file **must** contain
the following four columns (in any order, but with these exact headers):

| Column      | Description                                       |
|-------------|---------------------------------------------------|
| `Site`      | Site identifier (string or integer)               |
| `Longitude` | Modern longitude in decimal degrees (WGS-84)      |
| `Latitude`  | Modern latitude in decimal degrees (WGS-84)       |
| `Age (Ma)`  | Reconstruction age in millions of years (Ma)      |

**Example** (`Template_Site982.csv`):

| Site | Longitude  | Latitude  | Age (Ma) |
|------|------------|-----------|----------|
| 982  | -15.854183 | 57.512667 | 0.0001   |
| 982  | -15.854183 | 57.512667 | 0.5      |
| 982  | -15.854183 | 57.512667 | 1.0      |
| 982  | -15.854183 | 57.512667 | 1.5      |

> **Note:** Multiple sites can be included in the same file. Each row is
> treated independently, so a single site may appear at many different ages.

---

## Usage

1. Place your input CSV file in `test_dataset/`.
2. Open `extract_paleo_coordinates.py` and, in the `if __name__ == "__main__":`
   block, set `INPUT_CSV` to your file and uncomment the plate-model block
   that matches your application (see [Plate Models](#plate-models)).
3. Run the script from within your `pygplates` conda environment:

   ```bash
   conda activate pygplates
   python extract_paleo_coordinates.py
   ```

The function can also be imported and called programmatically:

```python
from extract_paleo_coordinates import reconstruct_paleo_coordinates

reconstruct_paleo_coordinates(
    input_path      = "test_dataset/Template_Site982.csv",
    rotation_path   = "Plate_model_rotation_and_StaticPolygons_files/zahirovic2022/Rotations/CombinedRotations.rot",
    polygon_path    = "Plate_model_rotation_and_StaticPolygons_files/zahirovic2022/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.shp",
    anchor_plate_id = 701701,   # paleomagnetic reference frame
)
```

---

## Plate Models

Two plate models are included in this repository. Choose the one that matches
your scientific application.

### Zahirovic et al. (2022) *(default)*

Set the anchor plate ID to `701701` to use paleomagnetic reference frame. The paleomagnetic frame is tied to Earth's spin axis and is preferred for recovering accurate paleolatitudes and paleolongitudes. As the DeepMIP community 
plans future simulations using paleogeography based on the Zahirovic et al. (2022) model, 
this model should be used when reconstructing paleo-coordinates to ensure consistency when comparing proxy data with those simulations.

### Müller et al. (2008)

Use `anchor_plate_id=0`. This plate model is recommended for reconstructing paleo-coordinates, when performing point-wise model-data comparisons with published **MioMIP1** and **DeepMIP-Eocene** simulations, as most of those simulations
use paleogeography based on the Müller et al. (2008) model.

---

## Output

The script writes a new CSV file to the same directory as the input, appending
`_with_paleo_coordinates` to the filename:

```
test_dataset/Template_Site982_with_paleo_coordinates.csv
```

Two columns are appended to the original data:

| Column           | Description                                        |
|------------------|----------------------------------------------------|
| `paleolatitude`  | Reconstructed latitude (decimal degrees, 3 d.p.)  |
| `paleolongitude` | Reconstructed longitude (decimal degrees, 3 d.p.) |

Rows that fall outside the plate-polygon coverage are assigned `NaN` and a
warning is printed to the console.

---

## References

If you use results generated by this tool in a publication, please cite the
relevant plate model:

- **Zahirovic et al. (2022)**
  Subduction kinematics and carbonate platform interactions.
  *Geoscience Data Journal.*
  <https://doi.org/10.1002/gdj3.146>


- **Merdith et al. (2021)** *(for the paleomagnetic reference frame)*
  Extending full-plate tectonic models into deep time: Linking the
  Neoproterozoic and the Phanerozoic.
  *Earth-Science Reviews*, 214, 103477.
  <https://doi.org/10.1016/j.earscirev.2020.103477>
