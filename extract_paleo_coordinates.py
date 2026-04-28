"""
Input: 
  - input_path   (input csv file)
  - rotation_path  (rotation file)
  - polygon_path    (static polygon file)
  - anchor_plate_id 
  
Rotation model Zahirovic_et al. (2022) provides global plate reconstructions with multiple absolute reference frames. 
These are selected via the anchor plate ID, including an optimized mantle reference frame (0), a no-net-rotation frame (777777), and a paleomagnetic reference frame (701701).
This script uses 701701, which places the reconstruction in a paleomagnetic reference frame based on Merdith et al. (2021). This frame is tied to Earth’s spin axis 
and is therefore preferred for recovering accurate paleolatitudes, which is important for paleoclimate applications.


Output:
  - output paleo coordinates in a csv file

"""
import pygplates
import pandas as pd
import numpy as np
from pathlib import Path

def reconstruct_paleo_coordinates(input_path, rotation_path, polygon_path, anchor_plate_id=701701):
    input_path = Path(input_path)
    if not input_path.exists():
        print(f"Error: Could not find input file at {input_path}")
        return

    df = pd.read_csv(input_path)
    rotation_model = pygplates.RotationModel(str(rotation_path), default_anchor_plate_id=anchor_plate_id)
    
    point_features = []
    for index, row in df.iterrows():
        point = pygplates.PointOnSphere(float(row['Latitude']), float(row['Longitude']))
        feature = pygplates.Feature()
        feature.set_geometry(point)
        feature.set_name(str(index)) 
        point_features.append(feature)

    partitioned_features = pygplates.partition_into_plates(str(polygon_path), rotation_model, point_features)
    reconstructed_lats, reconstructed_lons = [np.nan] * len(df), [np.nan] * len(df)

    for age, group in df.groupby('Age (Ma)'):
        indices = group.index.tolist()
        features_to_reconstruct = [f for f in partitioned_features if int(f.get_name()) in indices]
        reconstructed_geometries = []
        pygplates.reconstruct(features_to_reconstruct, rotation_model, reconstructed_geometries, age, anchor_plate_id=anchor_plate_id)

        for geom in reconstructed_geometries:
            idx = int(geom.get_feature().get_name())
            lat, lon = geom.get_reconstructed_geometry().to_lat_lon()
            reconstructed_lats[idx], reconstructed_lons[idx] = round(lat, 3), round(lon, 3)

    df['paleolongitude'], df['paleolatitude'] = reconstructed_lons, reconstructed_lats
    output_path = input_path.parent / f"{input_path.stem}_with_paleo_coordinate.csv"
    df.to_csv(output_path, index=False)
    print(f"Successfully saved results to: {output_path}")

if __name__ == '__main__':
    BASE_DIR = Path(__file__).parent
    INPUT_CSV = BASE_DIR / 'test_dataset' / 'Template_Site982.csv'
    ROT_FILE = BASE_DIR / 'GPlates_features' / 'Rotations' / 'Zahirovic_etal_2022_CombinedRotations.rot'
    POLY_FILE = BASE_DIR / 'GPlates_features' / 'StaticPolygons' / 'Zahirovic_etal_2022_Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpmlz'
    reconstruct_paleo_coordinates(INPUT_CSV, ROT_FILE, POLY_FILE)
