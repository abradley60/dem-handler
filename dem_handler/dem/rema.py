from pathlib import Path
from shapely import Polygon, box
import geopandas as gpd
from rasterio.merge import merge

from dem_handler.utils.spatial import BoundingBox, transform_polygon
from dem_handler.download.aws import download_REMA_tiles

# Create a custom type that allows use of BoundingBox or tuple(xmin, ymin, xmax, ymax)
BBox = BoundingBox | tuple[float | int, float | int, float | int, float | int]

DATA_DIR = Path(__file__).parents[1] / Path('data')
REMA_GPKG_PATH = DATA_DIR / Path('REMA_Mosaic_Index_v2.gpkg')

def get_rema_dem_for_bounds(
    bounds,
    save_path: Path,
    bounds_crs: int = 4326,
    resolution: int = 32,
    #ellipsoid_heights: bool = True,
    rema_index_path: Path = REMA_GPKG_PATH,
    rema_folder_path: Path = 'TMP',
    #geoid_tif_path: Path = 'egm_08_geoid.tif',
    #download_dem_tiles: bool = True,
    #download_geoid: bool =  True,
):
    if bounds_crs == 4326:
        # convert to 3031 for rema
        bounds_poly = (box(*bounds))
        bounds_poly_3031 = transform_polygon(bounds_poly, 4326, 3031)

    # res must be one of [2, 10, 32, 100, 500, 1000]
    
    # load into gpdf
    rema_index_df = gpd.read_file(rema_index_path)
    
    # find the intersecting tiles
    intersecting_rema_files = rema_index_df[
        rema_index_df.geometry.intersects(bounds_poly_3031)]
    s3_url_list = intersecting_rema_files['s3url'].to_list()
    print(f'{len(s3_url_list)} intersecting tiles found')
    dem_paths = download_REMA_tiles(s3_url_list[0:], resolution, save_folder=rema_folder_path)

    # note logic here should crop to bounds simillar to cop30 logic
    print('combining DEMS')
    merge(
            sources=dem_paths,
            dst_path=save_path,
            bounds=bounds_poly_3031.bounds
            )