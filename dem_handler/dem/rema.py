from pathlib import Path
from shapely import box
import geopandas as gpd
import rasterio as rio
from rasterio.merge import merge
from rasterio.profiles import Profile
import numpy as np
import shutil

from dem_handler.utils.spatial import BoundingBox, transform_polygon
from dem_handler.download.aws import download_rema_tiles

from dem_handler.dem.geoid import remove_geoid
from dem_handler.download.aws import download_egm_08_geoid
    

# Create a custom type that allows use of BoundingBox or tuple(xmin, ymin, xmax, ymax)
BBox = BoundingBox | tuple[float | int, float | int, float | int, float | int]

DATA_DIR = Path(__file__).parents[1] / Path('data')
REMA_GPKG_PATH = DATA_DIR / Path('REMA_Mosaic_Index_v2.gpkg')
REMA_VALID_RESOLUTIONS = [
    2,
    10,
    32,
]  # [2, 10, 32, 100, 500, 1000] It seems there are no higher resolutions in the new index

def get_rema_dem_for_bounds(
    bounds: BBox,
    save_path: str,
    rema_index_path: str = REMA_GPKG_PATH,
    resolution: int = 2,
    bbox_src_crs: int = 3031,
    bbox_dst_crs: int = 3031,
    ellipsoid_heights: bool = True,
    geoid_tif_path: Path = 'egm_08_geoid.tif',
    download_geoid: bool =  False,
) -> tuple[np.ndarray, Profile]:

    TEMP_SAVE_FOLDER = "rema_dems_temp_folder"
    GEOID_CRS = 4326

    assert (
        resolution in REMA_VALID_RESOLUTIONS
    ), f"resolution must be in {REMA_VALID_RESOLUTIONS}"

    if type(bounds) != BoundingBox:
        bounds = BoundingBox(*bounds)

    if bbox_src_crs != bbox_dst_crs:
        bounds_poly = transform_polygon(box(*bounds), bbox_src_crs, bbox_dst_crs)
    else:
        bounds_poly = box(*bounds)

    rema_layer = f"REMA_Mosaic_Index_v2_{resolution}m"
    rema_index_df = gpd.read_file(rema_index_path, layer=rema_layer)

    intersecting_rema_files = rema_index_df[
        rema_index_df.geometry.intersects(bounds_poly)
    ]
    s3_url_list = intersecting_rema_files["s3url"].to_list()
    print(f"{len(s3_url_list)} intersecting tiles found")
    dem_paths = download_rema_tiles(s3_url_list[0:], TEMP_SAVE_FOLDER)

    print("combining found DEMS")
    merge(dem_paths, dst_path=save_path)
    shutil.rmtree(TEMP_SAVE_FOLDER, ignore_errors=True)
    dem_raster = rio.open(save_path)
    dem_profile = dem_raster.profile

    if ellipsoid_heights:
        print(
            f"Subtracting the geoid from the DEM to return ellipsoid heights"
        )
        if not download_geoid and not Path(geoid_tif_path).exists():
            raise FileExistsError(f'Geoid file does not exist: {geoid_tif_path}. '\
                                    'correct path or set download_geoid = True'
                                    )
        elif download_geoid and not Path(geoid_tif_path).exists():
            print(f'Downloading the egm_08 geoid')
            geoid_bounds = bounds
            if bbox_src_crs != GEOID_CRS:
                geoid_bounds = transform_polygon(box(*bounds), bbox_src_crs, GEOID_CRS).bounds
            download_egm_08_geoid(geoid_tif_path, geoid_bounds)
        
        print(f"Using geoid file: {geoid_tif_path}")
        dem_array = remove_geoid(
            dem_array=dem_array,
            dem_profile=dem_profile,
            geoid_path=geoid_tif_path,
            buffer_pixels=2,
            save_path=save_path,
        )
    dem_raster = rio.open(save_path)
    dem_profile = dem_raster.profile

    return dem_raster.read(1), dem_profile
