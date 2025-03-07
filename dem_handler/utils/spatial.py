from dataclasses import dataclass
import pyproj
from shapely import segmentize
from shapely.geometry import Polygon, box
from pyproj.database import query_utm_crs_info
from pyproj.aoi import AreaOfInterest
from pyproj import CRS
from osgeo import gdal
import shapely
import rasterio
from rasterio.io import DatasetReader
from rasterio.profiles import Profile
import json
import logging
import numpy as np
import os
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)


# Construct a dataclass for bounding boxes
@dataclass
class BoundingBox:
    xmin: float | int
    ymin: float | int
    xmax: float | int
    ymax: float | int

    @property
    def bounds(self) -> tuple[float | int, float | int, float | int, float | int]:
        return (self.xmin, self.ymin, self.xmax, self.ymax)

    @property
    def top_left(self) -> tuple[float | int, float | int]:
        return (self.xmin, self.ymax)

    @property
    def bottom_right(self) -> tuple[float | int, float | int]:
        return (self.xmax, self.ymin)

    # Run checks on the bounding box values
    def __post_init__(self):
        if self.ymin >= self.ymax:
            raise ValueError(
                "The bounding box's ymin value is greater than or equal to ymax value. Check ordering"
            )
        if self.xmin >= self.xmax:
            raise ValueError(
                "The bounding box's x_min value is greater than or equal to x_max value. Check ordering"
            )


# Create a custom type that allows use of BoundingBox or tuple(xmin, ymin, xmax, ymax)
BBox = BoundingBox | tuple[float | int, float | int, float | int, float | int]


def transform_polygon(
    geometry: Polygon, src_crs: int, dst_crs: int, always_xy: bool = True
):
    src_crs = pyproj.CRS(f"EPSG:{src_crs}")
    dst_crs = pyproj.CRS(f"EPSG:{dst_crs}")
    transformer = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=always_xy)
    # Transform the polygon's coordinates
    if isinstance(geometry, Polygon):
        # Transform exterior
        exterior_coords = [
            transformer.transform(x, y) for x, y in geometry.exterior.coords
        ]
        # Transform interiors (holes)
        interiors_coords = [
            [transformer.transform(x, y) for x, y in interior.coords]
            for interior in geometry.interiors
        ]
        # Create the transformed polygon
        return Polygon(exterior_coords, interiors_coords)

    # Handle other geometry types as needed
    raise ValueError("Only Polygon geometries are supported for transformation.")


def adjust_bounds(
    bounds: BoundingBox | tuple[float | int, float | int, float | int, float | int],
    src_crs: int,
    ref_crs: int,
    segment_length: float = 0.1,
) -> tuple:
    """_summary_

    Parameters
    ----------
    bounds : BoundingBox | tuple[float | int, float | int, float | int, float | int],
        Bounds to adjust.
    src_crs : int
        Source EPSG. e.g. 4326
    ref_crs : int
        Reference crs to create the true bbox. i.e. 3031 in southern
        hemisphere and 3995 in northern (polar stereographic)
    segment_length : float, optional
        distance between generation points along the bounding box sides in
        src_crs. e.g. 0.1 degrees in lat/lon, by default 0.1

    Returns
    -------
    BoundingBox
        A polygon bounding box expanded to the true min max
    """
    if isinstance(bounds, tuple):
        bounds = BoundingBox(*bounds)

    geometry = box(*bounds.bounds)
    segmentized_geometry = segmentize(geometry, max_segment_length=segment_length)
    transformed_geometry = transform_polygon(segmentized_geometry, src_crs, ref_crs)
    transformed_box = box(*transformed_geometry.bounds)
    corrected_geometry = transform_polygon(transformed_box, ref_crs, src_crs)
    return BoundingBox(*corrected_geometry.bounds)


def get_local_utm(
    bounds: BoundingBox | tuple[float | int, float | int, float | int, float | int],
    antimeridian: bool = False,
) -> int:
    """_summary_

    Parameters
    ----------
    bounds : BoundingBox | tuple[float  |  int, float  |  int, float  |  int, float  |  int]
        The set of bounds (min_lon, min_lat, max_lon, max_lat)
    antimeridian : bool, optional
        Whether the bounds cross the antimerdian, by default False

    Returns
    -------
    int
        The CRS in integer form (e.g. 32749 for WGS 84 / UTM zone 49S)
    """
    if bounds.isinstance(tuple):
        bounds = BoundingBox(*bounds)

    logger.info("Finding best crs for area")
    centre_lat = (bounds.ymin + bounds.ymax) / 2
    centre_lon = (bounds.xmin + bounds.xmax) / 2
    if antimeridian:
        # force the lon to be next to antimeridian on the side with the scene centre.
        # e.g. (-177 + 178)/2 = 1, this is > 0 more data on -'ve side
        centre_lon = 179.9 if centre_lon < 0 else -179.9
    utm_crs_list = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=centre_lon - 0.01,
            south_lat_degree=centre_lat - 0.01,
            east_lon_degree=centre_lon + 0.01,
            north_lat_degree=centre_lat + 0.01,
        ),
    )
    crs = CRS.from_epsg(utm_crs_list[0].code)
    crs = str(crs).split(":")[-1]  # get the EPSG integer
    return int(crs)


def bounds_to_geojson(bounds, save_path=""):
    """
    Convert a bounding box to a GeoJSON Polygon.

    Parameters:
        bounds (tuple): A tuple of (min_lon, min_lat, max_lon, max_lat).
        save_path (str): path to save the geojson

    Returns:
        dict: A GeoJSON FeatureCollection with a single Polygon feature.
    """
    min_lon, min_lat, max_lon, max_lat = bounds

    # Define the polygon coordinates
    coordinates = [
        [
            [min_lon, min_lat],
            [min_lon, max_lat],
            [max_lon, max_lat],
            [max_lon, min_lat],
            [min_lon, min_lat],  # Closing the polygon
        ]
    ]

    # Create the GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": coordinates},
                "properties": {},
            }
        ],
    }

    if save_path:
        with open(save_path, "w") as f:
            f.write(json.dumps(geojson))

    return geojson


def check_s1_bounds_cross_antimeridian(bounds: BBox, max_scene_width: int = 20) -> bool:
    """Check if the s1 scene bounds cross the antimeridian. The bounds of a sentinel-1
    are valid at the antimeridian, just very large. By setting a max scene width, we
    can determine if the antimeridian is crossed. Alternate scenario is a bounds
    with a very large width (i.e. close to the width of the earth).

    Parameters
    ----------
    bounds : BoundingBox
        the set of bounds (xmin, ymin, xmax, ymax)
    max_scene_width : int, optional
        maximum allowable width of the scene bounds in degrees, by default 20

    Returns
    -------
    bool
        if the bounds cross the antimeridian
    """

    antimeridian_xmin = -180
    bounding_xmin = antimeridian_xmin + max_scene_width  # -160 by default

    antimeridian_xmax = 180
    bounding_xmax = antimeridian_xmax - max_scene_width  # 160 by default

    if (bounds.xmin < bounding_xmin) and (bounds.xmin > antimeridian_xmin):
        if bounds.xmax > bounding_xmax and bounds.xmax < antimeridian_xmax:
            return True
    return False


def get_target_antimeridian_projection(bounds: BoundingBox) -> int:
    """depending where were are on the earth, the desired
    crs at the antimeridian will change. e.g. polar stereo
    is desired at high and low lats, local utm zone elsewhere
    (e.g. at the equator).

    Parameters
    ----------
    bounds : BoundingBox
        The set of bounds (min_lon, min_lat, max_lon, max_lat)

    Returns
    -------
    int
        The CRS in integer form (e.g. 3031 for Polar Stereographic)
    """
    min_lat = min(bounds.ymin, bounds.ymax)
    target_crs = (
        3031
        if min_lat < -50
        else 3995 if min_lat > 50 else get_local_utm(bounds.bounds, antimeridian=True)
    )
    logger.warning(f"Data will be returned in EPSG:{target_crs} projection")
    return target_crs


def split_s1_bounds_at_am_crossing(
    bounds: BBox, lat_buff: float = 0
) -> tuple[BoundingBox]:
    """Split the s1 bounds at the antimeridian, producing one set of bounds for the
    Eastern Hemisphere (left of the antimeridian) and one set for the Western
    Hemisphere (right of the antimeridian)

    Parameters
    ----------
    bounds : BBox (BoundingBox | tuple[float | int, float | int, float | int, float | int])
        The set of bounds (xmin, ymin, xmax, ymax)
    lat_buff : float, optional
        An additional buffer to subract from lat, by default 0.

    Returns
    -------
    tuple[BoundingBox]
        A tuple containing two sets of bounds, one for the Eastern Hemisphere, one for
        the Western Hemisphere.
    """
    if isinstance(bounds, tuple):
        bounds = BoundingBox(*bounds)

    eastern_hemisphere_x = min([x for x in [bounds.xmin, bounds.xmax] if x > 0])
    if eastern_hemisphere_x > 180:
        raise ValueError(
            f"Eastern Hemisphere coordinate of {eastern_hemisphere_x} is more than 180 degrees, but should be less."
        )

    western_hemisphere_x = max([x for x in [bounds.xmin, bounds.xmax] if x < 0])
    if western_hemisphere_x < -180:
        raise ValueError(
            f"Western Hemisphere coordinate of {western_hemisphere_x} is less than -180 degrees, but should be greater."
        )

    min_y = max(-90, bounds.ymin - lat_buff)
    max_y = min(90, bounds.ymax + lat_buff)

    bounds_western_hemisphere = BoundingBox(-180, min_y, western_hemisphere_x, max_y)
    bounds_eastern_hemisphere = BoundingBox(eastern_hemisphere_x, min_y, 180, max_y)

    logger.info(f"Eastern Hemisphere bounds: {bounds_eastern_hemisphere.bounds}")
    logger.info(f"Western Hemisphere bounds: {bounds_western_hemisphere.bounds}")

    return (bounds_eastern_hemisphere, bounds_western_hemisphere)


def adjust_bounds_at_high_lat(bounds: BBox) -> tuple:
    """Expand the bounds for high lattitudes. The
    provided bounds sometimes do not contain the full scene due to
    warping at high latitudes. Solve this by converting bounds to polar
    steriographic, getting bounds, converting back to 4326. At high
    latitudes this will increase the longitude range.

    Parameters
    ----------
    bounds : BBox (BoundingBox | tuple[float | int, float | int, float | int, float | int])
        The set of bounds (min_lon, min_lat, max_lon, max_lat)

    Returns
    -------
    BoundingBox
        The expanded bounds (min_lon, min_lat, max_lon, max_lat)
    """
    if isinstance(bounds, tuple):
        bounds = BoundingBox(*bounds)

    if bounds.ymin < -50:
        logging.info(f"Adjusting bounds at high sourthern latitudes")
        bounds = adjust_bounds(bounds, src_crs=4326, ref_crs=3031)
    if bounds.ymin > 50:
        logging.info(f"Adjusting bounds at high northern latitudes")
        bounds = adjust_bounds(bounds, src_crs=4326, ref_crs=3995)

    return bounds


def resize_bbox(bbox: BoundingBox, scale_factor: float = 1.0) -> BoundingBox:
    """Resizes a bounding box

    Parameters
    ----------
    bbox : BoundingBox
        BoundingBox object.
    scale_factor : float, optional
        Factor to scale up or down the bounding box, by default 1.0

    Returns
    -------
    BoundingBox
        Resized bounding box.
    """
    x_dim = bbox.xmax - bbox.xmin
    y_dim = bbox.ymax - bbox.ymin

    dx = ((scale_factor - 1) * x_dim) / 2
    dy = ((scale_factor - 1) * y_dim) / 2

    return BoundingBox(bbox.xmin - dx, bbox.ymin - dy, bbox.xmax + dx, bbox.ymax + dy)


def crop_datasets_to_bounds(
    dem_rasters: list[str], bounds: BBox, save_path: str | Path | None = None
) -> tuple[np.ndarray, Profile]:
    """Merges a list of datasets and crops the merged tiles to a given bounding box.

    Parameters
    ----------
    dem_rasters : list[str]
        List of dataset paths or list of open datasets.
    bounds : BBox
        BoundingBox object or tuple of coordinates.
    save_path : str | Path | None, optional
        Local path to save the output merged tile, by default None

    Returns
    -------
    tuple[np.ndarray, Profile]
        tuple of array of merged tile and its profile.
    """

    if type(bounds) == BoundingBox:
        bounds = bounds.bounds

    logger.info(f"Creating VRT")
    vrt_path = (
        str(save_path).replace(".tif", ".vrt") if save_path else "temp_dem_file.vrt"
    )  # Temporary VRT file path
    logger.info(f"VRT path = {vrt_path}")
    VRT_options = gdal.BuildVRTOptions(
        resolution="highest",
        outputBounds=bounds,
        VRTNodata=0,
    )
    ds = gdal.BuildVRT(vrt_path, dem_rasters, options=VRT_options)
    ds.FlushCache()

    with rasterio.open(vrt_path, "r", count=1) as src:
        dem_array, dem_transform = rasterio.mask.mask(
            src,
            [shapely.geometry.box(*bounds)],
            all_touched=True,
            crop=True,
        )
        # Using the masking adds an extra dimension from the read
        # Remove this by squeezing before writing
        dem_array = dem_array.squeeze()
        logger.info(f"Dem array shape = {dem_array.shape}")

        dem_profile = src.profile
        dem_profile.update(
            {
                "driver": "GTiff",
                "height": dem_array.shape[0],
                "width": dem_array.shape[1],
                "transform": dem_transform,
                "count": 1,
                "nodata": np.nan,
            }
        )
        os.remove(vrt_path)

        if save_path:
            with rasterio.open(save_path, "w", **dem_profile) as dst:
                dst.write(dem_array, 1)
            # shutil.rmtree(
            #     dem_rasters[0].split("/")[0].split("\\")[0], ignore_errors=True
            # )

    return dem_array, dem_profile
