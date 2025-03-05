# dem-handler
Utilility package for handling various Digital Elevation Models (DEMs). The package enables access to data stored in cloud as well as local copies of dem datasets. 

## Supported DEMS
- Copernicus Global 30m (cop_glo30)
- REMA (2m,10m, ...)

## Usage
### Download tiles

```python
import os
from dem_handler.dem.cop_glo30 import get_cop30_dem_for_bounds
from dem_handler.dem.rema import get_rema_dem_for_bounds

import logging
logging.basicConfig(level=logging.INFO)

#Set the bounds and make a directory for the files to download

bounds = (72,-70, 73, -69)
save_dir = 'TMP'
os.makedirs(save_dir, exist_ok=True)

#The copernicus Global 30m DEM 

get_cop30_dem_for_bounds(
    bounds = bounds,
    save_path = f'{save_dir}/cop_glo30.tif',
    adjust_at_high_lat=False,
    cop30_folder_path = save_dir,
    ellipsoid_heights = False,
    download_dem_tiles = True
)

# The REMA DEM (32m)

get_rema_dem_for_bounds(
    bounds=bounds,
    save_path=f'{save_dir}/rema.tif',
    resolution=32,
    bounds_crs=4326,
    rema_folder_path='TMP'
)
```

## Developer Setup

```bash
git clone ...
conda create --file environment.yaml
pip install -e .
```

Test install

```bash
pytest
```

## Contact
...