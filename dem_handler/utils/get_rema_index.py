import os
import zipfile
from urllib.request import urlretrieve

REMA_INDEX_URL = "https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/REMA_Mosaic_Index_latest_gdb.zip"

def get_rema_index_file(save_folder: str) -> str:

    filename = "REMA_Mosaic_Index_latest_gdb.zip"
    # download and store locally

    os.makedirs(save_folder, exist_ok=True)
    zip_save_path = os.path.join(save_folder, filename)
    urlretrieve(REMA_INDEX_URL, zip_save_path)
    # unzip
    with zipfile.ZipFile(zip_save_path, "r") as zip_ref:
        zip_ref.extractall(save_folder)
        files = zip_ref.infolist()
        rema_index_file = os.path.dirname(files[0].filename)
    rema_index_path = os.path.join(save_folder, rema_index_file)
    os.remove(zip_save_path)
    return rema_index_path
