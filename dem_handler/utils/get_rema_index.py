import os
import zipfile
from urllib.request import urlretrieve

def get_REMA_index_file(save_folder):
    rema_index_url = 'https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/REMA_Mosaic_Index_latest_gdb.zip'
    filename = 'REMA_Mosaic_Index_latest_gdb.zip'
    # download and store locally
    zip_save_path = os.path.join(save_folder, filename)
    urlretrieve(rema_index_url, zip_save_path)
    #unzip 
    with zipfile.ZipFile(zip_save_path, 'r') as zip_ref:
        zip_ref.extractall(save_folder)
        files=zip_ref.infolist()
        rema_index_file = '/'.join(files[0].filename.split('/')[0:-1])
    rema_index_path = os.path.join(save_folder, rema_index_file)
    os.remove(zip_save_path)
    return rema_index_path
