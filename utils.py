import json
import os
import math
import time
import shutil
import geojson
import numpy
import zipfile
import platform
import ctypes
from typing import List
from collections import namedtuple

import shapely.ops
from shapely import wkt
from pyproj import Transformer
from shapely.geometry import Point
import data_convert as data_convert
# from helper.transform import wgs84togcj02

# convert nds coordnate to wgs84 coordinate
nds_to_wgs84_coordinate = lambda x: x * 360.0 / 4294967296
# convert wgs84 coordinate to nds coordinate
wgs84_to_nds_coordinate = lambda x: int((x / 360.0) * 4294967296)

MEASURED_LAYERS = ['Boundary', 'Line', 'Mark', 'Sign', 'TL', 'Overhead', 'OverHead', 'Arrow',
                   'HADRoadDivider', 'HADLaneDivider', 'LandMark', 'LocSign', 'LocTrafficLight', 'HADLaneNode']
TRAJ_PROTO_DIR = 'trajectory'
TRAJ_GEOJSON_DIR = 'TrajJsonData'


def wgs84_get_utm_epsg(lon, lat):
    epsg = 32700 - numpy.round((45 + lat) / 90, 0) * 100 + numpy.round((183 + lon) / 6, 0)
    return epsg


epsg_cache = {}


def wgs84_2_utm(longitude, latitude):
    epsg = wgs84_get_utm_epsg(longitude, latitude)
    if epsg in epsg_cache:
        wgs84_to_utm, utm_to_wgs84 = epsg_cache[epsg]
    else:
        wgs84_to_utm = Transformer.from_crs('epsg:4326', 'epsg:' + str(epsg))
        utm_to_wgs84 = Transformer.from_crs('epsg:' + str(epsg), 'epsg:4326')
        epsg_cache[epsg] = (wgs84_to_utm, utm_to_wgs84)

    x, y = wgs84_to_utm.transform(latitude, longitude)
    return x, y


def calc_distance(lon1, lat1, lon2, lat2):
    x1, y1 = wgs84_2_utm(lon1, lat1)
    x2, y2 = wgs84_2_utm(lon2, lat2)
    distance = math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
    return distance


def ogr_to_shapely(feature):
    geometry = feature.geometry()
    if not geometry:
        return None
    geometry.CloseRings()
    geometry = geometry.ExportToWkt()
    geometry = wkt.loads(geometry)
    # geometry.rowid = id
    return geometry


def extract_zip_file(file_path):
    # minio下载文件是zip格式,解压到当前目录
    zipfile_path = file_path + '.zip'
    f = zipfile.ZipFile(zipfile_path, 'r')
    for file in f.namelist():
        f.extract(file, file_path)
    f.close()
    os.remove(zipfile_path)


def split_bucket_path(path):
    bucket = path.split('/')[0]
    obj_list = path.split('/')[1:]
    key_path = '/'.join(obj_list)
    return bucket, key_path


def get_tile_ids(client, bucket: str, key_path: str):
    # 获取到当前bucket指定路径下的tildes
    objects = client.list_objects(bucket, key_path, recursive=True)
    file_detail_list: List[namedtuple] = objects
    tile_ids = []
    for detail in file_detail_list:
        tile_id = detail.object_name.split('/')[:-1]
        path = '/'.join(tile_id)
        if tile_id not in tile_ids:
            tile_ids.append(path)
    return tile_ids


def get_base_path(strtime):
    if platform.system() == 'Windows':
        groundtruth_path = 'pretreatment' + '\\' + 'data' + '\\' + 'data' + strtime
    else:
        groundtruth_path = 'pretreatment' + '/' + 'data' + '/' + 'data' + strtime
    return groundtruth_path


def get_truth_path(base_path):
    if platform.system() == 'Windows':
        groundtruth_path = base_path + '\\' + 'groundtruth'
    else:
        groundtruth_path = base_path + '/' + 'groundtruth'
    return groundtruth_path


def get_target_path(base_path, product_type):
    if platform.system() == 'Windows':
        target_path = base_path + '\\' + product_type
    else:
        target_path = base_path + '/' + product_type
    return target_path


def get_truth_area_path(base_path):
    if platform.system() == 'Windows':
        target_path = base_path + '\\' + 'truth_area.geojson'
    else:
        target_path = base_path + '/' + 'truth_area.geojson'
    return target_path


def get_target_area_path(base_path):
    if platform.system() == 'Windows':
        target_path = base_path + '\\' + 'target_area.geojson'
    else:
        target_path = base_path + '/' + 'target_area.geojson'
    return target_path


def get_import_result_path(base_path):
    if platform.system() == 'Windows':
        import_result_path = base_path + '\\' + 'import_result'
    else:
        import_result_path = base_path + '/' + 'import_result'
    return import_result_path


def get_filter_path(base_path):
    if platform.system() == 'Windows':
        target_path = base_path + '\\' + 'filter'
    else:
        target_path = base_path + '/' + 'filter'
    return target_path


def get_car_result_path(base_path):
    if platform.system() == 'Windows':
        car_path = base_path + '\\' + 'car_result' + '\\' + 'car.geojson'
    else:
        car_path = base_path + '/' + 'car_result' + '/' + 'car.geojson'
    return car_path


def get_car_trajectory_file_path(base_path):
    if platform.system() == 'Windows':
        car_path = base_path + '\\' + 'car_result' + '\\' + 'trajectory-all.json'
    else:
        car_path = base_path + '/' + 'car_result' + '/' + 'trajectory-all.json'
    return car_path


def get_car_trajectory_area_path(base_path):
    if platform.system() == 'Windows':
        car_path = base_path + '\\' + 'car_result' + '\\' + 'trajectory-area.geojson'
    else:
        car_path = base_path + '/' + 'car_result' + '/' + 'trajectory-area.geojson'
    return car_path


def get_car_truth_area_path(base_path):
    if platform.system() == 'Windows':
        car_path = base_path + '\\' + 'car_result' + '\\' + 'truth-area.geojson'
    else:
        car_path = base_path + '/' + 'car_result' + '/' + 'truth-area.geojson'
    return car_path


def get_area_path(root, index):
    if platform.system() == 'Windows':
        area_path = root + '\\' + 'area' + str(index) + ".geojson"
    else:
        area_path = root + '/' + 'area' + str(index) + ".geojson"
    return area_path


def is_truth_dat(truth_path):
    """
    judge is truth file .dat ?
    """
    if os.path.isfile(truth_path):
        file_type = os.path.basename(truth_path).split('.')[-1]
        if file_type == 'dat':
            return True
    else:
        for root, _, files in os.walk(truth_path):
            for file in files:
                file_type = file.split('.')[-1]
                if file_type == 'dat':
                    return True


def calc_runtime(desc):
    def wrapper(func):
        def inner(*args, **kwargs):
            start = time.time()
            res = func(*args, **kwargs)
            end = time.time()
            print(f'{desc}运行耗时: {end - start} 秒')
            return res

        return inner

    return wrapper


@calc_runtime(desc='获取轨迹范围面耗时')
def get_cardata_traj_area(cardata_path, traj_file_name):
    """检查待测数据下是否包含车端轨迹数据，若包含，则进行解析数据并计算轨迹面后清空目录，将polygon存成geojson文件"""
    # traj_x = []
    # traj_y = []
    traj_areas = []
    traj_feature_lst = []
    for root, _, files in os.walk(cardata_path):
        if files and root.split(os.sep)[-1] == TRAJ_PROTO_DIR:
            converter = data_convert.CarProtoConvert()
            for file in files:
                traj_file_path = os.path.join(root, file)
                traj_info = converter.read_location(traj_file_path)
                traj_points = []
                for traj in traj_info:
                    traj_features = traj[1]
                    geom = Point(traj_features.get('geometry')['coordinates'])
                    # traj_points.append(Point(wgs84togcj02(geom.x, geom.y)))
                    traj_points.append(Point(geom.x, geom.y))
                traj_area = shapely.geometry.LineString(traj_points).buffer(8.983152841195234e-06 * 6)
                traj_points.clear()
                os.remove(traj_file_path)  # 删除record
                feature = geojson.Feature(geometry=None, properties={})
                feature['geometry'] = json.loads(geojson.dumps(traj_area))
                traj_feature_lst.append(feature)
        elif root != cardata_path:
            shutil.rmtree(root)
    traj_area_file = os.path.join(cardata_path, traj_file_name)
    with open(traj_area_file, "w", encoding='utf-8') as fp:
        fp.write(json.dumps(geojson.FeatureCollection(traj_feature_lst)))
    # # 根据轨迹点拟合线
    # line_equation = numpy.polyfit(traj_x, traj_y, 3)
    # line1 = np.poly1d(line_equation)
    # y_list = line1(traj_x)
    # traj_line = list(zip(traj_x, y_list))


def transform_wgs84_to_gcj02(lon, lat):
    if platform.system() == 'Windows':
        shared = ctypes.cdll.LoadLibrary(r'pretreatment\\lib\\transform.dll')
    else:
        shared = ctypes.cdll.LoadLibrary('pretreatment/lib/TransformCoordinate.so')

    shared.transform_84_to_02.argtypes = [ctypes.c_double, ctypes.c_double]
    shared.transform_84_to_02.restype = ctypes.POINTER(ctypes.c_double * 2)

    r = shared.transform_84_to_02(lon, lat)
    new_lon, new_lat = r.contents[0], r.contents[1]

    return new_lon, new_lat


if __name__ == '__main__':
    get_cardata_traj_area(r'I:\ca data\02_车端矢量\2023-06-19\LS6A2E161NA505442_L17_2s',\
                           r'E:\dataca\temp')
    # get_cardata_traj_area(r'C:\Users\202207817\Downloads\车端路网')
    # get_cardata_traj_area(r'F:\BaiduSyncdisk\博士在读期间\02_质量子任务\10_07.11张志军提供数据\车端成图_L17_5s\xxxx_L17_5s',\
    #                       r'F:\BaiduSyncdisk\博士在读期间\02_质量子任务\10_07.11张志军提供数据\车端成图_L17_5s\xxxx_L17_5s\TrajJsonData2')
