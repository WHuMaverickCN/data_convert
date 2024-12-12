import os
import zlib
import json
import decimal
import geojson
import numpy as np
# from protobuf_to_dict import protobuf_to_dict
from protobuf_to_dict import protobuf_to_dict
from shapely import geometry as shapelyGeo
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import LineString
# import pretreatment.proto.DataCollection_pb2 as DataCollection_pb2
# import pretreatment.proto.TileData_pb2 as tile_data_pb
# import pretreatment.proto.Tencent_TileData_pb2 as tencent_tile_data_pb
# import pretreatment.utils as utils
import proto.DataCollection_pb2 as DataCollection_pb2
import proto.TileData_pb2 as tile_data_pb
import proto.Tencent_TileData_pb2 as tencent_tile_data_pb
import utils as utils
import time

cur_tmp = int(time.time() * 1000000)

_LONG_INT_NAX = 2 ** 63 - 1

class CarProtoConvert(object):
    """
        车端数据proto转geojson类
    """
    def __init__(self):
        self.precision1 = 100000000
        self.precision2 = 1
        self.max_x = 135.083333
        self.min_x = 73.55
        self.max_y = 53.55
        self.min_y = 3.85
        self.max_z = 8844450
        self.min_z = -154320

        self.detection_temp = []
        self.senmatic_list = []

    def get_rect_polygon(self, bounds):
        """
        根据元组(最小x，最小y，最大x，最大y)获取矩形框
        """
        min_x, min_y, min_z, max_x, max_y, max_z = bounds
        min_z = min_z / 1000
        max_z = max_z / 1000
        polygon_list = [[min_x, min_y, min_z], [max_x, min_y, max_z],
                        [max_x, max_y, max_z], [min_x, max_y, min_z],
                        [min_x, min_y, min_z]]

        return polygon_list

    def xcode_to_bit32(self, t):
        if '2' <= t <= '9':
            return ord(t) - ord('2')
        if 'A' <= t <= 'H':
            return ord(t) - ord('A') + 8
        if 'J' <= t <= 'N':
            return ord(t) - ord('J') + 16
        if 'P' <= t <= 'Z':
            return ord(t) - ord('P') + 21

        return 255

    def geo_range_check(self, x_point):
        if self.max_x < x_point.x or x_point.x < self.min_x:
            return False
        if self.max_y < x_point.y or x_point.y < self.min_y:
            return False
        if self.max_z < x_point.z or x_point.z < self.min_z:
            return False
        return True

    def convert_code_to_geo(self, m_code):
        if not m_code or 18 != len(m_code):
            return False
        x = ((((self.xcode_to_bit32(m_code[4]) & 0x01) << 2) | (
                self.xcode_to_bit32(m_code[5]) >> 3)) * 32 * 32 * 32 * 32 * 32 * 32) + \
            ((self.xcode_to_bit32(m_code[12])) * 32 * 32 * 32 * 32 * 32) + \
            ((self.xcode_to_bit32(m_code[13])) * 32 * 32 * 32 * 32) + \
            ((self.xcode_to_bit32(m_code[14])) * 32 * 32 * 32) + \
            ((self.xcode_to_bit32(m_code[15])) * 32 * 32) + \
            ((self.xcode_to_bit32(m_code[16])) * 32) + \
            (self.xcode_to_bit32(m_code[17]))
        y = (self.xcode_to_bit32(m_code[5]) & 0x07) * 32 * 32 * 32 * 32 * 32 * 32 + \
            (self.xcode_to_bit32(m_code[6])) * 32 * 32 * 32 * 32 * 32 + \
            (self.xcode_to_bit32(m_code[7])) * 32 * 32 * 32 * 32 + \
            (self.xcode_to_bit32(m_code[8])) * 32 * 32 * 32 + \
            (self.xcode_to_bit32(m_code[9])) * 32 * 32 + \
            (self.xcode_to_bit32(m_code[10])) * 32 + \
            (self.xcode_to_bit32(m_code[11]))
        z = (self.xcode_to_bit32(m_code[0])) * 32 * 32 * 32 * 32 * 32 + (
            self.xcode_to_bit32(m_code[1])) * 32 * 32 * 32 * 32 + (
                self.xcode_to_bit32(m_code[2])) * 32 * 32 * 32 + (
                self.xcode_to_bit32(m_code[3])) * 32 * 32 + (self.xcode_to_bit32(m_code[4])) * 32 + (
                self.xcode_to_bit32(m_code[5]))

        z >>= 6
        x = x / self.precision1 + self.min_x
        y = y / self.precision1 + self.min_y
        z = z / self.precision2 + self.min_z
        x_point = Point((x, y, z))
        if self.geo_range_check(x_point):
            return x_point
        else:
            return None

    def convert_feature(self, geometry, properties):
        feature = geojson.Feature(geometry=None, properties=properties)
        feature['geometry'] = json.loads(geojson.dumps(geometry))
        return feature

    def get_semantic_property(self, obj_data):

        property_list = {}

        # match obj_data.type:
        #     case DataCollection_pb2.ObjectType.ROAD_SURFACE_LINE:
        #         subtype = obj_data.lineType
        #     case DataCollection_pb2.ObjectType.ROAD_SURFACE_ARROW:
        #         subtype = obj_data.arrowType
        #     case DataCollection_pb2.ObjectType.ROAD_SURFACE_MARK:
        #         subtype = obj_data.markType
        #     case DataCollection_pb2.ObjectType.ROAD_SIGN:
        #         subtype = obj_data.signType
        #     case DataCollection_pb2.ObjectType.ROAD_TRAFFIC_LIGHT:
        #         subtype = obj_data.trafficLightType
        #     case DataCollection_pb2.ObjectType.ROAD_BOUNDARY:
        #         subtype = obj_data.boundaryType
        #     case DataCollection_pb2.ObjectType.ROAD_OVERHEAD:
        #         subtype = obj_data.overheadType
        #     case DataCollection_pb2.ObjectType.ROAD_POLE:
        #         subtype = obj_data.poleType
        #     case _:
        #         subtype = 1000
        if obj_data.type == DataCollection_pb2.ObjectType.ROAD_SURFACE_LINE:
            subtype = obj_data.lineType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_SURFACE_ARROW:
            subtype = obj_data.arrowType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_SURFACE_MARK:
            subtype = obj_data.markType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_SIGN:
            subtype = obj_data.signType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_TRAFFIC_LIGHT:
            subtype = obj_data.trafficLightType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_BOUNDARY:
            subtype = obj_data.boundaryType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_OVERHEAD:
            subtype = obj_data.overheadType
        elif obj_data.type == DataCollection_pb2.ObjectType.ROAD_POLE:
            subtype = obj_data.poleType
        else:
            subtype = 1000
        property_list['oid'] = str(obj_data.id)
        property_list['color'] = obj_data.color
        property_list['type'] = obj_data.type
        property_list['type_confidence'] = obj_data.type_confidence
        property_list['subtype_confidence'] = obj_data.subtype_confidence
        property_list['value'] = obj_data.value
        property_list['value_confidence'] = obj_data.value_confidence
        property_list['subtype'] = subtype
        property_list['longitudinal_typ'] = obj_data.lineLongitudinalType

        return property_list

    def convert_semantic(self, semantic_map_list, dir_name):
        for semantic_map in semantic_map_list.semantic_map_info_list:
            vec_fs = []
            obj_fs = []
            base_info = semantic_map.base_info
            start_time = base_info.start_time
            # end_time = start_time + base_info.work_time
            end_time_list = [obj_data.time_stamp.start_time + obj_data.time_stamp.work_time
                             for obj_data in semantic_map.objdata_list]
            if not end_time_list:
                continue
            end_time = max(end_time_list)
            if start_time < 0 or end_time < 0 or start_time > cur_tmp or end_time > cur_tmp:
                continue
            # print(semantic_map.base_info)

            for obj_data in semantic_map.objdata_list:
                point_list = []
                geometry = None

                if obj_data.HasField('outline'):
                    for m_point in obj_data.outline.point:
                        point = self.convert_code_to_geo(m_point.mcode)
                        if not point:
                            continue
                        point_list.append([point.x, point.y, point.z/1000])
                    geometry = Polygon(point_list)
                elif obj_data.HasField('line'):
                    for m_point in obj_data.line.point:
                        point = self.convert_code_to_geo(m_point.mcode)
                        point_list.append([point.x, point.y, point.z/1000])
                    geometry = LineString(point_list)
                elif obj_data.HasField('rect'):
                    lb = self.convert_code_to_geo(obj_data.rect.lb.mcode)
                    rt = self.convert_code_to_geo(obj_data.rect.rt.mcode)
                    if not lb and not rt:
                        continue
                    bounds = (lb.x, lb.y, lb.z,
                                rt.x, rt.y, rt.z)
                    point_list = self.get_rect_polygon(bounds)
                    geometry = Polygon(point_list)
                elif obj_data.HasField('center_point'):
                    mm_geometry = self.convert_code_to_geo(obj_data.center_point.mcode)
                    geometry = Point(mm_geometry.x, mm_geometry.y, mm_geometry.z/1000)
                # elif obj_data.HasField('solid'):
                #     top_point = self.convert_code_to_geo(obj_data.solid.point_top.mcode)
                #     bottom_point = self.convert_code_to_geo(obj_data.solid.point_bottom.mcode)
                #     point_list.append([top_point.x, top_point.y, top_point.z])
                #     point_list.append([bottom_point.x, bottom_point.y, bottom_point.z])
                #     geometry = MultiPoint(point_list)

                if not geometry:
                    continue

                property_list = self.get_semantic_property(obj_data)
                # if obj_data.source == DataCollection_pb2.Source.DETECTION:
                property_list['start_time'] = start_time
                property_list['end_time'] = end_time
                feature = self.convert_feature(geometry=geometry, properties=property_list)

                if obj_data.source == DataCollection_pb2.Source.DETECTION:
                    obj_fs.append(feature)
                elif obj_data.source == DataCollection_pb2.Source.SEMANTIC_RECOGNITION:
                    vec_fs.append(feature)

            # file_name = str(start_time) + '_' + str(end_time) + '.geojson'
            # obj_path = os.path.join(dir_name, 'ObjJsonData')
            # vec_path = os.path.join(dir_name, 'VecJsonData')

            # if len(obj_fs) > 0:
            #     print('obj: ' + str(start_time)  + ' ' + str(end_time))
            #     fc = geojson.FeatureCollection(obj_fs)
            #     file_path = os.path.join(obj_path, file_name)
            #     # print(file_path)
            #     with open(file_path, "w", encoding='utf-8') as file:
            #         file.write(json.dumps(fc))

            if len(obj_fs) > 0:
                self.detection_temp.append((start_time, end_time, obj_fs))

            if len(vec_fs) > 0:
                # print('vec: ' + str(start_time)  + ' ' + str(end_time))
                # fc = geojson.FeatureCollection(obj_fs + vec_fs)
                # file_path = os.path.join(vec_path, file_name)

                # obj_fs.clear()
                # with open(file_path, "w", encoding='utf-8') as file:
                #     file.write(json.dumps(fc))

                self.senmatic_list.append((start_time, end_time, vec_fs))

    def convert_semantic_data(self, path, dir_path):
        with open(path, "rb") as file:
            buf = file.read()
        collection_info = DataCollection_pb2.DataCollectInfo()
        collection_info.ParseFromString(buf)

        if collection_info.data_type == DataCollection_pb2.SEMANTIC_MAP:
            semantic_map_list = DataCollection_pb2.SemanticMapList()
            semantic_map_list.ParseFromString(collection_info.collect_data)
            self.convert_semantic(semantic_map_list, dir_path)

    def get_trajectory_property(self, location_info):
        property_list = {}

        property_list['pitch'] = location_info.pitch
        property_list['roll'] = location_info.roll
        property_list['timestamp'] = location_info.timestamp
        property_list['ve'] = location_info.neg_speed.east_velocity
        property_list['vn'] = location_info.neg_speed.north_velocity
        property_list['vg'] = location_info.neg_speed.ground_velocity
        property_list['azimuth'] = location_info.yaw

        return property_list

    def convert_trajectory(self, location_data, common_data, dir_name):
        fs = []
        for location_info in location_data.location_info_list:
            property_list = self.get_trajectory_property(location_info)

            mm_geometry = self.convert_code_to_geo(location_info.geo_point.mcode)
            geometry = Point(mm_geometry.x, mm_geometry.y, mm_geometry.z/1000)
            property_list['altitude'] = mm_geometry.z/1000

            feature = self.convert_feature(geometry=geometry, properties=property_list)
            fs.append(feature)

        base_info = location_data.base_info

        fc = geojson.FeatureCollection(fs)
        fc['car_id'] = common_data.id
        fc['model_id'] = common_data.model_id
        fc['start_time'] = base_info.start_time
        fc['duration_time'] = base_info.start_time + base_info.work_time
        traj = {}
        # traj['format_version'] = common_data.protocol_version
        traj[common_data.id] = fc

        file_name = 'trajectory_' + str(base_info.start_time) + '.geojson'
        file_path = os.path.join(dir_name, file_name)
        # print(file_path)

        with open(file_path, "w", encoding='utf-8') as file:
            file.write(json.dumps(traj))

    def convert_location_data(self, location_data, common_data, dir_name):

        location_list = []

        # base_info = location_data.base_info
        # print(base_info)

        for location_info in location_data.location_info_list:

            property_list = self.get_trajectory_property(location_info)

            mm_geometry = self.convert_code_to_geo(location_info.geo_point.mcode)
            geometry = Point(mm_geometry.x, mm_geometry.y, mm_geometry.z/1000)
            property_list['altitude'] = mm_geometry.z/1000

            # print(mm_geometry.x, mm_geometry.y)

            feature = self.convert_feature(geometry=geometry, properties=property_list)
            location_list.append((location_info.timestamp, feature, common_data))

        return location_list

    def read_all_location(self, trajectory_path):
        all_traj_list = []
        for top_path, _, files in os.walk(trajectory_path):
            file_num = len(files)
            sep_num = file_num // 100
            sep_list = [100 * sep for sep in range(0, sep_num + 1)] + [file_num]
            for n in range(len(sep_list) - 1):
                sub_traj_list = []
                sub_files = files[sep_list[n]:sep_list[n + 1]]
                for file in sub_files:
                    file_path = os.path.join(top_path, file)
                    with open(file_path, "rb") as file:
                        buf = file.read()
                    collection_info = DataCollection_pb2.DataCollectInfo()
                    collection_info.ParseFromString(buf)
                    if collection_info.data_type == DataCollection_pb2.LOCATION_INFO:
                        location_data = DataCollection_pb2.LocationData()
                        location_data.ParseFromString(collection_info.collect_data)

                        traj_list = self.convert_location_data(location_data, collection_info.common_data, file_path)

                        sub_traj_list = sub_traj_list + traj_list
                        traj_list.clear()
                all_traj_list += sub_traj_list
                sub_traj_list.clear()
        return all_traj_list

    def read_location(self, trajectory_file):
        with open(trajectory_file, "rb") as file:
            buf = file.read()
        collection_info = DataCollection_pb2.DataCollectInfo()
        collection_info.ParseFromString(buf)

        if collection_info.data_type == DataCollection_pb2.LOCATION_INFO:
            location_data = DataCollection_pb2.LocationData()
            location_data.ParseFromString(collection_info.collect_data)
            return self.convert_location_data(location_data, collection_info.common_data, trajectory_file)

    def export_traj_file(self, common_data, traj_list, traj_path, start, end):
        fc = geojson.FeatureCollection(traj_list)
        fc['car_id'] = common_data.id
        fc['model_id'] = common_data.model_id
        fc['start_time'] = start
        fc['duration_time'] = end - start

        traj = {}
        # traj['format_version'] = common_data.protocol_version
        traj[common_data.id] = fc

        file_name = 'trajectory_' + str(start) + '.geojson'
        file_path = os.path.join(traj_path, file_name)

        with open(file_path, "w", encoding='utf-8') as file:
            file.write(json.dumps(traj))

    def rename_semantic_file(self, start, end, new_start, root):
        file_name = str(start) + '_' + str(end) + '.geojson'
        file_path = os.path.join(root, file_name)
        new_name = str(new_start) + '_' + str(end) + '.geojson'
        new_file_path = os.path.join(root, new_name)

        print(file_path + ' -> ' + new_file_path)

        os.rename(file_path, new_file_path)

    def convert_data(self, path, dir_path):
        with open(path, "rb") as file:
            buf = file.read()
        collection_info = DataCollection_pb2.DataCollectInfo()
        collection_info.ParseFromString(buf)
        if collection_info.data_type == DataCollection_pb2.SEMANTIC_MAP:
            semantic_map_list = DataCollection_pb2.SemanticMapList()
            semantic_map_list.ParseFromString(collection_info.collect_data)

            self.convert_semantic(semantic_map_list, dir_path)
        elif collection_info.data_type == DataCollection_pb2.LOCATION_INFO:
            location_data = DataCollection_pb2.LocationData()
            location_data.ParseFromString(collection_info.collect_data)

            traj_path = os.path.join(dir_path, 'TrajJsonData')
            self.convert_trajectory(location_data, collection_info.common_data, traj_path)
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


    @calc_runtime(desc='原始数据转换至geojson耗时')
    def convert_all_data(self, dir_path, download_path):
        # root_path = os.path.dirname(dir_path)
        root_path = dir_path

        if not os.path.exists(root_path):
            os.mkdir(root_path)

        # obj_path = os.path.join(root_path, 'ObjJsonData')
        vec_path = os.path.join(root_path, 'VecJsonData')
        traj_path = os.path.join(root_path, 'TrajJsonData')

        # if not os.path.exists(obj_path):
        #     os.mkdir(obj_path)

        if not os.path.exists(vec_path):
            os.mkdir(vec_path)

        if not os.path.exists(traj_path):
            os.mkdir(traj_path)

        # for root, _, files in os.walk(download_path):
        #     for file in files:
        #         file_path = os.path.join(root, file)
        #         convert_data(file_path, root_path)

        semantic_path = os.path.join(download_path, 'semantic')
        if not os.path.exists(semantic_path):
            print("error: no semantic data" + semantic_path)
            return

        trajectory_path = os.path.join(download_path, 'trajectory')
        if not os.path.exists(trajectory_path):
            print("error: no trajectory data")
            return

        for top_path, _, files in os.walk(semantic_path):
            for name in files:
                file_path = os.path.join(top_path, name)
                # print("序列化数据解析中..",name,end='\r')
                print(name,end='\r')
                self.convert_semantic_data(file_path, root_path)

        for _, end_time, obj_fs in self.detection_temp:
            for senmatic_start, senmatic_end, feature_list in self.senmatic_list:
                if end_time >= senmatic_start and end_time <= senmatic_end:
                    # print(str(start_time) + " in " + str(senmatic_start) + " - " + str(senmatic_end))
                    feature_list += obj_fs
                    break

        for senmatic_start, senmatic_end, feature_list in self.senmatic_list:
            print("写入车端语义矢量："+senmatic_start.__str__(),end="\r")
            fc = geojson.FeatureCollection(feature_list)
            file_name = str(senmatic_start) + '_' + str(senmatic_end) + '.geojson'
            # vec_path = os.path.join(dir_path, 'VecJsonData')
            file_path = os.path.join(vec_path, file_name)
            # print(file_path)
            with open(file_path, "w", encoding='utf-8') as handle:
                handle.write(json.dumps(fc))

        vec_timestamp_list = []
        for top_path, _, files in os.walk(vec_path):
            for file_name in files:
                pos =  file_name.rfind('.')
                name = file_name[:pos]
                items = name.split('_')
                vec_timestamp_list.append((int(items[0]), int(items[1])))

        vec_timestamp_list.sort(key=lambda x:x[0])

        all_locations = self.read_all_location(trajectory_path)
        all_locations.sort(key=lambda x:x[0])

        index = 0
        common_data = None
        fs = []
        traj_start = 0
        for key in vec_timestamp_list:
            start, end = key
            if start < 0 or end < 0 or start > cur_tmp or end > cur_tmp:
                continue
            for item in all_locations:
                timestamp = item[0]
                traj = item[1]
                if common_data is None:
                    common_data = item[2]

                if timestamp >= start and timestamp <= end:  # 表示轨迹点时间戳在语义时间戳范围内
                    if len(fs) == 0:
                        if timestamp > start:
                            # 表示第一个轨迹点时间戳大于第一个semantic start timestamp，修改semantic文件名
                            self.rename_semantic_file(start, end, timestamp, vec_path)
                            print(str(start) + ' -> ' + str(timestamp))
                            traj_start = timestamp
                        else:
                            traj_start = start
                    fs.append(traj)

            if len(fs) > 0:  # 语义时间戳范围内轨迹点数量大于0时输出对应轨迹文件
                self.export_traj_file(common_data, fs, traj_path, traj_start, end)
                print("写入车端轨迹矢量：",traj_start.__str__(),end = "\r")
            fs.clear()
            common_data = None
            index = index + 1
        ## 合并全量要素
        '''
        all_semantic_json_path = os.path.join(root_path, 'all_semantic_feature.geojson')
        all_traj_json_path = os.path.join(root_path, 'all_traj_feature.geojson')
        # 创建一个空的semantic GeoJSON 特征集
        all_semantic_features = {
            "type": "FeatureCollection",
            "features": []
        }

        # 遍历文件夹中的每个 GeoJSON 文件
        for file in os.listdir(vec_path):
            # 构建完整的文件路径
            file_path = os.path.join(vec_path, file)

             # 打开并读取当前 GeoJSON 文件
            semantic_info = json.loads(open(file_path).read())
            semantic_features = semantic_info['features']

            # 将当前 GeoJSON 文件中的每个特征添加到合并后的特征集中
            all_semantic_features['features'].extend(semantic_features)
        with open(all_semantic_json_path, 'w') as f1:
            json.dump(all_semantic_features, f1, indent=2)

        
        # 创建一个空的trajectory GeoJSON 特征集
        all_traj_features = {
            "type": "FeatureCollection",
            "features": []
        }

        # 遍历文件夹中的每个 GeoJSON 文件
        for file in os.listdir(traj_path):
            # 构建完整的文件路径
            file_path = os.path.join(traj_path, file)

            # 打开并读取当前 GeoJSON 文件
            traj_info = json.loads(open(file_path).read())
            traj_features = list(traj_info.values())[0]['features']

            # 将当前 GeoJSON 文件中的每个特征添加到合并后的特征集中
            all_traj_features['features'].extend(traj_features)
        with open(all_traj_json_path, 'w') as f2:
            json.dump(all_traj_features, f2, indent=2)

        '''

        '''
        # collapse
        # 下述部分采用将输出的geojson高程置为0
        # 读取 GeoJSON 文件
        with open(all_traj_json_path, 'r') as geojson_file:
            geojson_data = json.load(geojson_file)

        # 遍历每个要素，并将高程设置为 0
        for feature in geojson_data['features']:
            # 获取要素的几何信息
            geometry = feature['geometry']

            # 检查是否包含 'coordinates' 字段
            if 'coordinates' in geometry:
                # 获取坐标列表
                coordinates = geometry['coordinates']
                # 遍历坐标列表并将高程设置为 0
                coordinates[2] = 0
        # 将修改后的数据保存回 GeoJSON 文件
        with open(all_traj_json_path.replace('.geojson','_noH.geojson'), 'w') as output_file:
            json.dump(geojson_data, output_file, indent=2)
        print("已将语义高程设置为 0，并保存到 _noH.geojson 文件中")

        with open(all_semantic_json_path, 'r') as geojson_file:
            geojson_data = json.load(geojson_file)

        # 遍历每个要素，并将高程设置为 0
        for feature in geojson_data['features']:
            # 获取要素的几何信息
            geometry = feature['geometry']

            # 检查是否包含 'coordinates' 字段
            if 'coordinates' in geometry:
                # 获取坐标列表
                coordinates = geometry['coordinates']
                if geometry['type']=="Polygon":
                # 遍历坐标列表并将高程设置为 0
                    for i in range(len(coordinates[0])):
                        if len(coordinates[0][i]) >= 3:
                            coordinates[0][i][2] = 0
                elif geometry['type']=="LineString":
                    for i in range(len(coordinates)):
                        if len(coordinates[i]) >= 3:
                            coordinates[i][2] = 0
                elif geometry['type']=="Point":
                    coordinates[2] = 0

        # 将修改后的数据保存回 GeoJSON 文件
        with open(all_semantic_json_path.replace('.geojson','_noH.geojson'), 'w') as output_file:
            json.dump(geojson_data, output_file, indent=2)
        print("已将语义高程设置为 0，并保存到 _noH.geojson 文件中")

        '''
def get_cartographic_type(layer_name):

    layer_type_table = {'HADLane':1, 'HADLaneNode':2, 'HADRoadDivider':3, 'HADLaneDivider':4,
                        'LandMark':5, 'LocSign1':6, 'LocSign2':7, 'LocTrafficLight':8}

    if layer_name in layer_type_table:
        return layer_type_table[layer_name]
    else:
        return 0

def convert_cartographic_type_in_feature(feature, layer_name):
    #获取制图结果类别
    layer_type = get_cartographic_type(layer_name)

    if layer_name == 'HADLaneDivider' or layer_name == 'HADLaneNode'\
        or layer_name == 'HADRoadDivider' or layer_name == 'LandMark':
        if 'type' in feature['properties']:
            feature['properties']['subtype'] = feature['properties']['type']
        feature['properties']['type'] = layer_type
    elif layer_name == 'LocSign1' or layer_name == 'LocSign2':
        if 'subtype' in feature['properties']:
            feature['properties']['additional_linetype'] = feature['properties']['subtype']
        if 'type' in feature['properties']:
            feature['properties']['subtype'] = feature['properties']['type']
        feature['properties']['type'] = layer_type

def convert_cartographic_type(file_path, layer_name):
    with open(file_path, 'r', encoding='utf-8') as handle:
        fc = json.load(handle, parse_float=decimal.Decimal())

    if layer_name != 'HADLaneNode':
        for feature in fc['features']:
            convert_cartographic_type_in_feature(feature, layer_name)
    else:
        fs = []
        for feature in fc['features']:
            node_type = feature['properties']['type']
            if node_type == '4' or node_type == '5':
                fs.append(feature)
        fc = geojson.FeatureCollection(fs)

    with open(file_path, 'w', encoding='utf-8') as handle:
        handle.write(json.dumps(fc))

class ComplierProtoConvert(object):
    """
        编译数据proto转geojson类
    """
    def __init__(self):
        pass

    def convert(self, d_file_path, base_path = None):
        tile_data_header_struct = np.dtype([
            ("md5", (bytes, 16)),
            ("z_lib_level", "<u1"),
            ("main_proto_version", "<u1"),
            ("sub_proto_version", "<u1"),
            ("reserved", "<u1"),
            ("version", "<u4"),
            ("uncompressed_length", "<u4"),
            ("compressed_length", "<u4"),
        ])

        tile_data_dict = {}
        with open(d_file_path, "rb") as fp:
            header = np.fromfile(fp, dtype=tile_data_header_struct, count=1)[0]
            tile_data = tile_data_pb.TileData()
            unziped_tile_data = zlib.decompressobj().decompress(fp.read())
            tile_data.ParseFromString(unziped_tile_data)
            tile_data_dict = protobuf_to_dict(tile_data)

            # tile_data_dict = protobuf_to_dict(tile_data, including_default_value_fields=True)

        tile_time = tile_data_dict['TileTime']
        if not base_path:
            base_path = os.getcwd()

        result_root = os.path.join(base_path, str(tile_time['id']))
        if not os.path.exists(result_root):
            os.mkdir(result_root)

        for item, value in tile_data_dict.items():
            if item == 'TileTime':
                continue

            features = []

            layer_name = item[:-4]
            result_path = os.path.join(result_root, layer_name + '.geojson')
            for n in value:
                geo = {}
                if 'geometry' in n:
                    compiler_geo = eval(str(n['geometry']))
                    geo = self.convert_compiler_geometry_to_geojson(layer_name, compiler_geo)
                feature = {
                    "type": "Feature",
                    "geometry": self.nds_to_wcg(geo) if geo else geo,
                    "properties": {}
                }

                for key in n:
                    if key != 'geometry':
                        try:
                            value = n[key]
                            if key == 'id':
                                feature['properties'][key] = str(value)
                            else:
                                feature['properties'][key] = value
                        except Exception as e:
                            print(e)

                # self.write_type_in_feature(feature, layer_name)
                features.append(feature)

            with open(result_path, 'w', encoding='utf-8') as handle:
                data = {
                    "type": "FeatureCollection",
                    "features": features
                }
                handle.write(json.dumps(data, indent=4))

    def convert_compiler_geometry_to_geojson(
        self,
        layer_name,
        compiler_geo):

        geometry = None
        if (layer_name == 'HADLane' or
            layer_name == 'HADLaneDivider' or
            layer_name == 'HADRoadDivider'):
            geo_line = []
            for point in compiler_geo:
                geo_line.append(list(point.values()))
            geometry = shapelyGeo.LineString(geo_line)
        elif (layer_name == 'HADLaneAdas'or
                layer_name == 'HADLaneNode' or
                layer_name == 'LocSign' or
                layer_name == 'LocTrafficLight' or
                layer_name == 'HADLaneZlevelNode'):
            geometry = shapelyGeo.Point(list(compiler_geo.values()))
        elif (layer_name == 'LandMark' or
              layer_name == 'LocObj'):
            geo_coords = []
            for point in compiler_geo:
                geo_coords.append(tuple(point.values()))
            geometry = shapelyGeo.Polygon(tuple(geo_coords))

        if geometry is not None:
            geometry = eval(geojson.dumps(geometry))
        return geometry

    def nds_to_wcg(self, geometry):
        result_dict = {}
        if geometry.get("type") == "Polygon":
            point_polygon_list = []
            for coords in geometry.get("coordinates"):
                coords_list = []
                for coord in coords:
                    coords_list.append(
                        [coord[0] * 360 / 4294967296, coord[1] * 360 / 4294967296, 0.0])
                point_polygon_list.append((coords_list))
            result_dict["type"] = "Polygon"
            result_dict["coordinates"] = point_polygon_list
        if geometry.get("type") == "LineString":
            coords_list = []
            for coord in geometry.get("coordinates"):
                coords_list.append(
                    [coord[0] * 360 / 4294967296, coord[1] * 360 / 4294967296, 0.0])
            result_dict["type"] = "LineString"
            result_dict["coordinates"] = coords_list
        if geometry.get("type") == "Point":
            result_dict["type"] = "Point"
            coord = geometry.get("coordinates")
            result_dict["coordinates"] = [coord[0] * 360 / 4294967296, coord[1] * 360 / 4294967296, 0.0]
        return result_dict

class TencentProtoConvert(object):
    """
        编译数据proto转geojson类
    """
    def __init__(self):
        pass

    def convert(self, d_file_path, base_path = None):

        tile_data_header_struct = np.dtype([
            ("md5", (bytes, 16)),
            ("z_lib_level", "<u1"),
            ("main_proto_version", "<u1"),
            ("sub_proto_version", "<u1"),
            ("reserved", "<u1"),
            ("version", "<u4"),
            ("uncompressed_length", "<u4"),
            ("compressed_length", "<u4"),
        ])

        tile_data_dict = {}
        with open(d_file_path, "rb") as fp:
            header = np.fromfile(fp, dtype=tile_data_header_struct, count=1)[0]
            unziped_tile_data = zlib.decompressobj().decompress(fp.read())
            tile_data = tencent_tile_data_pb.TileData()
            tile_data.ParseFromString(unziped_tile_data)
            tile_data_dict = protobuf_to_dict(tile_data, including_default_value_fields=True)

        tile_time = tile_data_dict['TileTime']
        if not base_path:
            base_path = os.getcwd()

        result_root = os.path.join(base_path, str(tile_time['id']))
        if not os.path.exists(result_root):
            os.mkdir(result_root)

        for item, value in tile_data_dict.items():
            if item == 'TileTime':
                continue

            features = []
            layer_name = item[:-4]
            result_path = os.path.join(result_root, layer_name + '.geojson')
            for n in value:
                geo = {}
                if 'geometry' in n:
                    compiler_geo = eval(str(n['geometry']))
                    geo = self.convert_compiler_geometry_to_geojson(layer_name, compiler_geo)
                feature = {
                    "type": "Feature",
                    "geometry": geo,
                    "properties": {}
                }

                for key in n:
                    if key != 'geometry':
                        try:
                            value = n[key]
                            if key == 'id':
                                feature['properties'][key] = str(value)
                            else:
                                feature['properties'][key] = value
                        except Exception as e:
                            print(e)

                # self.write_type_in_feature(feature, layer_name)
                features.append(feature)

            with open(result_path, 'w', encoding='utf-8') as handle:
                data = {
                    "type": "FeatureCollection",
                    "features": features
                }
                handle.write(json.dumps(data, indent=4))

    def convert_compiler_geometry_to_geojson(
        self,
        layer_name,
        compiler_geo):

        geometry = None
        if (layer_name == 'HADLane' or
            layer_name == 'HADLaneDivider' or
            layer_name == 'HADRoadDivider' or
            layer_name == 'HADLink'):
            geo_line = []
            for point in compiler_geo:
                values = list(point.values())
                x, y = utils.nds_to_wgs84_coordinate(values[0]), utils.nds_to_wgs84_coordinate(values[1])
                geo_line.append([x, y, 0.0])
                # geo_line.append(list(point.values()))
            geometry = shapelyGeo.LineString(geo_line)
        elif (layer_name == 'HADLaneAdas'or
                layer_name == 'HADLaneNode' or
                layer_name == 'HADLaneShapePoint' or
                layer_name == 'LocSign' or
                layer_name == 'LocTrafficLight' or
                layer_name == 'HADLaneZlevelNode' or
                layer_name == 'HADNode' or
                layer_name == 'LocObj'):
            values = list(compiler_geo.values())
            x, y = utils.nds_to_wgs84_coordinate(values[0]), utils.nds_to_wgs84_coordinate(values[1])
            geometry = shapelyGeo.Point([x, y, 0.0])
            # geometry = shapelyGeo.Point(list(compiler_geo.values()))
        elif (layer_name == 'LandMark'):
            geo_coords = []
            for point in compiler_geo:
                x, y = utils.nds_to_wgs84_coordinate(point['x']), utils.nds_to_wgs84_coordinate(point['y'])
                geo_coords.append((x, y))
            if len(geo_coords) > 2:
                geometry = shapelyGeo.Polygon(tuple(geo_coords))
            else:
                geometry = shapelyGeo.LineString(tuple(geo_coords))

        if geometry is not None:
            geometry = eval(geojson.dumps(geometry))
        return geometry

if __name__ == '__main__':
    root = r'D:\data\records'
    all_vecs = os.listdir(root)

    for target_vec in all_vecs:
        # data_path = os.path.join(root,'data','samples',矢量数据)
        data_path = os.path.join(root,target_vec)
        current_car_data_dir_list = os.listdir(data_path) 
        for dir in current_car_data_dir_list:
            target_path = os.path.join(data_path,dir)
            download_path = target_path
            print(download_path)
            converter = CarProtoConvert()
            # all_lst = converter.read_all_location(r'C:\Users\202207817\Downloads\车端路网\trajectory')
            print('当前处理目录：'+target_path)
            converter.convert_all_data(target_path, download_path)
