import requests
import geopandas as gpd
import numpy as np
from requests import ReadTimeout
from shapely.geometry import Polygon, Point, LineString, MultiPolygon, MultiLineString, MultiPoint, shape
import pandas as pd
import folium
import json
import random
import osmnx as ox
from osmnx._errors import *
from datetime import datetime


def load_park(id):
    data = requests.get(f'http://194.113.35.2:8000/api/v1/territories/{id}').json()['data']
    if data['type'] == 'Polygon':
        park_polygon = Polygon(data['coordinates'])
    elif data['type'] == 'MultiPolygon':
        park_polygon = MultiPolygon(data['coordinates'])
    return park_polygon


def find_three_equidistant_points(polygon):
    """
    Находит три точки на полигоне, расположенные примерно на равном расстоянии.
    """
    boundary = polygon.boundary
    total_length = boundary.length
    segment_length = total_length / 3
    points = []
    distance_covered = 0
    for i in range(1, 4):
        distance_covered += segment_length
        point = boundary.interpolate(distance_covered)
        points.append(point)
    return points


def create_polygons_around_points(coord, sizes, penalty=None, grid=None):
    """
    Создает квадраты вокруг точек заданных размеров или расширяет полигоны.
    """
    polygons = {size: [] for size in sizes}
    coord = prepare_for_draw(coord)
    for size in sizes:
        polygon = coord.buffer(size)
        polygons[size].append(polygon)
    gdf_small = gpd.GeoDataFrame(geometry=polygons[sizes[0]])
    gdf_big = gpd.GeoDataFrame(geometry=polygons[sizes[1]])
    gdf_small.set_crs(epsg=4326, inplace=True)
    gdf_big.set_crs(epsg=4326, inplace=True)
    for ter in gdf_small.geometry:
        intersecting_cellsf = grid[grid.geometry.intersects(ter)]
        grid.loc[intersecting_cellsf.index, 'index_value'] -= penalty
    for ter in gdf_big.geometry:
        intersecting_cellss = grid[(grid.geometry.intersects(ter))]
        intersecting_cellss = intersecting_cellss[~intersecting_cellss.index.isin(intersecting_cellsf.index)]
        grid.loc[intersecting_cellss.index, 'index_value'] -= penalty / 2


def ecological_function(beta, beta_k, alpha_k_lower, alpha_k_upper, zero_optimum, zero_optimum_upper):
    sum_term = 0
    zero_optimum_term = 0
    for i in range(len(beta)):
        if beta[i] >= beta_k[i]:
            sum_term += (beta[i] - beta_k[i]) ** 2 / (alpha_k_upper[i])
        else:
            sum_term += (beta[i] - beta_k[i]) ** 2 / (alpha_k_lower[i])
        print(f"sum term   {sum_term}")
    for i in range(len(zero_optimum)):
        zero_optimum_term += zero_optimum[i] ** 2 / (zero_optimum_upper[i])
    return 90 * np.exp(-sum_term) + 10 * np.exp(-zero_optimum_term)


# Функция для создания сетки
def create_grid(bounds, cell_size):
    xmin, ymin, xmax, ymax = bounds
    width = xmax - xmin
    height = ymax - ymin
    rows = int(np.ceil(height / cell_size))
    cols = int(np.ceil(width / cell_size))

    x = np.arange(xmin, xmax, cell_size)
    y = np.arange(ymin, ymax, cell_size)

    grid_cells = [Polygon([(x_left, y_top), (x_left + cell_size, y_top),
                           (x_left + cell_size, y_top - cell_size), (x_left, y_top - cell_size)])
                  for x_left in x for y_top in y[::-1]]

    grid = gpd.GeoDataFrame(grid_cells, columns=['geometry'])
    return grid


def prepare_for_draw(coord):
    if isinstance(coord, Point) and coord.has_z:
        x, y, _ = coord.coords[0]
        return Point(x, y)
    elif isinstance(coord, MultiPoint):
        new_points = []
        for point in coord.geoms:
            if point.has_z:
                x, y, _ = point.coords[0]
                new_points.append(Point(x, y))
            else:
                new_points.append(point)
        return MultiPoint(new_points)
    elif isinstance(coord, LineString) and coord.has_z:
        new_coords = [(x, y) for x, y, z in coord.coords]
        return LineString(new_coords)
    elif isinstance(coord, MultiLineString):
        new_linestrings = []
        for linestring in coord.geoms:
            if linestring.has_z:
                new_coords = [(x, y) for x, y, z in linestring.coords]
                new_linestrings.append(LineString(new_coords))
            else:
                new_linestrings.append(linestring)
        return MultiLineString(new_linestrings)
    elif isinstance(coord, Polygon) and coord.has_z:
        new_coords = [(x, y) for x, y, z in coord.exterior.coords]
        return Polygon(new_coords, [list(map(lambda p: (p[0], p[1]), interior.coords)) for interior in coord.interiors])
    elif isinstance(coord, MultiPolygon):
        new_polygons = []
        for polygon in coord.geoms:
            if polygon.has_z:
                new_coords = [(x, y) for x, y, z in polygon.exterior.coords]
                new_polygons.append(Polygon(new_coords,
                                            [list(map(lambda p: (p[0], p[1]), interior.coords)) for interior in
                                             polygon.interiors]))
            else:
                new_polygons.append(polygon)
        return MultiPolygon(new_polygons)
    return coord


def ensure_2d_coordinates(geometry):  # -> Any:
    """Convert all coordinates in a GeoJSON geometry to 2D by removing any Z values."""
    coords = geometry['coordinates']

    def convert_to_2d(coords):
        """Recursively convert coordinates to 2D."""
        if isinstance(coords[0], (list, tuple)):
            return [convert_to_2d(c) for c in coords]
        else:
            return coords[:2]

    geometry['coordinates'] = convert_to_2d(coords)
    return geometry


def json_to_gdf(json_data):
    features = json_data['features']
    for feature in features:
        feature['geometry'] = ensure_2d_coordinates(feature['geometry'])
    geometries = [shape(feature['geometry']) for feature in features]
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf.set_crs(epsg=4326, inplace=True)
    return gdf


def trail(data):
    gdf = json_to_gdf(data)
    for i in range(len(gdf.geometry)):
        pol = prepare_for_draw(gdf.geometry[i])
        create_polygons_around_points(pol, [0.001, 0.002], 20 * min(anthropogenic_pressure, 1), grid)
        create_polygons_around_points(pol, [0.001, 0.002], 3 * (sum(marks) / len(marks)), grid)


def get_applications(id):
    c = 0
    data = requests.get(f'http://194.113.35.2:8000/api/v1/applications/')
    if data.status_code == 200:
        data = data.json()
        for el in data:
            if el['territory_id'] == id and datetime.strptime(el['start_date'],
                                                              "%Y-%m-%d") < datetime.now() < datetime.strptime(
                    el['end_date'], "%Y-%m-%d"):
                c += len(el['visitors'])
        return c
    else:
        return 0


def get_feedbacks(id):
    feedback = 0
    data = requests.get(f'http://194.113.35.2:8000/api/v1/feedbacks/by-track/{id}')
    if data.status_code == 200:
        data = data.json()
        i = 0
        for el in data:
            if i < 100:
                feedback += (el["rate_painting"] + el["rate_facilities"] + el["rate_purity"] + el[
                    "rate_expectations"]) / 4
            i += 1
        return feedback / i
    else:
        return 0


def get_trails(territory_id):
    dic = {1: 1, 2: 4, 3: 2, 4: 3}
    trails = []
    ids = []
    data = requests.get(f'http://194.113.35.2:8000/api/v1/territories/with_tracks/{dic[territory_id]}').json()['tracks']
    for el in data:
        ids.append(el['id'])
        trail = requests.get(f"http://194.113.35.2:8000/api/v1/tracks/{el['id']}").json()
        trails.append(json_to_gdf(trail['data']))
    return trails, ids


def get_current_weather(lat, lon, api_key='3d5869b98df729ba78f03f53ace903e7'):
    call_url = f'https://api.openweathermap.org/data/2.5/weather?lat={lat}&lon={lon}&appid={api_key}'
    call_data = requests.get(call_url).json()
    return_data = dict()
    return_data['temp'] = call_data['main']['temp']
    return_data['humidity'] = call_data['main']['humidity']
    if 'sea_level' in call_data['main'].keys():
        return_data['sea_level'] = call_data['main']['sea_level']
    else:
        return_data['sea_level'] = call_data['main']['pressure']
    if 'rain' in call_data.keys():
        if '1h' in call_data['rain'].keys():
            return_data['rain_1h'] = call_data['rain']['1h']
        if '3h' in call_data['rain'].keys():
            return_data['rain_3h'] = call_data['rain']['3h']

    if 'snow' in call_data.keys():
        if '1h' in call_data['snow'].keys():
            return_data['snow_1h'] = call_data['snow']['1h']
        if '3h' in call_data['snow'].keys():
            return_data['snow_3h'] = call_data['snow']['3h']
    return_data['day_length'] = call_data['sys']['sunset'] - call_data['sys']['sunrise']
    # return_data['wind_gust'] = call_data['wind']['gust']
    return_data['wind_speed'] = call_data['wind']['speed']
    return_data['counter'] = 1
    for el in ['rain_1h', 'rain_3h', 'snow_1h', 'snow_3h']:
        if not el in list(return_data.keys()):
            return_data[el] = 0
    return return_data


def get_current_air_pollution(lat, lon, api_key='3d5869b98df729ba78f03f53ace903e7'):
    call_url = f'http://api.openweathermap.org/data/2.5/air_pollution?lat={lat}&lon={lon}&appid={api_key}'
    call_data = requests.get(call_url).json()

    return_data = dict()
    return_data['co'] = call_data['list'][0]['components']['co']
    return_data['nh3'] = call_data['list'][0]['components']['nh3']
    return_data['no'] = call_data['list'][0]['components']['no']
    return_data['no2'] = call_data['list'][0]['components']['no2']
    return_data['o3'] = call_data['list'][0]['components']['o3']
    return_data['pm10'] = call_data['list'][0]['components']['pm10']
    return_data['pm2_5'] = call_data['list'][0]['components']['pm2_5']
    return_data['so2'] = call_data['list'][0]['components']['so2']
    return_data['counter'] = 1

    return return_data


koef = 111320
territory_id = 1
park_polygon = load_park(territory_id)
gdf_park = gpd.GeoDataFrame(index=[0], geometry=[park_polygon])
gdf_park.set_crs(epsg=4326, inplace=True)

points = find_three_equidistant_points(park_polygon)
weather = {'temp': 0,
           'humidity': 0,
           'sea_level': 0,
           'rain_1h': 0,
           'rain_3h': 0,
           'snow_1h': 0,
           'snow_3h': 0,
           'day_length': 0,
           'wind_speed': 0,
           'counter': 0}
air_p = {'co': 0,
         'nh3': 0,
         'no': 0,
         'no2': 0,
         'o3': 0,
         'pm10': 0,
         'pm2_5': 0,
         'so2': 0,
         'counter': 0}
try:
    for point in points:
        cur_w = get_current_weather(point.y, point.x)
        cur_ap = get_current_air_pollution(point.y, point.x)
        for param in weather.keys():
            weather[param] += cur_w[param]
        for param in cur_ap.keys():
            air_p[param] += cur_ap[param]
    for k, _ in weather.items():
        weather[k] = round(weather[k] / weather['counter'], 2)
    for k, _ in air_p.items():
        air_p[k] = round(air_p[k] / air_p['counter'], 2)
except ReadTimeout:
    print('API broke down')
    weather = {'temp': 293, 'humidity': 50, 'sea_level': 1013, 'rain_1h': 4, 'rain_3h': 12, 'snow_1h': 0, 'snow_3h': 0,
               'day_length': 43200, 'wind_speed': 3, 'counter': 1}
    air_p = {'co': 0, 'nh3': 0, 'no': 0, 'no2': 0, 'o3': 0, 'pm10': 0, 'pm2_5': 0, 'so2': 0, 'counter': 1}
# 1. Оценка здоровья парка в целом
# Параметры
betas = [
    weather['temp'],  # температура
    weather['humidity'],  # влажность
    weather['sea_level'],  # давление
    weather['wind_speed'],  # скорость ветра
    weather['rain_1h'] + weather['snow_1h'],  # осадки час
    weather['rain_3h'] + weather['snow_3h'],  # осадки 3 часа
]
beta_ks = [  # оптимум экологических факторов
    293,
    50,
    1013,
    4,
    4,
    3,
]
upper_bound = [
    250,
    20000,
    15000,
    500,
    300,
    900,
]
lower_bound = [
    10000,
    23000,
    8000,
    400,
    1000,
    10000,
]
zero_optimum = [
    air_p['co'],  # CO
    air_p['nh3'],  # NH3
    air_p['no'],  # NO
    air_p['no2'],  # NO2
    air_p['o3'],  # O3
    air_p['pm10'],  # PM10
    air_p['pm2_5'],  # PMPM2.5
    air_p['so2']  # SO2
]
zero_optimum_upper_bound = [
    15400,
    200,
    100,
    200,
    180,
    200,
    75,
    350
]
park_health = ecological_function(betas, beta_ks, lower_bound, upper_bound, zero_optimum, zero_optimum_upper_bound)
print(f"Общая оценка здоровья парка: {park_health:.2f}")

# 2. Деление парка на сетку 100x100 метров
cell_size = 0.001  # Примерный размер клетки 100x100 метров в градусах
grid = create_grid(gdf_park.total_bounds, cell_size)
grid.set_crs(epsg=4326, inplace=True)

sindex = grid.sindex
possible_matches_index = list(sindex.intersection(park_polygon.bounds))
possible_matches = grid.iloc[possible_matches_index]
precise_matches = possible_matches[possible_matches.intersects(park_polygon)]
grid = precise_matches
grid['index_value'] = park_health

# 3. Штрафы для ячеек с антропогенными объектами
north, south, east, west = gdf_park.total_bounds[3], gdf_park.total_bounds[1], gdf_park.total_bounds[2], \
gdf_park.total_bounds[0]
tags_with_penalty_and_area = {
    ('waste', 'landfill'): (15, [200 / koef, 1000 / koef]),  # свалка
    ('waste', 'dump_site'): (8, [200 / koef, 1000 / koef]),  # место для свалки
    ('waste', 'incinerator'): (20, [300 / koef, 2000 / koef]),  # мусоросжигательный завод
    ('industrial', 'asbestos'): (10, [200 / koef, 1000 / koef]),  # завод и дальше уточнение какоц
    ('industrial', 'cement'): (15, [200 / koef, 1000 / koef]),
    ('industrial', 'chemical'): (10, [200 / koef, 1000 / koef]),
    ('industrial', 'coal'): (20, [300 / koef, 2000 / koef]),
    ('industrial', 'oil'): (20, [300 / koef, 1500 / koef]),
    ('industrial', 'steel'): (15, [200 / koef, 1000 / koef]),
    ('man_made', 'works'): (10, [100 / koef, 500 / koef]),  # завод или фабрика
    ('man_made', 'pipeline'): (6, [100 / koef, 200 / koef]),  # трубопровод
    ('man_made', 'mineshaft'): (7, [200 / koef, 500 / koef]),  # шахта
    ('landuse', 'quarry'): (15, [200 / koef, 1000 / koef]),  # карьер
    ('landuse', 'landfill'): (25, [100 / koef, 500 / koef]),  # полигон для захоронения отходов
    ('landuse', 'industrial'): (7, [200 / koef, 1000 / koef]),
    ('power', 'plant'): (5, [200 / koef, 1000 / koef])  # электростанция
}
for (key, value), (penalty, sizes) in tags_with_penalty_and_area.items():
    try:
        gdf = ox.features_from_bbox(bbox=[north, south, east, west], tags={key: value})
        if gdf.empty:
            continue
        else:
            for el in gdf.geometry:
                create_polygons_around_points(el, sizes, penalty, grid)
    except InsufficientResponseError as e:
        print(f"Ошибка при загрузке данных для тега {key}={value}: {e}")

# # 4. Дополнительно для троп

trails = {}
tracks, ids = get_trails(territory_id)
intersecting_cells1 = set()
intersecting_cells2 = set()

for i in range(len(tracks)):
    ic1 = set()
    ic2 = set()
    print(ids[i])
    tr_id = ids[i]
    if tr_id in [27, 28]:
        continue
    marks = get_feedbacks(tr_id)
    ppl = get_applications(tr_id)
    anthropogenic_pressure = 0.5  # ppl / (10 * (1 + 2 * (weather['day_length'] / 3600 - 10)))
    gdf = tracks[i]
    combined_geom = gdf.geometry.unary_union
    sizes = [0.001, 0.002]
    polygons = {size: combined_geom.buffer(size) for size in sizes}
    for size, poly in polygons.items():
        gdf_buffered = gpd.GeoDataFrame(geometry=[poly], crs='epsg:4326')
        sindex = grid.sindex
        possible_matches_index = list(sindex.intersection(poly.bounds))
        possible_matches = grid.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(poly)]
        if size == 0.001:
            cells1 = precise_matches
            ic1.update(cells1.index)
        else:
            cells2 = precise_matches[~precise_matches.index.isin(ic1)]
            ic2.update(cells2.index)
    grid.loc[list(ic1 - intersecting_cells1), 'index_value'] *= (1 - 0.2 * min(anthropogenic_pressure, 1))
    grid.loc[list(ic2 - intersecting_cells2), 'index_value'] *= (1 - (0.2 * min(anthropogenic_pressure / 2, 1)))
    grid.loc[list(ic1 - intersecting_cells1), 'index_value'] *= max(marks / 5 * 1.3, 0.7)
    grid.loc[list(ic2 - intersecting_cells2), 'index_value'] *= max(marks / 5 * 1.3 / 2, 0.7)
    intersecting_cells1.update(ic1)
    intersecting_cells2.update(ic2)
    trails[tr_id] = grid.loc[list(ic1 - intersecting_cells1), 'index_value'].min()

grid = grid.dissolve(by='index_value', as_index=False)
grid['date'] = datetime.today().strftime("%d.%m")

grid.to_file(f"park{territory_id}.geojson", driver="GeoJSON")

# создание карты
# import folium.plugins
# def add_geometry_to_map(geometry, map_obj, color='blue'):
#     if isinstance(geometry, LineString):
#         folium.PolyLine(locations=[(y, x) for x, y in geometry.coords], color=color).add_to(map_obj)
#     elif isinstance(geometry, MultiLineString):
#         for geom in geometry.geoms:
#             folium.PolyLine(locations=[(y, x) for x, y in geom.coords], color=color).add_to(map_obj)
#     elif isinstance(geometry, MultiPolygon):
#         for polygon in geometry.geoms:
#             for poly in [polygon.exterior] + list(polygon.interiors):
#                 folium.Polygon(locations=[(y, x) for x, y in poly.coords], color=color).add_to(map_obj)
#     elif isinstance(geometry, Polygon):
#         for poly in [geometry.exterior] + list(geometry.interiors):
#             folium.Polygon(locations=[(y, x) for x, y in poly.coords], color=color).add_to(map_obj)
#     elif isinstance(geometry, MultiPoint):
#         for point in geometry.geoms:
#             folium.Marker(location=[point.y, point.x]).add_to(map_obj)
#     elif isinstance(geometry, Point):
#         folium.Marker(location=[geometry.y, geometry.x]).add_to(map_obj)
#     else:
#         raise ValueError(f"Unsupported geometry type: {type(geometry)}")


# def get_color(index_value):
#     index_value = index_value
#     colors = ['black', 'gray', 'white']
#     return colors[0 if index_value < 30 else 1 if index_value < 68 else 2]
# # Функция стиля для раскрашивания ячеек
# def style_function(feature):
#     return {
#         'fillOpacity': 0.5,
#         'weight': 2,
#         'fillColor': get_color(feature['properties']['index_value'])
#     }
# def park_style_function(feature):
#     return {
#         'fillOpacity': 0.1,
#         'weight': 2,
#         'color': None,
#         'fillColor': None
#     }

# m = folium.Map(location=[52,150], zoom_start=5)
# m.add_child(folium.plugins.MeasureControl())
# # folium.GeoJson(gdf_park).add_to(m)
# minx, miny, maxx, maxy = gdf_park.total_bounds
# # folium.GeoJson(Polygon([[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]])).add_to(m)
# for el in tracks:
#     for geom in el.geometry:
#         add_geometry_to_map(geom, m, 'red')
# folium.GeoJson(grid, style_function=style_function).add_to(m)

# m.save('park_grid_map.html')
# m
