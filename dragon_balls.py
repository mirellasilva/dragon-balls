from scipy.spatial import distance_matrix
from scipy.optimize import minimize
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon
from math import radians, cos, sin, asin, sqrt
import fiona

def calculateDistancePair(lat1, lon1, lat2, lon2):
    # The math module contains a function named
    # radians which converts from degrees to radians.
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    
    # Haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2

    c = 2 * asin(sqrt(a)) 
    
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
    
    # calculate the result
    return(c * r)

def calculate_total_distance(coordinates):
  half_len = int(len(coordinates)/2)
  total_distance = 0
  for i in range(half_len):
    for j in range(i + 1, half_len):
      # Calculate distance between each pair of coordinates and add to total
      total_distance += calculateDistancePair(coordinates[i*2], coordinates[i*2+1], coordinates[j*2], coordinates[j*2+1])
  return total_distance

def isLand(gdf, latitude, longitude):
    #point = Point(-11.037601, -4.8417559)  # water
    point = Point(longitude, latitude)  # land
    # Create an empty list to store Polygon objects
    polygons = []
    point_contained = False

    # Iterate through each row (shape) in the GeoDataFrame
    for index, row in gdf.iterrows():
        # Access the geometry of the shape
        geometry = row['geometry']
        
        # Check if the geometry is a Polygon or MultiPolygon
        if isinstance(geometry, Polygon):
          #  print("is polygon") 
            # Create a Polygon object and append it to the list
            polygons.append(geometry)
            if geometry.contains(point):
              point_contained = True
        elif isinstance(geometry, MultiPolygon):
          #  print("is multi-polygon") 
            # If the geometry is a MultiPolygon, iterate through its constituent polygons
            for polygon in geometry.geoms:
                # Create a Polygon object for each constituent polygon and append it to the list
                polygons.append(polygon)
                if polygon.contains(point):
                  point_contained = True
        #else:
        #    print("neither") 


    if point_contained:
        #print("The point is contained in at least one polygon within the MultiPolygon.")
        # its land
        return True
    else:
        # print("The point is not contained in any polygon within the MultiPolygon.")
        # its water
        return False


def total_distance_land(positions):
    gdf = []
    try:
      gdf = gpd.read_file('ne_10m_land')
    except OSError as error:
      sys.exit(error)
      return
    # Reshape positions to separate the fixed first point
    fixed_point = [-12.3084, -55.452]  # Separate first point (longitude, latitude)
    merged_positions = fixed_point + positions.tolist()
    variable_positions = positions[2:].reshape(-1, 2)  # Reshape remaining points as (6, 2)
    #print("\npositions", positions, "\nfixed point", fixed_point, "\nvariable_positions:", variable_positions, "\nteste: ", merged_list)
    
    # Calculate distances using variable positions and the fixed first point
    # distance_matrix = distance_matrix(variable_positions, np.array([fixed_point]))
    distance_sum = calculate_total_distance(merged_positions)
    distance_sum_original = distance_sum
    all_land = gdf.unary_union
    distance = 0
    # Check for land and impose penalty only for variable positions
    for i, position in enumerate(variable_positions):
      point = Point(position[1], position[0])
      distance = point.distance(all_land)
      # if isLand(gdf, position[1], position[0]) == False:
      distance_sum += distance * (-200000)
      if distance > 0:
        break

    print("\ndistance_sum original e diferença: ", distance_sum_original, distance_sum - distance_sum_original)
    if distance_sum_original > 22000 and distance_sum - distance_sum_original > -1:
      print("\nsolution: ", variable_positions)
    return -distance_sum


# Define the number of balls
num_balls = 7

# Define initial guess (can be random or based on heuristics)
#initial_guess = [51.5553, -105.8189, 2.710676, 24.786559, 61.38768, 69.271557, -24.20024, 129.457487, 35.87418, 138.210087, -80.8930102, 62.33244603]
#initial_guess = [2.02173223, 24.38317631, 61.85332409, 70.09937284, -24.61042638, 130.27410599, 34.61215877, 138.85872568, -81.04121963, 62.33010985]
initial_guess = [-0.12716669704727976, 25.441026173546152] * 6

# You can replace this with your preferred initial guess generation method
# For example, placing points on a cube inscribed within the sphere

# Set bounds to restrict positions to the sphere's surface
bounds = [(-90, 90), (0, 180)] * (num_balls - 1)  # latitude, longitude limits

# Minimize the total distance function
result = minimize(total_distance_land, initial_guess, method='SLSQP', bounds=bounds)

# Extract the optimized positions
optimized_positions = result.x.reshape(num_balls-1, 2)  # Reshape to (longitude, latitude)
print(f"Ball origin: (-12.3084,-55.452)")

# Print the optimized positions (longitude, latitude)
for i in range(num_balls-1):
  print(f"Ball {i+1}: ({optimized_positions[i, 0]:.4f},{optimized_positions[i, 1]:.4f})")


fixed_point = [-12.3084, -55.452]
merged_positions = fixed_point + result
final_distance = calculate_total_distance(merged_positions)
print(f"final_distance: ", final_distance)