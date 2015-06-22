import shapely.geometry as geom
import matplotlib.pyplot as plt
from copy import copy
from numpy import inf
import math
import munkres

class PolygonInterpolator:
  def __init__(self, p1, p2):
    self.p1 = p1
    self.p2 = p2
    self.compute_interpolation()

  def compute_interpolation(self):
    done = set([])
    pstrt = [geom.Point(p) for p in self.p1.exterior.coords[:-1]]
    pdest = [geom.Point(p) for p in self.p2.exterior.coords[:-1]]

    self.pairs = []

    if len(pdest) > len(pstrt):
      pstrt += midpoints(pstrt)

    #Use hungarian algorithm to find the best set of match points
    mat = [[p1.distance(p2) for p2 in pstrt] for p1 in pdest]
    m = munkres.Munkres()
    indexes = m.compute(mat)

    for i in indexes:
      closest = pstrt[i[1]]
      p = pdest[i[0]]
      done.add(closest)
      self.pairs.append((closest, p))

    #Match remaining start points to closest destination point
    for p in done.symmetric_difference(pstrt):
      closest, _ = project_point_points(p, pdest, 0)
      self.pairs.append((p, closest))

  def interpolate(self, percent):
    poly = []
    for p1, p2 in self.pairs:
      line = geom.LineString((as_tuple(p1), as_tuple(p2)))
      interp = line.interpolate(percent, normalized=True)
      poly.append(as_tuple(interp))

    return geom.Polygon(poly).convex_hull

def plot_polygons(polys):
  for p in polys:
    x, y = p.exterior.xy
    plt.plot(x, y)
  plt.show()

def plot_poly_normals(poly):
  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                       xlim=(-0.5, 2.5), ylim=(-0.5, 2.5))
  x, y = poly.exterior.xy
  ax.plot(x, y)
  n, _ = normals_offset(poly)
  mid = midpoints([geom.Point(p) for p in zip(*poly.exterior.coords.xy)[:-1]])
  print len(n), len(mid)
  x, y = zip(*[(p.xy[0], p.xy[1]) for p in mid])
  u, v = zip(*n)
  ax.quiver(x, y, u, v)

def as_tuple(point):
  return (point.xy[0][0], point.xy[1][0])

def midpoints(points):
  l = []
  for i, point in enumerate(points):
    prev = points[i-1]
    mid = geom.Point((point.xy[0][0]+prev.xy[0][0])/2,
                     (point.xy[1][0]+prev.xy[1][0])/2)
    l.append(mid)
  return l

#Deprecated, better to use hungarian algorithm
def interpolate_point_points(point, points, percent, exclude=set()):
    dists = [point.distance(p) if p not in exclude else inf for p in points]
    dmin = min(dists)
    closest = points[dists.index(dmin)]
    if point.xy == closest.xy:
      return closest, copy(closest)
    else:
      line = geom.LineString((as_tuple(point), as_tuple(closest)))
      interp = line.interpolate(percent, normalized=True)
      return closest, interp

def project_point_points(point, points, percent):
    dists = [point.distance(p) for p in points]
    dmin = min(dists)
    imin = dists.index(dmin)
    closest = points[imin]
    if imin < len(points) - 1:
      next_p = points[imin+1]
    else:
      next_p = points[-1]

    prev_p = points[imin-1]

    if point.distance(prev_p) > point.distance(next_p):
      edge = geom.LineString([as_tuple(closest), as_tuple(next_p)])
    else:
      edge = geom.LineString([as_tuple(prev_p), as_tuple(closest)])

    if point.xy == closest.xy:
      return closest, copy(closest)
    else:
      dist = edge.project(point)
      projection = edge.interpolate(dist)
      line = geom.LineString((as_tuple(point), as_tuple(projection)))
      interp = line.interpolate(percent, normalized=True)
      return projection, interp

#Deprecated : better to use Interpolator and hungarian algorithm
def interpolate_poly(ps, pd, percent=0.1):
  new_p = []
  done = set([])
  pstrt = [geom.Point(p) for p in ps.exterior.coords[:-1]]
  pdest = [geom.Point(p) for p in pd.exterior.coords[:-1]]

  if len(pdest) > len(pstrt):
    pstrt += midpoints(pstrt)

  #Match every destination point to closest start point
  for p in pdest:
   closest, interpolate = interpolate_point_points(p, pstrt, 1-percent, done)
   new_p.append(as_tuple(interpolate))
   done.add(closest)

  #Match remaining start points to closest destination point
  for p in done.symmetric_difference(pstrt):
    projection, interpolate = project_point_points(p, pdest, percent)
    new_p.append(as_tuple(interpolate))

  return geom.Polygon(new_p).convex_hull

def normals_offset(polygon):
  n = []
  o = []
  coords = polygon.exterior.coords[:-1]
  for i, p in enumerate(coords):
    normal = (p[1] - coords[i-1][1], -(p[0] - coords[i-1][0]))
    norm = math.sqrt(normal[0]**2+normal[1]**2)
    normal = (normal[0]/norm, normal[1]/norm)
    n.append(normal)
    o.append(-(normal[0]*p[0] + normal[1]*p[1]))
  return n, o
