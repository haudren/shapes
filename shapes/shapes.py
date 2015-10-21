import shapely.geometry as geom
import matplotlib.pyplot as plt
from copy import copy
from numpy import inf, dot
import math
import munkres
import scipy.spatial as spatial
import numpy as np

class PolygonInterpolator:
  def __init__(self, p1, p2):
    self.p1 = p1
    self.p2 = p2

    self.compute_interpolation()
    self.compute_vertex_order()

  def compute_interpolation(self):
    done = set([])

    pstrt = [geom.Point(p) for p in self.p1.exterior.coords[:-1]]
    pdest = [geom.Point(p) for p in self.p2.exterior.coords[:-1]]

    self.pairs, self.tuple_pairs = [], []

    while len(pdest) > len(pstrt):
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
      self.tuple_pairs.append((as_tuple(closest), as_tuple(p)))

    #Match remaining start points to closest destination point
    for p in done.symmetric_difference(pstrt):
      closest, _ = project_point_points(p, pdest, 0)
      self.pairs.append((p, closest))
      self.tuple_pairs.append((as_tuple(p), as_tuple(closest)))

  def compute_vertex_order(self):
    points = np.vstack(self.fast_interpolate_pairs(0.5))
    hull = spatial.ConvexHull(points)
    #bar = (sum([p[0][0] for p in self.tuple_pairs]),
    #       sum([p[0][1] for p in self.tuple_pairs]))

    #order_ccw = lambda p: math.atan2(p[0][1] - bar[1], p[0][0] - bar[0])
    #order_ccw_p = lambda p: math.atan2(p[0].xy[1][0] - bar[1], p[0].xy[0][0] - bar[0])

    #self.tuple_pairs = sorted(self.tuple_pairs, key=order_ccw)
    #self.pairs = sorted(self.pairs, key=order_ccw_p)
    self.order = hull.vertices
    self.tuple_pairs = [self.tuple_pairs[i] for i in self.order]
    self.pairs = [self.pairs[i] for i in self.order]

  def interpolate(self, percent):
    poly = []
    for p1, p2 in self.pairs:
      line = geom.LineString((as_tuple(p1), as_tuple(p2)))
      interp = line.interpolate(percent, normalized=True)
      poly.append(as_tuple(interp))

    return geom.Polygon(poly).convex_hull

  def fast_interpolate_pairs(self, percent):
    perc = max(min(percent, 1.), -1.)
    if perc < 0:
      perc = 1 + perc
    return [((t1[0]*(1-perc) + t2[0]*perc,
              t1[1]*(1-perc) + t2[1]*perc))
            for t1, t2 in self.tuple_pairs]

  def fast_interpolate(self, percent):
    pairs = self.fast_interpolate_pairs(percent)
    return geom.Polygon(pairs).convex_hull

  def point_derivative(self, epsilon_derivative):
    return [point_derivative(p1, p2, epsilon_derivative)
            for p1, p2 in self.pairs]

  def pairs_derivative(self, epsilon_derivative):
    return [tuple_derivative(p1, p2, epsilon_derivative)
            for p1, p2 in self.tuple_pairs]

  def midpoint_derivative(self, epsilon_derivative):
    ps, pd = zip(*self.pairs)
    m_strt = midpoints(ps)
    m_dest = midpoints(pd)
    return [point_derivative(m, n, epsilon_derivative)
            for m, n in zip(m_strt, m_dest)]

  def point_dist_derivative(self, point, percent, epsilon_derivative):
    spd_points = self.pairs_derivative(epsilon_derivative)
    #cur_poly = self.fast_interpolate(percent).exterior.coords
    cur_poly = self.fast_interpolate_pairs(percent)
    spd_segment = [None]*len(spd_points)

    for i, p in enumerate(cur_poly):
      p2 = cur_poly[i-1]
      vec = (p[0] - p2[0], p[1] - p2[1])
      vec_point = (point[0] - p2[0], point[1] - p2[1])
      length = dot(vec, vec)
      dist = dot(vec, vec_point)

      spd_x = (spd_points[i-1][0]*dist+spd_points[i][0]*(length-dist))/length
      spd_y = (spd_points[i-1][1]*dist+spd_points[i][1]*(length-dist))/length
      #spd_y = spd_points[i-1][1]+spd_points[i][1]*dist/length

      #print vec
      #print vec_point
      #print spd_points[i-1], spd_points[i]
      #print length, dist
      #print spd_x, spd_y
      #print '---'
      spd_segment[i] = (spd_x, spd_y)

    return spd_segment

  def projections(self, point, percent):
    #cur_poly = self.fast_interpolate(percent)
    cur_poly = self.fast_interpolate_pairs(percent)
    #proj = [None]*(len(cur_poly.exterior.coords)-1)
    proj = [None]*len(cur_poly)

    for i, p in enumerate(cur_poly):
      p_prec = cur_poly[i-1]
      vec = (p[0] - p_prec[0],
             p[1] - p_prec[1])
      vec_point = (point[0] - p_prec[0],
                   point[1] - p_prec[1])

      dist = dot(vec_point, vec)
      length = dot(vec, vec)

      x = p_prec[0]+(dist*vec[0]/length)
      y = p_prec[1]+(dist*vec[1]/length)

      proj[i] = (x, y)

    return proj

  def normals_offset(self, percent):
    points = self.fast_interpolate_pairs(percent)
    normals = [None]*len(points)
    offsets = [None]*len(points)
    norms = [None]*len(points)

    for i, p in enumerate(points):
      p2 = points[i-1]
      normal = (-(p[1] - p2[1]), p[0] - p2[0])
      norm = math.sqrt(normal[0]**2+normal[1]**2)
      norms[i] = norm
      if norm > 9e-3:
        normal = (normal[0]/norm, normal[1]/norm)
        offset = -(normal[0]*p[0] + normal[1]*p[1])
      else:
        normal = (0.0, 0.0)
        offset = 0
      normals[i] = normal
      offsets[i] = offset
    print norms
    return normals, offsets

  def perc_normal_derivative(self, epsilon_derivative):
    n_strt, _ = self.normals_offset(0.)
    n_dest, _ = self.normals_offset(1.)
    return [tuple_derivative(n, m, epsilon_derivative)
            if n != (0., 0.) and n_dest != (0., 0.)
            else (0., 0.)
            for n, m in zip(n_strt, n_dest)]

  def normal_derivative(self, epsilon_derivative):
    ps, pd = zip(*self.pairs)
    poly_s = geom.Polygon([p.coords[0] for p in ps])
    poly_d = geom.Polygon([p.coords[0] for p in pd])
    n_strt, _ = normals_offset(poly_s)
    n_dest, _ = normals_offset(poly_d)
    return [tuple_derivative(n, m, epsilon_derivative)
            if n != (0., 0.) and n_dest != (0., 0.)
            else (0., 0.)
            for n, m in zip(n_strt, n_dest)]

def point_derivative(p1, p2, e_n):
  return ((p2.x-p1.x)/e_n, (p2.y-p1.y)/e_n)

def tuple_derivative(p1, p2, e_n):
  return ((p2[0]-p1[0])/e_n, (p2[1]-p1[1])/e_n)

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
    if norm > 0:
      normal = (normal[0]/norm, normal[1]/norm)
    else:
      normal = (0.0, 0.0)
    n.append(normal)
    o.append(-(normal[0]*p[0] + normal[1]*p[1]))
  return n, o
