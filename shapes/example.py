from shapes import PolygonInterpolator, normals_offset, midpoints
import matplotlib.pyplot as plt
import shapely.geometry as geom
import numpy as np
import matplotlib.animation as animation

p1s = geom.Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
#p2s = geom.Polygon([(0, 0), (0, 1), (0.5, 1.5), (1, 1.5),
#                    (1.5, 1.5), (1.5, 1.), (1, 0)])
#p2s = geom.Polygon([(0, 0), (0, 1), (0.5, 1.5),
#                    (1.5, 1.), (1, 0)])
p2s = geom.Polygon([(0, 0), (0, 1), (0.5, 1.5), (1, 1), (1, 0)])
#
#lrpoly = geom.Polygon(zip([-0.07050964684556796, -0.07050964825832298, 0.10949034957097684, 0.10949034894607837, -0.07050964684556796], [-0.15499999638048872, 0.15499999203700351, 0.1549999968283623, -0.15499999600629497, -0.15499999638048872]))
#
#rpoly = geom.Polygon(zip([0.10949026646949628, -0.0705096371815371, -0.07050963781963734, 0.1094903420613558, 0.10949026646949628], [0.05500000293727045, 0.05500001448651814, 0.15499998281677235, 0.15499999497076253, 0.05500000293727045]))
#
#lpoly = geom.Polygon(zip([0.10949026646949628, -0.0705096371815371, -0.07050963781963734, 0.1094903420613558, 0.10949026646949628], [0.05500000293727045, 0.05500001448651814, 0.15499998281677235, 0.15499999497076253, 0.05500000293727045]))
#
#p1s = lpoly
#p2s = lrpoly

def save_polys(epsilons):
  inter = PolygonInterpolator(p1s, p2s)
  all_points = None

  for i, eps in enumerate(epsilons):
    poly = inter.fast_interpolate(eps)
    arr = np.array(poly.exterior.coords)
    np.savetxt('polygon_{}'.format(i), arr, header='x y', comments='')

    normals, offsets = normals_offset(poly)

    if eps > 0 and eps < 1:
      speeds = inter.midpoint_derivative(1.)
    else:
      speeds = [(0, 0)]*len(normals)

    #  if prev_offsets is None:
    #    speeds = [0]*len(normals)
    #  else:
    #    speeds = [n-p for n, p in zip(offsets, prev_offsets)]
    #  prev_offsets = copy(offsets)

    #dots = [n[0]*s[0]+n[1]*s[1] for n, s in zip(normals, speeds)]

    #dots = [d if abs(d) > 0.0001 else 1 for d in dots]

    #reg_dots = [d/(10*(max(dots)-min(dots))) for d in dots]

    #if eps > 0 and eps < 1:
    #  #normals = [(d*n[0]/abs(d), d*n[1]/abs(d)) for d, n in zip(dots, normals)]
    #  normals = [(s*n[0], s*n[1]) for s, n in zip(speeds, normals)]

    #normals = [(s[0]*n[0], s[1]*n[1]) for s, n in zip(speeds, normals)]

    mid = midpoints([geom.Point(pt) for pt in zip(*poly.exterior.coords.xy)[:-1]])
    x, y = zip(*[(pt.xy[0][0], pt.xy[1][0]) for pt in mid])
    u, v = zip(*normals)
    dx, dy = zip(*speeds)

    if eps > 0 and eps < 1:
      if all_points is None:
        all_points = [[np.array([[xi, yi]])] for xi, yi in zip(x, y)]
      else:
        for l, xi, yi in zip(all_points, x, y):
          l.append(np.array([[xi, yi]]))

    np.savetxt('polygon_normals_{}'.format(i), np.vstack([x, y, u, v, dx, dy]).T,
               header='x y u v dx dy', comments='')

  for i, l in enumerate(all_points):
    np.savetxt('trajectory_{}'.format(i), np.vstack(l), header='x y', comments='')


def check_speeds():
  perc = 0.3
  points = [(0.5, 0.5)]#, (0.5, 0), (0.5, 1.)]

  interp = PolygonInterpolator(p1s, p2s)
  poly = interp.fast_interpolate(perc)
  pairs = interp.fast_interpolate_pairs(perc)

  fig = plt.figure()
  ax = fig.add_subplot('111', aspect='equal')

  for p, c in zip([poly, p1s, p2s], ['r', 'b--', 'g--']):
    x, y = p.exterior.xy
    plt.plot(x, y, c)

  pair_x, pair_y = zip(*pairs)

  points_qhull = np.vstack(pairs)
  qhull = spatial.ConvexHull(points_qhull)

  spd_x, spd_y = zip(*interp.point_derivative(1.))
  ax.quiver(pair_x, pair_y, spd_x, spd_y, color='r')

  for i, xy in enumerate(zip(pair_x, pair_y)):
    ax.annotate(str(i), xy)

  normals, _ = interp.normals_offset(perc)
  nx, ny = zip(*normals)

  for point in points:
    proj = interp.projections(point, perc)
    speeds = interp.point_dist_derivative(point, perc, 1.)

    x, y = zip(*proj)
    spd_x, spd_y = zip(*speeds)

    px, py = point
    plt.plot(px, py, 'or')
    plt.plot(x, y, 'ob')
    ax.quiver(x, y, spd_x, spd_y, color="b")
    ax.quiver(x, y, nx, ny)

  plt.show()


def launch_animation():
  global interp, quiv, quiv_spd, quiv_dot
  ##-- for animation purposes
  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                       #xlim=(-0.2, 0.2), ylim=(-0.2, 0.2))
                       xlim=(-1.2, 2.2), ylim=(-1.2, 2.2))
  ax.grid()
  line1, = ax.plot([], [], 'r')
  line2, = ax.plot([], [], 'b')
  line3, = ax.plot([], [], 'g--')

  quiv = ax.quiver([], [], [], [], color='g')
  quiv_spd = ax.quiver([], [], [], [], color='r')
  quiv_dot = ax.quiver([], [], [], [], color='b')

  nr_steps = 500

  interp = PolygonInterpolator(p1s, p2s)

  def init():
    global quiv, quiv_spd
    x, y = p1s.exterior.coords.xy
    line1.set_data(x, y)
    x, y = p2s.exterior.coords.xy
    line2.set_data(x, y)
    line3.set_data([], [])
    normals, _ = normals_offset(p2s)
    mid = midpoints([geom.Point(pt) for pt in zip(*p2s.exterior.coords.xy)[:-1]])
    x, y = zip(*[(pt.xy[0], pt.xy[1]) for pt in mid])
    u, v = zip(*normals)
    quiv.set_offsets(np.hstack([x, y]))
    quiv.set_UVC(u, v)

    speeds = interp.midpoint_derivative(1.)
    u, v = zip(*speeds)
    quiv_spd.set_offsets(np.hstack([x, y]))
    quiv_spd.set_UVC(u, v)

    dots = [n[0]*s[0]+n[1]*s[1] for n, s in zip(normals, speeds)]
    spd = [(d*n[0], d*n[1]) for n, d in zip(normals, dots)]

    u, v = zip(*spd)
    quiv_dot.set_offsets(np.hstack([x, y]))
    quiv_dot.set_UVC(u, v)

    return line1, line2, line3, quiv, quiv_spd, quiv_dot

  def animate(i):
    global interp, quiv
    perc = float(i)/float(nr_steps)
    p = interp.interpolate(perc)
    normals, offset = normals_offset(p)
    #mid = midpoints(geom.Point(map(as_tuple, p.exterior.coords)))
    line3.set_data(*p.exterior.coords.xy)
    #print p.exterior.coords.xy
    #normals = interp.midpoint_derivative(1)
    #normals = interp.normal_derivative(0.5)
    #interp.normal_derivative(1/float(nr_steps))
    #x = [point.xy[0][0] for point in mid]
    #y = [point.xy[1][0] for point in mid]
    mid = midpoints([geom.Point(pt) for pt in zip(*p.exterior.coords.xy)[:-1]])
    x, y = zip(*[(pt.xy[0], pt.xy[1]) for pt in mid])
    u, v = zip(*normals)
    quiv.set_offsets(np.hstack([x, y]))
    quiv.set_UVC(u, v)

    speeds = interp.point_dist_derivative((1., 0.), perc, 1.)
    u, v = zip(*speeds)
    quiv_spd.set_offsets(np.hstack([x, y]))
    quiv_spd.set_UVC(u, v)

    dots = [n[0]*s[0]+n[1]*s[1] for n, s in zip(normals, speeds)]
    spd = [(d*n[0], d*n[1]) for n, d in zip(normals, dots)]

    u, v = zip(*spd)
    quiv_dot.set_offsets(np.hstack([x, y]))
    quiv_dot.set_UVC(u, v)
    return line1, line2, line3, quiv, quiv_spd

  ani = animation.FuncAnimation(fig, animate, frames=nr_steps,
                                interval=0.2, blit=True, init_func=init,
                                repeat=True)
  plt.show()

if __name__ == '__main__':
  check_speeds()
