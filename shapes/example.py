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

def launch_animation():
  global interp, quiv
  ##-- for animation purposes
  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                       xlim=(-0.2, 0.2), ylim=(-0.2, 0.2))
  ax.grid()
  line1, = ax.plot([], [], 'r')
  line2, = ax.plot([], [], 'b')
  line3, = ax.plot([], [], 'g--')

  quiv = ax.quiver([], [], [], [], color='g')

  nr_steps = 500

  interp = PolygonInterpolator(p1s, p2s)

  def init():
    global quiv
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
    return line1, line2, line3, quiv

  def animate(i):
    global interp, quiv
    perc = float(i)/float(nr_steps)
    p = interp.interpolate(perc)
    normals, offset = normals_offset(p)
    #mid = midpoints(geom.Point(map(as_tuple, p.exterior.coords)))
    line3.set_data(*p.exterior.coords.xy)
    print p.exterior.coords.xy
    #x = [point.xy[0][0] for point in mid]
    #y = [point.xy[1][0] for point in mid]
    mid = midpoints([geom.Point(pt) for pt in zip(*p.exterior.coords.xy)[:-1]])
    x, y = zip(*[(pt.xy[0], pt.xy[1]) for pt in mid])
    u, v = zip(*normals)
    quiv.set_offsets(np.hstack([x, y]))
    quiv.set_UVC(u, v)
    return line1, line2, line3, quiv

  ani = animation.FuncAnimation(fig, animate, frames=nr_steps,
                                interval=0.2, blit=True, init_func=init,
                                repeat=True)
  plt.show()
