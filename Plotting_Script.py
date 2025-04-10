import numpy as np
import plotly.graph_objects as go
#import plotly.express as px
#from plotly.subplots import make_subplots

import FuncDef as fd
from importlib import reload
import itertools

savefolder = "/savedplots/"

FACTOR = 1653 #from unitless to n/cm-3

#load data from cloud-data.csv
data = fd.load_data("cloud-data.csv", sep=";")

#Example at Cloud 18
cloud = fd.Cloud(*data.loc["18"])

points = cloud.getPlottableDust(40, 40, 40)

#Query E+ dustmap in surrounding area
x, y, z = fd.cylinder(cloud.r, cloud.h, cloud.a, cloud.b, cloud.t, cloud.nr, cloud.nh, cloud.nt)

x_ext = np.linspace(np.min(x)-100, np.max(x)+100, 40)
y_ext = np.linspace(np.min(y)-100, np.max(y)+100, 40)
z_ext = np.linspace(np.min(z)-100, np.max(z)+100, 40)
x_ext, y_ext, z_ext = np.meshgrid(x_ext, y_ext, z_ext)
val = fd.query_region(x_ext.flatten(), y_ext.flatten(), z_ext.flatten())
val[np.where(fd.check_inside_arr(x_ext.flatten(), y_ext.flatten(), z_ext.flatten(), cloud.r, cloud.h, cloud.a, cloud.b, cloud.t))]=0

#create volume of cloud with inner structure visible
volume = go.Volume(
    x = cloud.points[:,0],
    y = cloud.points[:,1],
    z = cloud.points[:,2],
    value = FACTOR*cloud.points[:,3],
    isomin = FACTOR*1e-4,#cloud.v_f_min,
    isomax = FACTOR*3e-3,#cloud.v_f_max,
    opacity = 0.2,
    surface_count = 25,
    colorscale = "plasma",#[[0, "olive"], [1, "olive"]],
    opacityscale = [[0, 1], [1, 1]],
    showscale = True,
    name = "Cloud " + "18",
    showlegend = False,
    colorbar=dict(
        title="Dust density [n/cm³]",
        titleside="top",
    )
)

#create volume element of surrounding dust
volume_surr = go.Volume(
    x = x_ext.flatten(),
    y = y_ext.flatten(),
    z = z_ext.flatten(),
    value = FACTOR*val,
    isomin = FACTOR*1e-4,#cloud.v_f_min,
    isomax = FACTOR*3e-3,#cloud.v_f_max,
    opacity = 0.2,
    surface_count = 25,
    colorscale = "Greys",#[[0, "olive"], [1, "olive"]],
    opacityscale = [[0, 1], [1, 1]],
    showscale = True,
    name = "Density",
    showlegend = False,
    colorbar=dict(
        title="Dust density [n/cm³]",
        titleside="top",
    )
)

#a few constants useful for the layout, including calculating a camera position at a desired zoom level
xrange = [np.min(x_ext), np.max(x_ext)]
yrange = [np.min(y_ext), np.max(y_ext)]
zrange = [np.min(z_ext), np.max(z_ext)]

xscale = (xrange[1] - xrange[0])/(zrange[1] - zrange[0])
yscale = (yrange[1] - yrange[0])/(zrange[1] - zrange[0])

norm = np.sqrt(3)*1.25
cam_weights = [1, 1/8, 0]
cam_pos = np.sqrt(3/(np.sum(np.array(cam_weights)**2)))*np.array(cam_weights)

#define layout
layout = go.Layout(
    template="plotly_white",
    paper_bgcolor="white",
    plot_bgcolor="white",
    scene=dict(
        aspectmode="manual",
        aspectratio=dict(x=xscale, y=yscale, z=1),
        xaxis=dict(range=xrange),
        yaxis=dict(range=yrange),
        zaxis=dict(range=zrange),
        camera=dict(
            eye=dict(x=cam_pos[0], y=cam_pos[1], z=cam_pos[2]),
        ),
    ),
    width=1000,
    height=1000,
)

#create first figure
fig = go.Figure(data=[volume_surr, volume],layout=layout)
f1 = go.FigureWidget(fig)
fig.write_html(savefolder + "Cloud_selected.html")
print("done Cloud_selected")

#visualize cloud boundary 
volume2 = go.Volume(
    x = cloud.points[:,0],
    y = cloud.points[:,1],
    z = cloud.points[:,2],
    value = FACTOR*cloud.points[:,3],
    isomin = FACTOR*cloud.v_f_min,
    isomax = FACTOR*cloud.v_f_max,
    opacity = 0.75,
    surface_count = 1,
    colorscale = [[0, "olive"], [1, "olive"]],
    opacityscale = [[0, 0], [1, 1]],
    showscale = False,
    name = "Cloud " + "18",
    showlegend = False,
    hovertext = "Cutoff: " + str(FACTOR*(cloud.v_f_min+cloud.v_f_max)/2) + " n/cm^3"
)
fig = go.Figure(data=[volume_surr, volume2],layout=layout)
f2 = go.FigureWidget(fig)
fig.write_html(savefolder + "Dust_cutoff.html")
print("done Dust_cutoff")

# zoom in, query region again
x_ext = np.linspace(np.min(x)-15, np.max(x)+15, 40)
y_ext = np.linspace(np.min(y)-15, np.max(y)+15, 40)
z_ext = np.linspace(np.min(z)-15, np.max(z)+15, 40)
x_ext, y_ext, z_ext = np.meshgrid(x_ext, y_ext, z_ext)
val = fd.query_region(x_ext.flatten(), y_ext.flatten(), z_ext.flatten())
val[np.where(fd.check_inside_arr(x_ext.flatten(), y_ext.flatten(), z_ext.flatten(), cloud.r, cloud.h, cloud.a, cloud.b, cloud.t))]=0

#ranges and scale has changed
xrange = [np.min(x_ext), np.max(x_ext)]
yrange = [np.min(y_ext), np.max(y_ext)]
zrange = [np.min(z_ext), np.max(z_ext)]

xscale = (xrange[1] - xrange[0])/(zrange[1] - zrange[0])
yscale = (yrange[1] - yrange[0])/(zrange[1] - zrange[0])


layout = go.Layout(
    template="plotly_white",
    paper_bgcolor="white",
    plot_bgcolor="white",
    scene=dict(
        aspectmode="manual",
        aspectratio=dict(x=xscale, y=yscale, z=1),
        xaxis=dict(range=xrange),
        yaxis=dict(range=yrange),
        zaxis=dict(range=zrange),
        camera=dict(
            eye=dict(x=cam_pos[0], y=cam_pos[1], z=cam_pos[2]),
        ),
    ),
    width=1000,
    height=1000,
)

#visualize surface of the cylinder
xc, yc, zc = cloud.getCylinderPts()
cyl = go.Surface(
    x = xc[:,-1,:],
    y = yc[:,-1,:],
    z = zc[:,-1,:],
    colorscale = [[0, "green"], [1, "green"]],
    opacity = 0.2,
    showscale = False,
    name="cyl"
)

#new surrounding dust
volume_surr = go.Volume(
    x = x_ext.flatten(),
    y = y_ext.flatten(),
    z = z_ext.flatten(),
    value = FACTOR*val,
    isomin = FACTOR*1e-4,#cloud.v_f_min,
    isomax = FACTOR*3e-3,#cloud.v_f_max,
    opacity = 0.13,
    surface_count = 25,
    colorscale = "Greys",#[[0, "olive"], [1, "olive"]],
    opacityscale = [[0, 1], [1, 1]],
    showscale = True,
    name = "Cloud " + "18",
    showlegend = False,
    colorbar=dict(
        title="Dust density [n/cm³]",
        titleside="top",
    )
)

fig = go.Figure(data=[volume_surr, volume2, cyl],layout=layout)
f3 = go.FigureWidget(fig)
fig.write_html(savefolder + "Cloud_selected_cylinder.html")
print("done Cloud_selected_cylinder")


#parameter variation (See FuncDef)
params = np.array(list(itertools.product(
            (0.9, 1, 1.1),          # ir (radius variation)
            (0.9, 1, 1.1),          # iv_f_max (cutoff variation)
            (-20, 0, 20),           # it1 (x-translation)
            (-20, 0, 20),           # it2 (y-translation)
            (-20, 0, 20),           # it3 (z-translation)
            (-0.25, 0, 0.25),       # ia (angle a variation)
            (-0.25, 0, 0.25)        # ib (angle b variation)
        )))

cyl_params = params[np.random.randint(0, len(params), 4)]

#4 random cylinder surfaces
cyls = []
for i in range(4):
    var = cyl_params[i]
    xc, yc, zc = fd.cylinder(
        cloud.r*var[0],
        cloud.h,
        cloud.a+(var[5]/18)*np.pi,
        cloud.b+(var[6]/18)*np.pi,
        [cloud.t[0]+var[2], cloud.t[1]+var[3], cloud.t[2]+var[4]])
    cyls.append(go.Surface(
        x = xc[:,-1,:],
        y = yc[:,-1,:],
        z = zc[:,-1,:],
        colorscale = [[0, "red"], [1, "red"]],
        opacity = 0.08,
        showscale = False,
        name=str("Cylinder " + str(i))
    ))

cyl["opacity"] = 0.3
volume_surr["opacity"] = 0.1

fig = go.Figure(data=[volume_surr, volume2, cyl, *cyls],layout=layout)
f4 = go.FigureWidget(fig)
fig.write_html(savefolder + "Cloud_error_cylinders.html")
print("done all")
