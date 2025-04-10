import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from dustmaps.edenhofer2023 import Edenhofer2023Query

import plotly.graph_objects as go
from concurrent.futures import ProcessPoolExecutor
import json
from tqdm.contrib import itertools
import time

#def cylinder(r, h, a =0, nt=100, nv =50):
#    """
#    parametrize the cylinder of radius r, height h, base point a
#    """
#    theta = np.linspace(0, 2*np.pi, nt)
#    v = np.linspace(a, a+h, nv )
#    theta, v = np.meshgrid(theta, v)
#    x = r*np.cos(theta)
#    y = r*np.sin(theta)
#    z = v
#    return x, y, z

def load_data(filename, sep=";"):
    """
    Load the cloud data from the file, with formatting to the correct types
    """
    data = pd.read_csv(filename, sep=sep)
    data.set_index("i", inplace=True)

    data["a"] = data["a"].apply(eval)
    data["b"] = data["b"].apply(eval)
    data["t"] = data["t"].apply(json.loads)

    data["v_f_min"] = data["v_f_min"].apply(float)
    data["v_f_max"] = data["v_f_max"].apply(float)

    for col in data.columns[:-1]:
        data[col] = pd.to_numeric(data[col])
    return data

def filled_cylinder(r, h, nr=100, nh=100, nt =90):
    """
    parametrize the cylinder of radius r, height h, base point a
    """
    rr = np.linspace(0, r, nr)
    theta = np.linspace(0, 2*np.pi, nt)
    v = np.linspace(-h/2, h/2, nh)
    rr, theta, v = np.meshgrid(rr, theta, v)
    y = rr*np.cos(theta)
    z = rr*np.sin(theta)
    x = v
    return x, y, z

def RotationMatrix(a, b, c=0):
    """
    Return the rotation matrix for the angles a, b, c
    """
    #invert b from the wikipedia page
    Ra = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rb = np.array([[np.cos(b), 0, -np.sin(b)], [0, 1, 0], [np.sin(b), 0, np.cos(b)]])
    Rc = np.array([[1, 0, 0], [0, np.cos(c), -np.sin(c)], [0, np.sin(c), np.cos(c)]])

    return Ra @ Rb @ Rc

def cylinder(r, h, a, b, t, nr=100, nh=100, nt=90):
    """
    parametrize the cylinder of radius r, height h, base point a
    """
    x, y, z = filled_cylinder(r, h, nr, nh, nt)
    S = np.stack((x, y, z), axis=2)
    rotmat = RotationMatrix(a, b)
    S = rotmat@S

    xx, yy, zz = S[:, :, 0]+t[0], S[:, :, 1]+t[1], S[:, :, 2]+t[2]
    return xx, yy, zz

def check_inside(x, y, z, r, h, a, b, t):
    """
    check if the point (x, y, z) is inside the cylinder
    """
    x, y, z = np.array([x, y, z]) - t
    #invert the rotation
    rotmat = RotationMatrix(a, b)
    x, y, z = rotmat.T @ np.array([x, y, z])
    return z**2 + y**2 <= r**2 and abs(x) <= h/2

def dv_i(i, r, h, nr, nh, nt):
    """
    Calculate the volume element at the i-th slice
    """
    return h/nh*(r/nr)**2*(2*np.pi)/nt*(i+1/2)

def deltadv_i(i, r, h, nr, nh, nt):
    """
    Calculate the error in the volume element at the i-th slice
    """
    dr = 2*(r*h/(nr**2*nh*nt))
    dh = (r**2/(nr**2*nh*nt))
    dnr = -2*(r**2*h/(nr**3*nh*nt))
    dnh = -2*(r**2*h/(nr**2*nh**2*nt))
    dnt = -2*(r**2*h/(nr**2*nh*nt**2))
    # Error assumed to be sqrt(x), so (dr*r)**2 = dr**2*r etc.
    return 2*np.pi*(i+1/2)*np.sqrt(dr**2*r + dh**2*r + dnr**2*nr + dnh**2*nh + dnt**2*nt)

def query_region(x, y, z, std=False):
    """
    Returns Dust density at the given points, x,y,z as meshgrid
    """
    sc = SkyCoord(
        x.flatten() * u.pc,
        y.flatten() * u.pc,
        z.flatten() * u.pc,
        frame='galactic',
        representation_type = 'cartesian'
    )

    sc.representation_type = 'spherical'
    if std:
        dust = Edenhofer2023Query().query(sc, mode='std')
    else:
        dust = Edenhofer2023Query().query(sc)

    return dust

def query_dust_in_cyl(r, h, a, b, t, nx=20, ny=20, nz=20):
    """
    Returns Dust density, 0 if outside the cylinder
    """
    x, y, z = cylinder(r, h, a, b, t)

    xrange = (np.min(x), np.max(x))
    yrange = (np.min(y), np.max(y))
    zrange = (np.min(z), np.max(z))

    xt = np.linspace(xrange[0], xrange[1], nx)
    yt = np.linspace(yrange[0], yrange[1], ny)
    zt = np.linspace(zrange[0], zrange[1], nz)

    xt_grid, yt_grid, zt_grid = np.meshgrid(xt, yt, zt)

    dust = query_region(xt_grid, yt_grid, zt_grid)

    points_in_space = np.stack((xt_grid.flatten(), yt_grid.flatten(), zt_grid.flatten(), dust), axis=1)
    points_in_cyl = np.stack((xt_grid.flatten(), yt_grid.flatten(), zt_grid.flatten(), np.zeros_like(dust)), axis=1)

    for i in range(len(points_in_space)):
        if check_inside(points_in_space[i][0], points_in_space[i][1], points_in_space[i][2], r, h, a, b, t):
            points_in_cyl[i][3] = points_in_space[i][3]
    
    return points_in_cyl

def massCalc_simple(r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t):
    """
    Calculate the mass of the dust in the cylinder
    """
    x, y, z = cylinder(r, h, a, b, t, nr, nh, nt)
    dust = query_region(x, y, z)
    dust = dust.reshape(x.shape)
    mass = 0
    for i in np.arange(nr):
        mass += np.sum(dust[:,i,:]*dv_i(i, r, h, nr, nh, nt), where=dust[:,i,:]>=(v_f_min+v_f_max)/2)
    return mass*55.86893061840122 #mass in solar mass, 1653cm^-3(to pc^-3)*m_p*1.37/1.989e33

def process_param(param, r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t):
    """
    Worker function for parallel parameter variation.
    """
    ir, iv_f_max, it1, it2, it3, ia, ib = param
    mass = massCalc_simple(
        r * ir,  # Scale radius
        h,        # Height (unchanged)
        a + (ia / 18) * np.pi,  # Angle a variation
        b + (ib / 18) * np.pi,  # Angle b variation
        nr, nh, nt,
        v_f_min,
        v_f_max * iv_f_max,     # Scale velocity
        [t[0] + it1, t[1] + it2, t[2] + it3]  # Translation
    )
    return (param, mass)

def massCalc(r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t, var=False, folder="Bachelor/Project/", cloud_id=None):
    dictofmass = {}
    
    if var:
        # Generate parameter combinations
        params = list(itertools.product(
            (0.9, 1, 1.1),          # ir (radius variation)
            (0.9, 1, 1.1),          # iv_f_max (velocity variation)
            (-20, 0, 20),           # it1 (x-translation)
            (-20, 0, 20),           # it2 (y-translation)
            (-20, 0, 20),           # it3 (z-translation)
            (-0.25, 0, 0.25),       # ia (angle a variation)
            (-0.25, 0, 0.25)        # ib (angle b variation)
        ))
        total_params = len(params)
        start_time_cloud = time.time()
        # Parallel execution
        with ProcessPoolExecutor(max_workers=6) as executor:
            # Submit jobs with ALL required parameters
            futures = [
                executor.submit(
                    process_param,
                    p,  
                    r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t  
                ) for p in params
            ]
            
            # Collect results
            for idx, future in enumerate(futures, 1):
                param, mass = future.result()
                dictofmass[param] = mass

                # Optional: Save progress periodically
                if idx % 243 == 0 or idx == total_params:
                    elapsed = time.time() - start_time_cloud
                    progress = idx / total_params
                    remaining = (elapsed / idx) * (total_params - idx) if idx > 0 else 0
                    
                    # Per-cloud ETA
                    with open(folder + "progress.txt", "a") as f:
                        f.write(
                            f"Cloud {cloud_id}: Progress: {progress*100:.1f}% | "
                            f"Elapsed: {elapsed/60:.1f}m | ETA: {remaining/60:.1f}m\n"
                        )
    else:
        # Non-variation case (unchanged)
        mass = massCalc_simple(r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t)
        dictofmass[0] = mass
    
    return dictofmass

class Cloud():
    def __init__(self, r, h, a, b, nr, nh, nt, v_f_min, v_f_max, t):
        self.r = r
        self.h = h
        self.a = a
        self.b = b
        self.t = t
        self.nr = nr
        self.nh = nh
        self.nt = nt
        self.v_f_min = v_f_min
        self.v_f_max = v_f_max
        self.mass = None
        self.points = None
    
    def getCartesian(self, nx=20, ny=20, nz=20):
        """
        Return cartesian coordinates for all points within bounding box of the cylinder coordinates
        """
        x, y, z = cylinder(self.r, self.h, self.a, self.b, self.t, self.nr, self.nh, self.nt)
        xt = np.linspace(np.min(x), np.max(x), nx)
        yt = np.linspace(np.min(y), np.max(y), ny)
        zt = np.linspace(np.min(z), np.max(z), nz)
        xt_grid, yt_grid, zt_grid = np.meshgrid(xt, yt, zt)
        return xt_grid, yt_grid, zt_grid
    
    def getCylinderPts(self):
        """
        Return the points in the cloud cylinder
        """
        return cylinder(self.r, self.h, self.a, self.b, self.t, self.nr, self.nh, self.nt)
    
    def getPlottableDust(self, nx=20, ny=20, nz=20):
        """
        Returns all points in the whole bounding box around the cylinder, but all dust points are set to 0 if outside the cylinder, save as self.points
        """
        if self.points == None:
            self.points = query_dust_in_cyl(self.r, self.h, self.a, self.b, self.t, nx, ny, nz)

    def getMass(self):
        """
        Calculate the mass of the cloud, save as self.mass
        """
        if self.mass is None:
            self.mass = massCalc_simple(self.r, self.h, self.a, self.b, self.nr, self.nh, self.nt, self.v_f_min, self.v_f_max, self.t)
    
    def getMassStatError(self):
        """
        Calculate the statistical error in the mass of the cloud
        """
        x, y, z = cylinder(self.r, self.h, self.a, self.b, self.t, self.nr, self.nh, self.nt)
        dust_std = query_region(x, y, z, std=True)
        dust_mat = query_region(x, y, z)
        dust_std = dust_std.reshape(x.shape)
        dust_mat = dust_mat.reshape(x.shape)
        sum = 0
        for i in np.arange(self.nr):
            sum += np.sum((dust_std[:,i,:]*dv_i(i, self.r, self.h, self.nr, self.nh, self.nt))**2, where=dust_mat[:,i,:]>=(self.v_f_min+self.v_f_max)/2)
            sum += np.sum((dust_mat[:,i,:]*deltadv_i(i, self.r, self.h, self.nr, self.nh, self.nt))**2, where=dust_mat[:,i,:]>=(self.v_f_min+self.v_f_max)/2)
        return np.sqrt(sum)*55.86893061840122

    def getVolume(self, color="red", num: str="100", nx=20, ny=20, nz=20):
        """
        Return the volume object for the cloud
        """
        if self.points is None:
            self.getPlottableDust(nx, ny, nz)
        volume = go.Volume(
            x = self.points[:,0],
            y = self.points[:,1],
            z = self.points[:,2],
            value = self.points[:,3],
            isomin = self.v_f_min,
            isomax = self.v_f_max,
            opacity = 0.3,
            surface_count = 1,
            colorscale = [[0, color], [1, color]],
            opacityscale = [[0,0], [1, 1]],
            showscale = False,
            name = "Cloud " + num,
            showlegend = True,
            spaceframe = dict(show=False)
        )

        return volume
