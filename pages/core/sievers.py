"""
                    Waveplot: An online wavefunction viewer
                    Copyright (C) 2023  Jon G. C. Kragskow

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import numpy as np
import numpy.linalg as la
import py3nj
import functools
import inorgqm.multi_electron as iqme
from scipy.spatial import Delaunay

from . import utils as ut

def cf_r(theta, phi, n, l, s, j, mj):

    rho = rho(n, j, mj, l, s)


    r = 


def sievers_r(theta, a_2, a_4, a_6):

    c_0 = 3./(4.*np.pi)
    c_2 = a_2 / np.sqrt(4. * np.pi / 5)
    c_4 = a_4 / np.sqrt(4. * np.pi / 9)
    c_6 = a_6 / np.sqrt(4. * np.pi / 13)

    # Calculate r, x, y, z values for each theta
    r = c_0
    r += c_2 * 0.25 * np.sqrt(5/np.pi) * (3 * np.cos(theta)**2 - 1)
    r += c_4 * 3/16 * np.sqrt(1/np.pi) * (35 * np.cos(theta)**4 - 30 * np.cos(theta)**2 + 3) # noqa
    r += c_6 * 1/32 * np.sqrt(13/np.pi) * (231 * np.cos(theta)**6 - 315 * np.cos(theta)**4 + 105 * np.cos(theta)**2 - 5) # noqa
    r = r**(1./3)

    return r


def tri_normal(vertices: list['Vector']):

    n = np.cross(
        vertices[1].pos-vertices[0].pos,
        vertices[2].pos-vertices[0].pos
    )
    vertices[0].normal += n
    vertices[1].normal += n
    vertices[2].normal += n

    return n

def compute_trisurf(a_2, a_4, a_6, scale=2):

    # Create angular grid, with values of theta
    phi = np.linspace(0, np.pi, 51)
    theta = np.linspace(0, np.pi*2, 51)
    u, v = np.meshgrid(phi, theta)
    u = u.flatten()
    v = v.flatten()
    r = sievers_r(v, a_2, a_4, a_6)*scale

    #r = walter_r(v, u)

    # Points on 2d grid
    points2D = np.vstack([u, v]).T
    tri = Delaunay(points2D)
    verts_to_simp = tri.simplices

    # coordinates of sievers surface
    x = r * np.sin(v)*np.cos(u)
    y = r * np.sin(v)*np.sin(u)
    z = r * np.cos(v)
    vertices = np.array([x, y, z]).T

    # Rotate spheroid
    vertices = ut.set_z_alignment(vertices, [0,1,0], [0,0,1])

    vertices = np.array([
        Vector(vertex)
        for vertex in vertices
    ])

    # Calculate norm of each triangle
    for simp in verts_to_simp:
        tri_normal(vertices[simp])

    normals = np.array(
        [
            vertex.normal 
            if la.norm(vertex.normal) < 1E-9 else vertex.normal
            for vertex in vertices
        ]
    )

    vertices = np.array(
        [
            vertex.pos
            for vertex in vertices
        ]
    )

    return vertices, verts_to_simp, normals


class Vector():
    def __init__(self, pos) -> None:
        self.pos = pos
        self.normal = np.zeros(3)
        pass


@functools.lru_cache(maxsize=32)
def compute_isosurface(a_2, a_4, a_6, n_x, n_y, n_z, scale=1, comment=""):

    x_min, x_max = -0.4, 0.4
    y_min, y_max = -0.4, 0.4
    z_min, z_max = -0.4, 0.4

    x = np.linspace(x_min, x_max, n_x)
    y = np.linspace(y_min, y_max, n_y)
    z = np.linspace(z_min, z_max, n_z)

    x_step = x[1]-x[0]
    y_step = y[1]-y[0]
    z_step = z[1]-z[0]

    iso = np.zeros([n_x, n_y, n_z])
    for xit, xval in enumerate(x):
        for yit, yval in enumerate(y):
            for zit, zval in enumerate(z):
                r_grid = np.sqrt(xval**2 + yval**2 + zval**2)
                if r_grid > 0:
                    theta = np.arccos(yval/r_grid)
                else:
                    theta = 0.
                phi = np.arctan2(xval, zval)
                r_func = sievers_r(theta, phi, a_2, a_4, a_6)
                if r_grid-r_func < 0:
                    iso[xit, yit, zit] = np.abs(r_grid-r_func)

    x_min *= scale
    y_min *= scale
    z_min *= scale
    x_step *= scale
    y_step *= scale
    z_step *= scale

    c = ""

    c += "Use isoval of 0.55 for visualisation\n"
    c += "{}\n".format(comment)
    c += "1     {:.6f} {:.6f} {:.6f}\n".format(x_min, y_min, z_min)
    c += "{:d}   {:.6f}    0.000000    0.000000\n".format(n_x, x_step)
    c += "{:d}   0.000000    {:.6f}    0.000000\n".format(n_y, y_step)
    c += "{:d}   0.000000    0.000000    {:.6f}\n".format(n_z, z_step)
    c += " 8   0.000000    0.000000   0.000000  -0.165000\n"

    a = 0

    for xit, xval in enumerate(x):
        for yit, yval in enumerate(y):
            for zit, zval in enumerate(z):
                a += 1
                c += "{:.5e} ".format(iso[xit, yit, zit])
                if a == 6:
                    c += "\n"
                    a = 0
            c += "\n"
            a = 0

    return c


@functools.lru_cache(maxsize=32)
def wigner3(a, b, c, d, e, f):

    a = int(2*a)
    b = int(2*b)
    c = int(2*c)
    d = int(2*d)
    e = int(2*e)
    f = int(2*f)

    return py3nj.wigner3j(a, b, c, d, e, f)


@functools.lru_cache(maxsize=32)
def wigner6(a, b, c, d, e, f):

    a = int(2*a)
    b = int(2*b)
    c = int(2*c)
    d = int(2*d)
    e = int(2*e)
    f = int(2*f)

    return py3nj.wigner6j(a, b, c, d, e, f)


def compute_a_vals(n, J, mJ, L, S):

    k_max = min(6, int(2*J+1))

    a_vals = []

    if n == 7:
        a_vals = [0., 0., 0.]
    elif n < 7:
        for k in range(2, k_max+2, 2):
            a_vals.append(_compute_light_a_val(n, J, mJ, L, S, k))
    else:
        for k in range(2, k_max+2, 2):
            a_vals.append(_compute_heavy_a_val(J, mJ, n, k))

    return a_vals


def _compute_light_a_val(n, J, mJ, L, S, k):

    a_k = np.sqrt(4*np.pi/(2*k+1))
    a_k *= (-1)**(2*J-mJ+L+S)
    a_k *= 7./(np.sqrt(4*np.pi)) * (2*J + 1) * np.sqrt(2*k+1)
    a_k *= wigner3(J, k, J, -mJ, 0, mJ)/wigner3(L, k, L, -L, 0, L)
    a_k *= wigner6(L, J, S, J, L, k)
    a_k *= wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n+1):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4-it), (it-4))
    a_k *= summa

    return a_k


def _compute_light_rho_k(n, J, mJ, L, S, k):

    rho_k = (-1)**(2*J-mJ+L+S)
    rho_k *= 7./(np.sqrt(4*np.pi)) * (2*J + 1) * np.sqrt(2*k+1)
    rho_k *= wigner6(L, J, S, J, L, k)
    rho_k *= wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n+1):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4-it), (it-4))
    rho_k *= summa

    return rho_k


def _compute_heavy_rho_k(n, k):

    rho_k = 7/(np.sqrt(4*np.pi))
    rho_k *= np.sqrt(2*k+1)*wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n-6):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4-it), (it-4))
    rho_k *= summa

    return    


def _compute_heavy_a_val(J, mJ, n, k):
    a_k = _compute_heavy_rho_k(n, k)
    a_k *= np.sqrt(4*np.pi/(2*k+1))
    a_k *= (-1)**(J-mJ)
    a_k *= wigner3(J, k, J, -mJ, 0, mJ)/wigner3(J, k, J, -J, 0, J)
    return a_k


def compute_CF_coeffs(J):

    cfps = np.loadtxt("cfps.dat")
    k_max = 6

    _, _, jz, jp, jm, _ = iqme.calc_ang_mom_ops(J)

    stev_ops = iqme.calc_stev_ops(k_max, jp, jm, jz)[1::2]

    hcf, vals, vecs = iqme.calc_HCF(J, cfps, stev_ops, kmax=k_max)

    k_vals = [2, 4, 6]
    for kit in range(k_max/2):
        k = k_vals[kit]
        for qit in range(2*k+1):
            stev_ops[kit, qit] = la.inv(vecs) @ stev_ops[kit, qit] @ vecs
