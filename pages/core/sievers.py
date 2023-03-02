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
import scipy.special as ssp
from . import utils as ut


def sievers_r(theta, a_2, a_4, a_6):

    c_0 = 3. / (4. * np.pi)
    c_2 = a_2 / np.sqrt(4. * np.pi / 5)
    c_4 = a_4 / np.sqrt(4. * np.pi / 9)
    c_6 = a_6 / np.sqrt(4. * np.pi / 13)

    # Calculate r, x, y, z values for each theta
    r = c_0
    r += c_2 * 0.25 * np.sqrt(5 / np.pi) * (3 * np.cos(theta)**2 - 1)
    r += c_4 * 3/16 * np.sqrt(1 / np.pi) * (35 * np.cos(theta)**4 - 30 * np.cos(theta)**2 + 3) # noqa
    r += c_6 * 1/32 * np.sqrt(13 / np.pi) * (231 * np.cos(theta)**6 - 315 * np.cos(theta)**4 + 105 * np.cos(theta)**2 - 5) # noqa
    r = r**(1. / 3)

    return r


def tri_normal(vertices: list['Vector']):

    n = np.cross(
        vertices[1].pos - vertices[0].pos,
        vertices[2].pos - vertices[0].pos
    )
    vertices[0].normal += n
    vertices[1].normal += n
    vertices[2].normal += n

    return n


def compute_trisurf(a_2, a_4, a_6, scale=2):

    # Create angular grid, with values of theta
    phi = np.linspace(0.0001, np.pi, 51)
    theta = np.linspace(0.0001, 2 * np.pi, 51)
    u, v = np.meshgrid(phi, theta)
    u = u.flatten()
    v = v.flatten()
    r = sievers_r(v, a_2, a_4, a_6) * scale

    # Points on 2d grid
    points2D = np.vstack([u, v]).T
    tri = Delaunay(points2D)
    verts_to_simp = tri.simplices

    # coordinates of sievers surface
    x = r * np.sin(v) * np.cos(u)
    y = r * np.sin(v) * np.sin(u)
    z = r * np.cos(v)
    vertices = np.array([x, y, z]).T

    # Rotate spheroid
    vertices = ut.set_z_alignment(vertices, [0, 1, 0], [0, 0, 1])

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
def wigner3(a, b, c, d, e, f):

    a = int(2 * a)
    b = int(2 * b)
    c = int(2 * c)
    d = int(2 * d)
    e = int(2 * e)
    f = int(2 * f)

    return py3nj.wigner3j(a, b, c, d, e, f)


@functools.lru_cache(maxsize=32)
def wigner6(a, b, c, d, e, f):

    a = int(2 * a)
    b = int(2 * b)
    c = int(2 * c)
    d = int(2 * d)
    e = int(2 * e)
    f = int(2 * f)

    return py3nj.wigner6j(a, b, c, d, e, f)


def compute_a_vals(n, J, mJ, L, S):

    k_max = min(6, int(2 * J + 1))

    a_vals = []

    if n == 7:
        a_vals = [0., 0., 0.]
    elif n < 7:
        for k in range(2, k_max + 2, 2):
            a_vals.append(_compute_light_a_val(n, J, mJ, L, S, k))
    else:
        for k in range(2, k_max + 2, 2):
            a_vals.append(_compute_heavy_a_val(J, mJ, n, k))

    return a_vals


def _compute_light_a_val(n, J, mJ, L, S, k):

    a_k = np.sqrt(4 * np.pi / (2 * k + 1))
    a_k *= (-1)**(2 * J - mJ + L + S)
    a_k *= 7. / (np.sqrt(4 * np.pi)) * (2 * J + 1) * np.sqrt(2 * k + 1)
    a_k *= wigner3(J, k, J, -mJ, 0, mJ) / wigner3(L, k, L, -L, 0, L)
    a_k *= wigner6(L, J, S, J, L, k)
    a_k *= wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n + 1):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4 - it), (it - 4))
    a_k *= summa

    return a_k


def _compute_heavy_a_val(J, mJ, n, k):
    a_k = _compute_heavy_rho_k(n, J, k)
    a_k *= np.sqrt(4 * np.pi / (2 * k + 1))
    a_k *= (-1)**(J - mJ)
    a_k *= wigner3(J, k, J, -mJ, 0, mJ)
    return a_k


def rho(n, J, mJ, L, S, k):

    if n < 7:
        rho = _compute_light_rho_k(n, J, mJ, L, S, k)
    else:
        rho = _compute_heavy_rho_k(n, J, k)

    return rho


def _compute_light_rho_k(n, J, mJ, L, S, k):

    rho_k = (-1)**(2 * J - mJ + L + S)
    rho_k *= 7. / (np.sqrt(4 * np.pi))
    rho_k *= 2 * J + 1
    rho_k *= np.sqrt(2 * k + 1)
    rho_k *= wigner6(L, J, S, J, L, k)
    rho_k *= wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n + 1):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4 - it), (it - 4))
    rho_k *= summa

    return rho_k


def _compute_heavy_rho_k(n, j, k):

    if k == 0:
        pre = 1
    else:
        pre = 0.

    rho_k = np.sqrt(2 * k + 1) * wigner3(k, 3, 3, 0, 0, 0)
    summa = 0
    for it in range(1, n - 6):
        summa += (-1)**it * wigner3(k, 3, 3, 0, (4 - it), (it - 4))
    rho_k *= summa

    rho_k += pre

    rho_k *= 7/(np.sqrt(4*np.pi))

    rho_k /= wigner3(j, k, j, -j, 0, j)

    return rho_k


def cf_density(n, l, s, j, mj, soi):

    # Create angular grid, with values of theta
    phi = np.linspace(-np.pi/2, np.pi/2, 51)
    theta = np.linspace(0., 2 * np.pi, 51)
    u, v = np.meshgrid(phi, theta)
    u = u.flatten()
    v = v.flatten()

    rho_k = [
        rho(n, j, mj, l, s, k)
        for k in [0, 2, 4, 6]
    ]

    # rho_k[0] = n/np.sqrt(np.pi*4)

    rho_k[0] = 3/(4*np.pi)

    cfps = np.loadtxt("cfps.dat")
    k_max = 6

    _, _, jz, jp, jm, _ = iqme.calc_ang_mom_ops(j)

    stev_ops = iqme.calc_stev_ops(k_max, j, jp, jm, jz)[1::2]

    _, _, vecs = iqme.calc_HCF(j, cfps, stev_ops, k_max=k_max)

    # State of interest

    mj_values = np.arange(j, -j - 1, -1)

    # Calculate a_kq values

    a = np.zeros([4, 13])

    for kit, k in enumerate([0, 2, 4, 6]):
        for qit, q in enumerate(range(-k, k + 1)):
            if q == 0:
                for m_ind, m in enumerate(mj_values):
                    _tmp = np.abs(vecs[m_ind, soi])**2
                    _tmp *= (-1.)**(j - m)
                    _tmp *= wigner3(j, k, j, -m, 0, m)
                    a[kit, qit] += _tmp
            elif q > 0:
                for m1_ind, m1 in enumerate(mj_values):
                    for m2_ind, m2 in enumerate(mj_values):
                        if m2 < m1:
                            _tmp = np.real(np.conj(vecs[m1_ind, soi]) * vecs[m2_ind, soi])
                            _tmp *= (-1.)**(j - m + q)
                            _tmp *= wigner3(j, k, j, -m1, q, m2)
                            a[kit, qit] += np.sqrt(2) * _tmp
            elif q < 0:
                for m1_ind, m1 in enumerate(mj_values):
                    for m2_ind, m2 in enumerate(mj_values):
                        if m2 < m1:
                            _tmp = np.imag(np.conj(vecs[m1_ind, soi]) * vecs[m2_ind, soi])
                            _tmp *= (-1.)**(j - m + q)
                            _tmp *= wigner3(j, k, j, -m1, q, m2)
                            a[kit, qit] += np.sqrt(2) * _tmp

    # Density

    density = 0.

    for kit, k in enumerate([0, 2, 4, 6]):
        for qit, q in enumerate(range(-k, k + 1)):
            if k == 0 and q == 0:
                density += rho_k[kit]
            elif q == 0:
                density += a[kit, qit] * rho_k[kit] * tess_h(k, q, v, u)
            elif q > 0:
                density += (a[kit, qit] * rho_k[kit] * tess_h(k, q, v, u))
                density += (a[kit, qit - 2 * np.abs(q)] * rho_k[kit] * tess_h(k, -q, v, u))

    # Points on 2d grid
    points2D = np.vstack([u, v]).T
    tri = Delaunay(points2D)
    verts_to_simp = tri.simplices
    # Radius as cube root of density
    r = np.real(density) ** (1. /3.)

    r *= 2

    # coordinates of sievers surface
    x = r * np.sin(v) * np.cos(u)
    y = r * np.sin(v) * np.sin(u)
    z = r * np.cos(v)
    vertices = np.array([x, y, z]).T

    # Rotate spheroid
    vertices = ut.set_z_alignment(vertices, [0, 1, 0], [0, 0, 1])

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


def tess_h(k, q, phi, theta):

    if q == 0:
        val = ssp.sph_harm(q, k, theta, phi)
    elif q > 0:
        val = ssp.sph_harm(-q, k, theta, phi)
        val += (-1)**q * ssp.sph_harm(q, k, theta, phi)
        val *= 1. / np.sqrt(2)
    elif q < 0:
        val = ssp.sph_harm(q, k, theta, phi)
        val -= (-1)**q * ssp.sph_harm(-q, k, theta, phi)
        val *= 1j / np.sqrt(2)

    return val
