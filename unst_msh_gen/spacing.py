
""" Mesh spacing utilities, for barotropic flows
"""

# Authors: Darren Engwirda, Ali Salimi-Tarazouj

# routines to compute mesh spacing functions that scale
# with shallow-water wave lengths, elev. gradients, etc 

import numpy as np
import jigsawpy
import netCDF4 as nc

from skimage.filters import gaussian, median
from skimage.measure import label, regionprops_table
from skimage.morphology import disk

def form_land_mask_connect(elev, edry=1):

    print("Forming connected land mask...")

    mask = label(elev>edry, background=1)
    prop = regionprops_table(
        mask, properties=["area", "label"])

    imax = np.argmax(prop["area"])
    
    land = np.zeros(mask.shape, dtype=np.uint8)
    land[mask == prop["label"][imax]] = 0
    land[mask != prop["label"][imax]] = 1

    return land


def setup_shoreline_pixels(hmat, land, hval):

    print("Computing shore adj. h(x)...")

    epos = np.logical_and.reduce((
        land[+1:-1, +1:-1]>=1, land[+1:-1, +2:]==0))
    wpos = np.logical_and.reduce((
        land[+1:-1, +1:-1]>=1, land[+1:-1, :-2]==0))
    npos = np.logical_and.reduce((
        land[+1:-1, +1:-1]>=1, land[:-2, +1:-1]==0))
    spos = np.logical_and.reduce((
        land[+1:-1, +1:-1]>=1, land[+2:, +1:-1]==0))

    mask = np.full(hmat.shape, False, dtype=bool)
    mask[+1:-1, +1:-1] = \
        np.logical_or.reduce((npos, epos, spos, wpos))

    hmat[mask] = hval * 1./1.  # mark shoreline min spacing

    return hmat


def coarsen_spacing_pixels(hmat, down):

    print("Coarsening mesh-spacing pixels...")

    rows = hmat.shape[0] // down
    cols = hmat.shape[1] // down

    htmp = np.full(
        (rows, cols), (np.amax(hmat)), dtype=hmat.dtype)

    for jpos in range(down):
        for ipos in range(down):

            iend = hmat.shape[0] - down + ipos + 1
            jend = hmat.shape[1] - down + jpos + 1

            htmp = np.minimum(
                htmp,
            hmat[ipos:iend:down, jpos:jend:down])

    return htmp


def filter_pixels_harmonic(hmat, exp=1):

    filt = remap_pixels_to_corner(hmat, exp)
    filt = remap_corner_to_pixels(filt, exp)

    return filt


def remap_pixels_to_corner(hmat, exp=1):

    R = hmat.shape[0]; C = hmat.shape[1]

    npos = np.arange(+0, hmat.shape[0] + 1)
    epos = np.arange(-1, hmat.shape[1] - 0)
    spos = np.arange(-1, hmat.shape[0] - 0)
    wpos = np.arange(+0, hmat.shape[1] + 1)
    
    npos[npos >= +R] = R - 1; spos[spos <= -1] = +0
    epos[epos <= -1] = C - 1; wpos[wpos >= +C] = +0

    npos, epos = np.meshgrid(
        npos, epos, sparse=True, indexing="ij")
    spos, wpos = np.meshgrid(
        spos, wpos, sparse=True, indexing="ij")

    htmp = (1. / hmat) ** exp
    hinv = htmp[npos, epos] + \
           htmp[npos, wpos] + \
           htmp[spos, epos] + \
           htmp[spos, wpos]

    return (4. / hinv) ** (1.0 / exp)


def remap_corner_to_pixels(hmat, exp=1):

    R = hmat.shape[0]; C = hmat.shape[1]

    npos = np.arange(+1, hmat.shape[0] + 0)
    epos = np.arange(+0, hmat.shape[1] - 1)
    spos = np.arange(+0, hmat.shape[0] - 1)
    wpos = np.arange(+1, hmat.shape[1] + 0)
    
    npos[npos >= +R] = R - 1; spos[spos <= -1] = +0
    epos[epos <= -1] = C - 1; wpos[wpos >= +C] = +0

    npos, epos = np.meshgrid(
        npos, epos, sparse=True, indexing="ij")
    spos, wpos = np.meshgrid(
        spos, wpos, sparse=True, indexing="ij")

    htmp = (1. / hmat) ** exp
    hinv = htmp[npos, epos] + \
           htmp[npos, wpos] + \
           htmp[spos, epos] + \
           htmp[spos, wpos]

    return (4. / hinv) ** (1.0 / exp)
    
    
def swe_wavelength_spacing(
        elev, land, nwav, hmin, hmax, grav=9.80665,
        T_M2=12.42*60.*60.):

    print("Computing wavelength heuristic...")

    vals = np.maximum(1, -elev)
    vals = T_M2 * (grav * vals) ** (1./2.) / nwav / 1000.

    vals[np.logical_and(elev >= -4., elev <= 4.)] = hmin

    vals = np.maximum(vals, hmin)
    vals = np.minimum(vals, hmax)

    vals = np.asarray(vals, dtype=np.float32)

    return vals


def elev_sharpness_spacing(
        xlon, ylat, 
        elev, land, nslp, hmin, hmax):

    print("Computing GRAD(elev) heuristic...")
    meters_per_degree = 111132.92
    dx = np.mean(np.diff(xlon)) * meters_per_degree
    dy = np.mean(np.diff(ylat)) * meters_per_degree

    by, bx = _earth_gradient(elev, dx, dy)
    bs = np.sqrt(bx**2 + by**2)  # get overall slope

    # Calculating the slope function
    eps = 1e-10  # small number to approximate derivative
    dp = np.clip(-elev, None, -50)
    vals = (2 * np.pi / nslp) * np.abs(dp) / (bs + eps)
    vals[vals > hmax] = hmax
    vals[vals < hmin] = hmin

    return vals


def _earth_gradient(F, dx, dy):
    """
    earth_gradient(F,HX,HY), where F is 2-D, uses the spacing
    specified by HX and HY. HX and HY can either be scalars to specify
    the spacing between coordinates or vectors to specify the
    coordinates of the points.  If HX and HY are vectors, their length
    must match the corresponding dimension of F.
    """
    Fy, Fx = np.zeros(F.shape), np.zeros(F.shape)

    # Forward diferences on edges
    Fx[:, 0] = (F[:, 1] - F[:, 0]) / dx
    Fx[:, -1] = (F[:, -1] - F[:, -2]) / dx
    Fy[0, :] = (F[1, :] - F[0, :]) / dy
    Fy[-1, :] = (F[-1, :] - F[-2, :]) / dy

    # Central Differences on interior
    Fx[:, 1:-1] = (F[:, 2:] - F[:, :-2]) / (2 * dx)
    Fy[1:-1, :] = (F[2:, :] - F[:-2, :]) / (2 * dy)

    return Fy, Fx


def scale_spacing_via_mask(args, vals):
    print("User-defined h(x) scaling...")
    if hasattr(args, 'mask_file') and args.mask_file:  # Check if mask_file exists and is not empty
        data = nc.Dataset(args.mask_file, "r")
        vals = np.asarray(data["val"][:], dtype=np.float32)
    else:
        print("No mask file provided. Scaling will not be applied.")
    return vals

