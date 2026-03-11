import os, logging, argparse
import json
import numpy as np

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')

logger = logging.getLogger('json2eqdsk')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

awd = os.getenv('AWD')
if awd is None:
    awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

flt = np.float32
coco_dpsi_sign = [1, 1, -1, -1, 1, 1, -1, -1]


def scatter_to_rectangular(r, z, data, Rmesh, Zmesh, fill=np.nan):
    """
    Interpolate scattered data points to rectangular grid
    :param r, y: r, y coordinates of data points
    :param data1, data2: input functions, shape as r and z
    :param R, Z: meshgrid of cartesian output grid
    :return: R, Z, interpolated_data (and cache if return_cache)
    """
    import scipy

    cache = scipy.spatial.Delaunay(np.vstack((r, z)).T)
    interpolant = scipy.interpolate.CloughTocher2DInterpolator(cache, data, fill_value=fill)
    interp_data = np.reshape(interpolant(np.vstack((Rmesh.flat, Zmesh.flat)).T), Rmesh.shape).T

    return interp_data


def json2eqdsk(f_json, nR=129, nZ=257, cocos_out=7):

    with open(f_json, 'r') as fjson:
        json_d = json.load(fjson)

    geq = {}
    geq['CASE2'] = 'ASTRA'
    
    ip_sgn = json_d['internal']['SGNIP']['data']
    bt_sgn = json_d['internal']['SGNBT']['data']
    dpsi_sign = coco_dpsi_sign[cocos_out-1]
    geq['CURRENT'] = 1e6*json_d['variables']['IPL']['data']
    geq['BCENTR'] = json_d['variables']['BTOR']['data']
    
# Contours
    rsurf = np.array(json_d['equil']['r']['data'], dtype=np.float32)
    zsurf = np.array( json_d['equil']['z']['data'], dtype=np.float32)
    geq['RMAXIS'] = rsurf[0, 0]
    geq['ZMAXIS'] = zsurf[0, 0]
    geq['RBBBS'] = rsurf[-1, :]
    geq['ZBBBS'] = zsurf[-1, :]

    n_rho, n_the = rsurf.shape

    Rmin = np.min(rsurf[-1, :]) - 0.03
    Rmax = np.max(rsurf[-1, :]) + 0.03
    Zmin = np.min(zsurf[-1, :]) - 0.03
    Zmax = np.max(zsurf[-1, :]) + 0.03
    geq['Rgrid'] = np.linspace(Rmin, Rmax, nR)
    geq['Zgrid'] = np.linspace(Zmin, Zmax, nZ)

    geq['RDIM'] = Rmax - Rmin
    geq['ZDIM'] = Zmax - Zmin
    geq['RCENTR'] = json_d['variables']['RTOR']['data']
    geq['RLEFT'] = Rmin
    geq['ZMID'] = 0.5*(Zmax + Zmin)
    geq['NW'] = nR
    geq['NH'] = nZ
    
# 1d profiles
    f_dia   = np.array(json_d['equil']['f_dia']['data']   , dtype=flt)
    ffprime = np.array(json_d['equil']['ffprime']['data'] , dtype=flt)
    pprime  = np.array(json_d['equil']['pprime']['data']  , dtype=flt)
    pres    = np.array(json_d['equil']['pressure']['data'], dtype=flt)
    psi     = np.array(json_d['equil']['psi']['data']     , dtype=flt)
    q       = np.array(json_d['equil']['q']['data']       , dtype=flt)

    geq['SIMAG'] = ip_sgn*dpsi_sign*psi[0]
    geq['SIBRY'] = ip_sgn*dpsi_sign*psi[-1]

# Interpolate on nR-grid
    psin = (psi - psi[0])/(psi[-1] - psi[0])
    psin_rect = np.linspace(0, 1, nR)

    geq['FPOL']   = np.interp(psin_rect, psin, f_dia)
    geq['FFPRIM'] = np.interp(psin_rect, psin, ffprime)
    geq['PRES']   = np.interp(psin_rect, psin, pres)
    geq['PPRIME'] = np.interp(psin_rect, psin, pprime)
    geq['QPSI']   = np.interp(psin_rect, psin, q)
    
# Psi(rho, theta) -> Psi(R, Z)
    X, Y = np.meshgrid(geq['Rgrid'], geq['Zgrid'])
    pf_in = np.repeat(psi, n_the)
    geq['PSIRZ'] = scatter_to_rectangular(rsurf.flat, zsurf.flat, pf_in.flat, X, Y)
    
    return geq


if __name__ == '__main__':

    import eqdsk
    import matplotlib.pylab as plt

    parser = argparse.ArgumentParser(description='Concatenate json files to 1 NetCDF along time')
    parser.add_argument('-e', '--expequ', help='<exp><equ>', required=True)
    parser.add_argument('-t', '--t_id', help='time ID', required=False, default=1)
    args = parser.parse_args()

    f_json = '%s/ncdf_out/%s-%s.json' %(awd, args.expequ, args.t_id)
    geq = json2eqdsk(f_json)

    f_eqdsk = '%s/ncdf_out/%s-%s.eqdsk' %(awd, args.expequ, args.t_id)
    eqd = eqdsk.EQDSK()
    eqd.write(f_eqdsk, geq=geq)
    eqd.plot()
    plt.show()
