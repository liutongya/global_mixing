## This script is used for correcting the AVISO velocity

import xarray as xr
from GeoApps.GridUtils import add_MITgcm_missing_metrics
from xmitgcm import open_mdsdataset
import numpy as np
from GeoApps.DiagnosticMethods import Dynamics

def write_field(fname, data):
    print('wrote to file: ' + fname)
    fid = open(fname, "wb")
    data.tofile(fid)
    fid.close()

path = '/data/home/liutongya/mixing/gcm/run_vel/'
dset = open_mdsdataset(path, prefix=['Diag_stat'])

dset, grid = add_MITgcm_missing_metrics(dset, periodic='X', boundary={'Y':'extend'})

length = 723


for num in np.arange(length):

    # calculate divergence and vorticity
    ug = dset.UVEL.where(dset.UVEL!=0)[num].load() # select a single timestep
    vg = dset.VVEL.where(dset.VVEL!=0)[num].load() # select a single timestep

    dyn = Dynamics(dset, grid=grid, arakawa='C')

    div = dyn.divg(ug, vg).load() # calculate divergence on C-grid
    #vor = dyn.curl(ug, vg).load() # calculate vorticity  on C-grid

    # load the streamfunction and velocity potential data output from the offline model
    psi = dset.PsiVEL.where(dset.PsiVEL!=0)[num].load() # select a single timestep

    # direct method for rotational and divergent components
    vs, us = dyn.grad(psi/100) # 100 is the depth of offline model
    us = -us

    u_out = '/data/home/liutongya/mixing/velocity/corr_vel/uvelcorr.' + str(num).zfill(10) + '.data'
    v_out = '/data/home/liutongya/mixing/velocity/corr_vel/vvelcorr.' + str(num).zfill(10) + '.data'
    
    utmp = us.fillna(0) * dset.maskW
    vtmp = vs.fillna(0) * dset.maskS
    
    write_field(u_out, utmp.values.astype('>f4'))
    write_field(v_out, vtmp.values.astype('>f4'))