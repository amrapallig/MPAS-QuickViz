
import _pickle as pickle
import numpy as np
import netCDF4





data = np.load('lagrangian_data.npz')
print ('done')
npart = data['Npart']
mux = data['mux']
muy = data['muy']
dxdx = data['dxdx_sum']
dxdy = data['dxdy_sum']
dydy = data['dydy_sum']
drdr = data['drdr_sum']
print(data.files)


print("npart",npart.shape)
print(mux.shape)
print(mux[:,:])
print(npart[:,:])
print(dxdx[:,:])
print(dxdy[:,:])
print(dydy[:,:])
print(drdr[:,:])
