#compute cluster nd realization mean and sd
#AG
#test cases


import _pickle as pickle
import numpy as np
import netCDF4





#data = np.load('lagrangian_data.npz')
##print ('done')
#npart = data['Npart']
#mux = data['mux']
#muy = data['muy']
#dxdx = data['dxdx_sum']
#dxdy = data['dxdy_sum']
#dydy = data['dydy_sum']
#drdr = data['drdr_sum']


#print(npart.shape)
#print(mux.shape)
#print(mux[:,1])


###############################

def spherical_bearing(phi1, phi2, lam1, lam2): #{{{
  """
  compute the spherical bearing
      http://www.movable-type.co.uk/scripts/latlong.html

      where lambda is longitude and phi is latitude

      spherical_bearing returned is from N (-pi, pi)

      Phillip Wolfram
      LANL
      08/07/2014
  """
  dphi = phi2 - phi1
  dlam = lam2 - lam1

  return np.arctan2(np.sin(dlam)*np.cos(phi2), np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dlam)) #}}}



def normalized_haversine_formula(phi1, phi2, lam1, lam2):  #{{{
    """
    compute the distance between two points via Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html

    where lambda is longitude and phi is latitude

    c returned in non-dimensional units (radians)

    Phillip Wolfram
    LANL
    07/18/2014
    """
    dphi = phi2 - phi1
    dlam = lam2 - lam1

    a = np.sin(dphi/2.0)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam/2.0)**2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0-a))

    return c # distance is c * Radius_of_earth #}}}

def compute_com_and_dispersion(plat, plon, r, indexOnCluster, maxIndex, indexToParticle): #{{{
    """
    compute particle disperison (2nd moment) for cluster, using lat / lon for basis of
    calculation to set coordinates for Kxx,Kxy,Kyy where x is zonal and y meridional

    dispersion units in m^2

    Phillip Wolfram
    LANL
    07/18/2014
    """

    N = maxIndex.shape[0]
    clat = np.zeros(N)
    clon = np.zeros(N)
    dxdx = np.zeros(N)
    dxdy = np.zeros(N)
    dydy = np.zeros(N)
    drdr = np.zeros(N)
    Npart = np.zeros(N)

    print ('Computing COM and dispersion for ', N, ' clusters')
    for aCluster, maxInd in enumerate(maxIndex):
        # get points of cluster
        particles = indexToParticle[indexOnCluster[aCluster,0:maxInd]]
        pslat = plat[particles]
        pslon = plon[particles]

        # compute center of mass
        clat[aCluster] = np.mean(pslat)
        clon[aCluster] = np.mean(pslon)

        # compute distances in m from lat / long  (use COM for simplicity, but there will be some error because of coordinate transform)
        dx = r * normalized_haversine_formula(clat[aCluster], clat[aCluster], clon[aCluster], pslon)
        dy = r * normalized_haversine_formula(clat[aCluster], pslat,          clon[aCluster], clon[aCluster])
        dr = r * normalized_haversine_formula(clat[aCluster], pslat,          clon[aCluster], pslon)

        # fix orientation of points
        bearing = spherical_bearing(clat[aCluster], pslat, clon[aCluster], pslon)
        # because arctan2 returns results from -pi to pi for bearing, flip values to get right sign
        dx -= 2*dx*(bearing < 0)
        dy -= 2*dy*(np.fabs(bearing) > np.pi/2.0)

        # store values
        Nparticles = len(particles)
        dxdx[aCluster] = sum(dx*dx)
        dxdy[aCluster] = sum(dx*dy)
        dydy[aCluster] = sum(dy*dy)
        drdr[aCluster] = sum(dr*dr)
        Npart[aCluster] = Nparticles

    return clon, clat, dxdx, dxdy, dydy, drdr, Npart #}}}




def proj_lat_long(x, y, z):  #{{{
    """
    compute the latitude and longitude from
    the x,y,z points (follow's Doug's pnt.h)
    """
    plat = np.arcsin(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    plong = np.arctan2(y, x)
    return plat, plong  #}}}



indexToParticle = pickle.load(open("verticesToParticle.p","rb"))
indexOnCluster = pickle.load(open("verticesOnCell.p","rb"))
maxIndex = pickle.load(open("maxVertices.p","rb"))


fname_in="/Users/garanaik/Desktop/LIGHT/Data_and_plots/mpaso.hist.am.lagrPartTrack.0001-01-01_00.00.00.nc"
f_in = netCDF4.Dataset(fname_in, 'r')
rEarth = f_in.sphere_radius
Ntime = len(f_in.dimensions['Time'])
Nclusters = maxIndex.shape[0]

mux = np.zeros((Ntime,Nclusters))
muy = np.zeros((Ntime,Nclusters))
dxdx_sum = np.zeros((Ntime,Nclusters))
dxdy_sum = np.zeros((Ntime,Nclusters))
dydy_sum = np.zeros((Ntime,Nclusters))
drdr_sum = np.zeros((Ntime,Nclusters))
Npart = np.zeros((Ntime,Nclusters))


print(indexToParticle.shape)

print(indexOnCluster.shape)
print(maxIndex[:])



for t in np.arange(Ntime):
    print ('Looping over step ', t, ' of ', Ntime)
    x = f_in.variables['xParticle'][t]
    y = f_in.variables['yParticle'][t]
    z = f_in.variables['zParticle'][t]
    plat, plon = proj_lat_long(x,y,z)
    #print(x,y,z,plat,plon)
 
    # compute center of mass and relative dispersion to this center of mass
    mux[t,:], muy[t,:], dxdx_sum[t,:], dxdy_sum[t,:], dydy_sum[t,:], drdr_sum[t,:], Npart[t,:] = \
      compute_com_and_dispersion(plat, plon, rEarth, indexOnCluster, maxIndex, indexToParticle)


print ('done')
