import netCDF4
import oc
import datetime as dt
import matplotlib.pyplot as plt
import geopy
plt.style.use('ggplot')

##


year = 1996

mldc = []
uc = []
vc = []
timec = []
for nn in range(1,13):
    data = netCDF4.Dataset('/run/media/georg/TRANSCEND/EIS2/results/eis3_test/test' + str(nn).zfill(2) + '.nc')
    if nn==1:
        print(data.variables['latitude'][0])
        print(data.variables['longitude'][0])
    mldc += list(data.variables['mld'][:,0,0])
    uc += list(data.variables['u'][:,0,0])
    vc += list(data.variables['v'][:,0,0])
    timec += list(data.variables['time'][:])

comp = oc.PCKL('/run/media/georg/TRANSCEND/EIS2/results/mimoc_1996/' +
               '{:.2f}'.format(data.variables['latitude'][0]) + '/' +
               '{:.2f}'.format(data.variables['longitude'][0]) + '/' +
               '20160122_p1.pckl')
mldp = comp.mld

mld = np.empty(12)
for nn in range(1,13):
    data = netCDF4.Dataset('/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month{:02d}_regrid.nc'.format(nn))
    mld[nn-1] = data.variables['DEPTH_MIXED_LAYER'][178, 1024]
    if nn==1:
        print(data.variables['lat'][178])
        print(data.variables['lon'][1024])

ref = dt.datetime(year=1900, month=1, day=1, hour = 0, minute=0, second=0)
date = np.empty(len(timec)).astype(object)
for nn in range(0, len(timec)):
    date[nn] = ref + dt.timedelta(hours=float(timec[nn]))

mldate = np.empty(12).astype(object)
for nn in range(1, 13):
    mldate[nn-1] = dt.datetime(year = year, month = nn, day=15, hour =0, minute = 0 )

# plt.plot(date, mldc, label='c interpolation')
# plt.plot(comp.time, mldp, label='python interpolation')
# plt.plot(mldate, mld, label='original')
#
# plt.legend()
# plt.show()

fig, axes = plt.subplots(2,1)
axes[0].plot(date, uc)
axes[0].plot(comp.time, comp.u)
# axes[0].plot(comp.time, comp.taux)
axes[1].plot(date, vc)
axes[1].plot(comp.time, comp.v)

plt.show()

##
# test correlation scale

plot = 2

nmonth = 1
thres = 0
nw = 150
nt = 100
nd = 40

data = netCDF4.Dataset('/run/media/georg/TRANSCEND/EIS2/results/eis3_test/test' + str(nmonth).zfill(2) + '.nc')

ww = data.variables['w'][:]
lat = data.variables['latitude'][:]
lon = data.variables['longitude'][:]

data.close()

long, latg = np.meshgrid(lon, lat)

origin = (33, 330)

nx = np.argmin(np.abs(lon - origin[1]))
ny = np.argmin(np.abs(lat - origin[0]))

origin = ww[:, ny, nx]

RR = np.empty((ww.shape[1]-2, ww.shape[2]-2))

for nn in range(1, ww.shape[1]-1):
    for mm in range(1, ww.shape[2]-1):
        RR[nn-1, mm-1] = np.corrcoef(origin, ww[:,nn,mm])[1,0]

aux = np.ma.masked_where(np.abs(RR)<thres, RR)

if plot==2:
    plt.pcolormesh(lon[1:-1], lat[1:-1], aux, vmin=-1, vmax=1)
elif plot==1:
    plt.pcolormesh(lon, lat, ww.mean(axis=0))
elif plot==0:
    plt.pcolormesh(lon, lat, ww[nw])
# plt.xlim(0,ww.shape[1]-2)
# plt.ylim(0,ww.shape[2]-2)
plt.colorbar()
plt.close()

# calcualte spatial correlation function

delta = np.linspace(0, 1200, nd)
distances = np.empty(latg.shape)
dr = delta[1] - delta[0]
dmin = np.argmin(np.abs(delta-200))


#MoranISum = np.zeros(delta.shape)
corrlength = np.empty(ww[:,1:-1,1:-1].shape)

avg = np.mean(np.mean(ww[:,1:-1,1:-1], axis=1), axis=1)
var = np.array([np.sum((ww[nn,1:-1,1:-1] - avg[nn])**2) / ww[nt,1:-1,1:-1].ravel().shape[0] for nn in range(0, len(ww))])
##
for nx in range(1, latg.shape[1]-2):
#for nx in range(1, 2):
    print(nx)
    for ny in range(1, latg.shape[0]-2):
    #for ny in range(1, 2):

        for nn in range(1, latg.shape[0]-1):
            for mm in range(1, latg.shape[1]-1):
                distances[nn,mm]  = geopy.distance.great_circle((lat[ny], lon[nx]), (lat[nn], lon[mm])).km
        d_distances = np.round(distances/dr) * dr

        MoranI = np.zeros((len(ww), ) + delta.shape)
        # for nt in range(0, len(ww)):
        for nr in range(0, len(delta)):
            index = np.where(d_distances==delta[nr])

            MoranI[:,nr] = (ww[:,ny,nx] - avg)*np.sum(np.array([ww[nn][index] for nn in range(0, len(ww))]).T - avg, axis=0) / var / index[0].shape[0]
            # MoranISum[nr] += np.sum((ww[nt,ny,nx] - avg)*(ww[nt][index] - avg)) / var / index[0].shape[0]

        corrlength[:, ny, nx] = delta[np.argmax(MoranI[:,dmin:], axis=-1) + dmin]
        # del MoranI

##

fig, ax = plt.subplots(1,3, figsize=(18,6))

# image = ax[0].pcolormesh(lon, lat, ww[0])

# image = ax[0].pcolormesh(lon[1:-1], lat[1:-1], aux, vmin=-1, vmax=1)
# image1 = ax[1].pcolormesh(lon[1:-1], lat[1:-1], corrlength)
image0 = ax[0].pcolormesh(lon, lat, ww[100])
image1 = ax[1].pcolormesh(lon[1:-1], lat[1:-1], corrlength.mean(axis=0))
# image1 = ax[1].pcolormesh(lon[1:-1], lat[1:-1], corrlength[500])

ax[0].set_xlim(lon[0], lon[-1])
ax[1].set_xlim(lon[0], lon[-1])
ax[0].set_ylim(lat[-1], lat[0])
ax[1].set_ylim(lat[-1], lat[0])
# ax[0].contour(lon, lat, d_distances, levels=delta, colors='k')
plt.colorbar(image0, ax=ax[0])
plt.colorbar(image1, ax=ax[1])
#
# plt.plot(index[0], index[1])
for nn in range(0, len(ww)):
    ax[2].plot(delta[1:], MoranI[nn,1:])
ax[2].plot(delta[1:], MoranI.mean(axis=0)[1:], color='k', linewidth=2)

plt.show()

##

max = np.abs(ww[0]).max()

hist = scipy.histogram(ww[0,1:-1,1:-1], 100, range=(-max/10, max/10))

oc.plotStepH(hist[1], hist[0])
plt.show()
