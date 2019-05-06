import numpy as np

from .h5py_gadget import gather_data
import h5py as h5
import ntpath

class Snapdata:
    S = None
    #particle data:
    npart = None
    time = None
    p_mass = None
    ID = None
    x = None
    y = None
    z = None
    vx = None
    vy = None
    vz = None
    v = None
    vr = None
    vphi = None
    p = None
    ax = None
    ay = None
    az = None
    a = None
    ax_p = None
    ay_p = None
    az_p = None
    a_p = None
    dt = None
    r = None
    rp = None
    ra = None
    E = None
    L = None
    Tr = None
    #potential approximation stuff:
    initialized_p_fit = None
    interp_x = None
    interp_y = None
    p_interpolator = None
    def __init__(self, fname, nfiles=1, rp_ra_fname=None,  Tr_fname=None):
        if fname[-5:] == ".hdf5":
            f = h5.File(fname, "r")
            self.npart = f.get("Header").attrs.get("NumPart_Total")[1]
            self.time = f.get("Header").attrs.get("Time")
            self.p_mass = f.get("Header").attrs.get("MassTable")[1]
            f.close()
            path, filename = ntpath.split(fname)
            filename = filename[:-5]
            self.ID = gather_data(path, filename, "PartType1/ParticleIDs")
            coords = gather_data(path, filename, "PartType1/Coordinates")
            self.x = coords[:,0]
            self.y = coords[:,1]
            self.z = coords[:,2]
            vels = gather_data(path, filename, "PartType1/Velocities")
            self.vx = vels[:,0]
            self.vy = vels[:,1]
            self.vz = vels[:,2]
            self.p = gather_data(path, filename, "PartType1/Potential")
            acc = gather_data(path, filename, "PartType1/Acceleration")
            self.ax = acc[:,0]
            self.ay = acc[:,1]
            self.az = acc[:,2]
            acc_p = gather_data(path, filename, "PartType1/PertAccel")
            if (acc_p.shape == (0,)):
                self.ax_p = np.zeros(self.ax.shape)
                self.ay_p = np.zeros(self.ay.shape)
                self.az_p = np.zeros(self.az.shape)
            else:
                self.ax_p = acc_p[:,0]
                self.ay_p = acc_p[:,1]
                self.az_p = acc_p[:,2]
            self.dt = gather_data(path, filename, "PartType1/TimeStep")
        else:
            import snap
            self.S = snap.cvar.S
            self.S.load_snapshot(fname,nfiles)
            self.npart = self.S.np()
            self.time = self.S.t()
            self.p_mass = self.S.mass(0)
            self.ID = np.array([self.S.ID(i) for i in range(self.npart)])
            self.x = np.array([self.S.pos(i,0) for i in range(self.npart)])
            self.y = np.array([self.S.pos(i,1) for i in range(self.npart)])
            self.z = np.array([self.S.pos(i,2) for i in range(self.npart)])
            self.vx = np.array([self.S.vel(i,0) for i in range(self.npart)])
            self.vy = np.array([self.S.vel(i,1) for i in range(self.npart)])
            self.vz = np.array([self.S.vel(i,2) for i in range(self.npart)])
            self.p = np.array([self.S.pot(i) for i in range(self.npart)])
            self.S.unload_snapshot()
            self.ax = np.zeros(self.npart)
            self.ay = np.zeros(self.npart)
            self.az = np.zeros(self.npart)
            self.dt = np.zeros(self.npart)
        self.r = np.sqrt(np.power(self.x,2) + np.power(self.y,2) + np.power(self.z,2))
        self.v = np.sqrt(np.power(self.vx,2) + np.power(self.vy,2) + np.power(self.vz,2))
        self.a = np.sqrt(np.power(self.ax,2) + np.power(self.ay,2) + np.power(self.az,2))
        if self.ax_p is not None:
            self.a_p = np.sqrt(np.power(self.ax_p,2) + np.power(self.ay_p,2) + np.power(self.az_p,2))
        self.E = .5 * np.power(self.v,2) + self.p
        self.vr = (self.vx * self.x + self.vy * self.y + self.vz * self.z) / self.r
        self.vphi = np.sqrt(np.power(self.v,2) - np.power(self.vr,2))
        self.L = self.r * self.vphi
        self.rp = np.zeros(self.npart)
        self.ra = np.zeros(self.npart)
        self.Tr = np.zeros(self.npart)
        try:
            from slvars import loadvars
            if rp_ra_fname != None:
                [tmp_ID,tmp_r,tmp_rp,tmp_ra] = loadvars(rp_ra_fname)
                self.sort(self.ID)
                tmp_isort = np.argsort(tmp_ID, kind='quicksort')
                assert (tmp_ID[tmp_isort] == self.ID).all()
                self.ra = np.abs(tmp_ra[tmp_isort])
                self.rp = np.abs(tmp_rp[tmp_isort])
                del tmp_isort
                del tmp_ID
                del tmp_r
                del tmp_ra
                del tmp_rp
            if Tr_fname != None:
                [tmp_ID, tmp_Tr] = loadvars(Tr_fname)
                self.sort(self.ID)
                tmp_isort = np.argsort(tmp_ID, kind='quicksort')
                assert (tmp_ID[tmp_isort] == self.ID).all()
                self.Tr = tmp_Tr[tmp_isort]
                del tmp_isort
                del tmp_ID
                del tmp_Tr
        except AssertionError:
            print "ID mismatch loading extra data!"
            exit()
        self.initialized_p_fit = False
        self.initialize_p_fit()
        self.sort(self.ID)
    def sort(self, key):
        isort = np.argsort(key, kind='quicksort')
        self.ID = np.copy(self.ID[isort])
        self.x = np.copy(self.x[isort])
        self.y = np.copy(self.y[isort])
        self.z = np.copy(self.z[isort])
        self.vx = np.copy(self.vx[isort])
        self.vy = np.copy(self.vy[isort])
        self.vz = np.copy(self.vz[isort])
        self.v = np.copy(self.v[isort])
        self.vr = np.copy(self.vr[isort])
        self.vphi = np.copy(self.vphi[isort])
        self.p = np.copy(self.p[isort])
        self.r = np.copy(self.r[isort])
        self.rp = np.copy(self.rp[isort])
        self.ra = np.copy(self.ra[isort])
        self.E = np.copy(self.E[isort])
        self.L = np.copy(self.L[isort])
        self.Tr = np.copy(self.Tr[isort])
        self.ax = np.copy(self.ax[isort])
        self.ay = np.copy(self.ay[isort])
        self.az = np.copy(self.az[isort])
        self.a = np.copy(self.a[isort])
        if self.ax_p is not None:
            self.ax_p = np.copy(self.ax_p[isort])
        if self.ay_p is not None:
            self.ay_p = np.copy(self.ay_p[isort])
        if self.az_p is not None:
            self.az_p = np.copy(self.az_p[isort])
        if self.a_p is not None:
            self.a_p = np.copy(self.a_p[isort])
        self.dt = np.copy(self.dt[isort])
        del isort
        return
    def recentre(self):
        dx = np.sum((-1 * self.p) * self.x) / np.sum(-1 * self.p)
        dy = np.sum((-1 * self.p) * self.y) / np.sum(-1 * self.p)
        dz = np.sum((-1 * self.p) * self.z) / np.sum(-1 * self.p)
        self.x = self.x - dx
        self.y = self.y - dy
        self.z = self.z - dz
        self.r = np.sqrt(np.power(self.x,2) + np.power(self.y,2) + np.power(self.z,2))
        return
    def initialize_p_fit(self, N=1000):
        self.sort(self.r)
        rmin = np.min(self.r)
        rmax = np.max(self.r)
        r_bins = np.logspace(np.log10(rmin),np.log10(rmax),N)
        from scipy.stats import binned_statistic
        self.interp_y = binned_statistic(self.r,self.p,statistic='median',bins=r_bins)[0]
        self.interp_x = np.power(10, .5 * (np.log10(r_bins[1:]) + np.log10(r_bins[:-1])))
        self.interp_y = np.concatenate((np.array([self.p[0]]),self.interp_y,np.array([self.p[-1]])))
        self.interp_x = np.concatenate((np.array([rmin]),self.interp_x,np.array([rmax])))
        self.p_interpolator = lambda r: np.interp(r,self.interp_x[np.isfinite(self.interp_y)],self.interp_y[np.isfinite(self.interp_y)])
        self.initialized_p_fit = True
        return
    def p_fit(self, r):
        if self.initialized_p_fit:
            return self.p_interpolator(r)
        else:
            print "Uninitialized call to p_fit!"
            exit()
