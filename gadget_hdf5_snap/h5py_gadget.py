import h5py as h5
import numpy as np
import multiprocessing
import os.path

class H5Data():
    def __init__(self, path, fbase):
        self.__path = path
        self.__fbase = fbase
        self.__parts = self.__find_parts(self.__path, self.__fbase)
        #self.__nb_cpu = multiprocessing.cpu_count() - 1
        self.__nb_cpu = 1
    def __subitem(self, name, parts, output):
        catalyst = []
        for part in parts:
            try:
                f = h5.File(part, 'r')
                catalyst.append(f[name].value.copy())
                f.close()
            except:
                f = 0
        output.put(catalyst)
        return 0
    def __getitem__(self, name):
        if self.__nb_cpu > 1:
            try:
                parts_split = np.array_split(self.__parts, self.__nb_cpu)
                procs = []
                outputs = []
                for parts in parts_split:
                    outputs.append(multiprocessing.Queue())
                    procs.append(multiprocessing.Process(target=self.__subitem, args=(name,parts.tolist(), outputs[-1])))
                    procs[-1].start()
                items = []
                for output in outputs:
                    items += output.get()
                for p in procs:
                    p.join()
            except IOError:
                self.__nb_cpu = 1
                return self[name]
        else:
            items = []
            for part in self.__parts:
                try:
                    f = h5.File(part,'r')
                    items.append(f[name].value.copy())
                    f.close()
                except:
                    f = 0
        if not(items):
            raise KeyError("Unable to open object (Object '" + name + "' doesn't exist in file with path '" + self.__path + "' and basename '" + self.__fbase + "')")
        else:
            return np.concatenate(items)
    def get_parts(self):
        return self.__parts
    def __find_parts(self,path,fbase):
        if os.path.exists(path + '/' + fbase + '.hdf5'):
            return [path + '/' + fbase + '.hdf5']
        elif os.path.exists(path+'/'+fbase+'.0.hdf5'):
            fcount = 0
            retval = []
            while os.path.exists(path + '/' + fbase + '.' + str(fcount) + '.hdf5'):
                retval.append(path + '/' + fbase + '.' + str(fcount) + '.hdf5')
                fcount += 1
            return retval
        else:
            raise IOError("Unable to open file (File with path '" + path + "' and basename'" + fbase + "' doesn't exist)")

def gather_data(path, fbase, hpath):
    """
    path: path of simulation data (can be particle or group)
    fbase: filename of data file (omit ".X.hdf5" portion)
    hpath: path of data table to gather, e.g. "PartType1/ParticleIDs"
    """
    f = H5Data(path,fbase)
    retval = f[hpath]
    del f
    return retval

def gather_attr(path, fbase, hpath, attr):
    """
    path: path of simulation data (can be particle or group)
    fbase: filename of data file (omit ".X.hdf5" portion)
    hpath: path of data table to gather, e.g. "PartType1/ParticleIDs"
    attr: name of attribute to fetch
    """

    for fname in H5Data(path, fbase).get_parts():
        with h5.File(fname, 'r') as f:
            try:
                return f[hpath].attrs[attr]
            except KeyError:
                continue
    raise KeyError("Unable to open attribute (One of object '" + hpath + "' or attribute '" + attr + "' doesn't exist in file with path '" + path + "' and basename '" + fbase + "')")
