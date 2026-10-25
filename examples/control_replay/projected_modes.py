"""Frozen independent replay functions; provenance is in the manifest."""
import numpy as np

def projected_modes(matrix,xyz,mass):
    centered=xyz-np.average(xyz,axis=0,weights=mass)
    vectors=[(np.sqrt(mass)[:,None]*axis).ravel() for axis in np.eye(3)]
    vectors += [(np.sqrt(mass)[:,None]*np.cross(centered,axis)).ravel() for axis in np.eye(3)]
    rigid=[]
    for vector in vectors:
        residual=vector.copy()
        for _ in range(2):
            for unit in rigid:residual-=unit*np.dot(unit,residual)
        if np.linalg.norm(residual)>1e-10*max(1.,np.linalg.norm(vector)):
            rigid.append(residual/np.linalg.norm(residual))
    basis=np.array(rigid).T
    values,vectors=np.linalg.eigh(np.eye(matrix.shape[0])-basis@basis.T)
    internal=vectors[:,values>.5]
    roots=np.repeat(np.sqrt(mass*(1.660539040e-27/9.10938356e-31)),3)
    weighted=(matrix+matrix.T)/(2*np.outer(roots,roots))
    roots=np.linalg.eigvalsh(internal.T@weighted@internal)
    return np.sign(roots)*np.sqrt(abs(roots))*219474.63067,len(rigid)
