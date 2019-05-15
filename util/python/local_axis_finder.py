import numpy as np
import numpy.linalg as LA

# This code is a tool that I used to find the SU(2) local axis associated
# with the SO(3) local axis (spatial local axis usually used in respack)
# The way to make it work is to define a vector n and an angle theta
# such that the SO(3) rotation matrix "R" resulting is the rotation
# that would be required to turn the local axis to the cartesian axis x-y-z.
# In the end, the local axis in real space is printed followed by the 
# local spinor axis necessary for the wannier calculation.
#



#####################################
#input parameters:
theta = 2*np.pi/3.0 #+ 1*2*np.pi/12.0 
nx = 0.
ny = 0.
nz = 1.
#####################################

n = np.array([nx, ny, nz])

n_unit =  n / (np.sqrt(np.dot(n,n)))
#print np.dot(n,n)
#print np.dot(n_unit,n_unit)

norm = 0.0+ np.sqrt(np.dot(n,n))
nx = nx/norm
ny = ny/norm
nz = nz/norm

cos = np.cos(theta)
sin = np.sin(theta)

nx2= nx*nx
ny2= ny*ny
nz2= nz*nz

#print nx, ny, nz
#print nx,ny,nz,  nx2+ny2+nz2

#print 'cos' , cos
#print 'sin' , sin


#find matrix R

R0 = np.array([[ nx2 +(ny2+nz2)*cos, (1.-cos)*ny*nx, (1.-cos)*nz*nx],
               [ (1.-cos)*nx*ny, ny2 +(nx2+nz2)*cos, (1.-cos)*nz*ny],
               [ (1.-cos)*nx*nz, (1.-cos)*ny*nz, nz2 +(ny2+nx2)*cos]])

R1 = np.array([[       0, -nz*sin,  ny*sin],
               [  nz*sin,       0, -nx*sin],
               [ -ny*sin,  nx*sin,       0]])

R = R0 + R1


#find matrix U

cos12 = np.cos(theta/2.0)
sin12 = np.sin(theta/2.0)

U = np.array([[ cos12-1.0j*nz*sin12, -(ny+1.0j*nx)*sin12],
              [  (ny-1.0j*nx)*sin12, cos12+1.0j*nz*sin12]])

print "\nLine to be appended to determine local axis in input.in:\n"

for ii in range(3):
  for jj in range(3):
    print '% 2.5f' % R[ii,jj],

print '  ',
for ii in range(2):
  for jj in range(2):
    print '(% 2.5f,% 2.5f) ' % (U[ii,jj].real,U[ii,jj].imag),


print "\n\nDo not forget to put the the flag FLG_INITIAL_GUESS_DIREC = 2"
print "in input.in to use the SU(2) part. Otherwise, if"
print "FLG_INITIAL_GUESS_DIREC = 1, only the SO(3) will be considered.\n"




