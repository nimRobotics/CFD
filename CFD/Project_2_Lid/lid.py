import numpy
from matplotlib import pyplot
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from math import pi
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def plot_3D(x, y, p):
    '''Creates 3D plot with appropriate limits and viewing angle
    Parameters:
    ----------
    x: array of float
        nodal coordinates in x
    y: array of float
        nodal coordinates in y
    p: 2D array of float
        calculated potential field
    '''
    fig = pyplot.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = numpy.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    ax.view_init(30,45)

def L2_error(p, pn):
	return numpy.sqrt(numpy.sum((p - pn)**2)/numpy.sum(pn**2))

def L1norm(new, old):
    norm = numpy.sum(numpy.abs(new-old))
    return norm

def initialise(nx, ny, xmax, xmin, ymax, ymin):
    '''Initialize the Poisson problem initial guess and other variables
    Parameters:
    ----------
    nx : int
        number of mesh points in x
    ny : int
        number of mesh points in y
    xmax: float
        maximum value of x in mesh
    xmin: float
        minimum value of x in mesh
    ymax: float
        maximum value of y in mesh
    ymin: float
        minimum value of y in mesh
    Returns:
    -------
    X  : 2D array of floats
        X-position of mesh
    Y  : 2D array of floats
        Y-position of mesh
    p_i: 2D array of floats
        initial guess of p
    dx : float
        mesh size in x direction
    dy : float
        mesh size in y direction
    '''

    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)

    # Mesh
    x  = numpy.linspace(xmin,xmax,nx)
    y  = numpy.linspace(ymin,ymax,ny)
    X,Y = numpy.meshgrid(x,y)

    # Initialize
    p_i  = numpy.zeros((ny,nx))

    return X, Y, x, y, p_i, dx, dy

def calculate(omega, psi, dx, dy, l1_target):
    '''Performs Jacobi relaxation

    Parameters:
    ----------
    omega : 2D array of floats
        Initial guess for omega
    psi : 2D array of floats
        Initial guess for psi
    dx: float
        Mesh spacing in x direction
    dy: float
        Mesh spacing in y direction
    l1_target: float
        Target difference between two consecutive iterates

    Returns:
    -------
    omega: 2D array of float
        Distribution after relaxation
    '''

    l1_norm_omega = 1.
    l1_norm_psi = 1.
    iterations = 0

    nx, ny = omega.shape
    u_j = 1.

    while l1_norm_omega > l1_target or l1_norm_psi > l1_target:

        # calculate psi

        psi_d = psi.copy()

        psi[1:-1,1:-1] = 1./(2.*(dx**2 + dy**2)) * \
                         ((psi_d[1:-1,2:]+psi_d[1:-1,:-2])*dy**2 +\
                         (psi_d[2:,1:-1] + psi_d[:-2,1:-1])*dx**2 -\
                         -omega[1:-1,1:-1]*dx**2*dy**2)

        # BCs for psi are automatically enforced

        l1_norm_psi = L1norm(psi_d, psi)

        # calculate omega

        omega_d = omega.copy()

        omega[1:-1,1:-1] = .25 * (omega_d[1:-1,2:] + omega_d[1:-1, :-2] \
                       + omega_d[2:, 1:-1] + omega_d[:-2, 1:-1])

        factor = (-1./(2.*(dy**2)));

        ##Neumann B.C.
        omega[-1, 1:-1] = (factor*((8.*psi[-2, 1:-1]) - psi[-3, 1:-1])) - ((3. * u_j) / dy)
        # bottom
        omega[0, 1:-1] = (factor*((8.*psi[1, 1:-1]) - psi[2, 1:-1]))
        # right
        omega[1:-1, -1] = (factor*((8.*psi[1:-1, -2]) - psi[1:-1, -3]))
        # left
        omega[1:-1, 0] = (factor*((8.*psi[1:-1, 1]) - psi[1:-1, 2]))

        l1_norm_omega = L1norm(omega_d, omega)

        iterations += 1

    print('Number of Jacobi iterations: {0:d}'.format(iterations))
    return omega, psi

##variable declarations
nx = 41
ny = 41

xmin = 0.
xmax = 1.
ymin = 0.
ymax = 1.
l1_target = 1e-6
X_omega, Y_omega, x_omega, y_omega, omega_i, dx, dy = initialise(nx, ny, xmax, xmin, ymax, ymin)
X_psi, Y_psi, x_psi, y_psi, psi_i, dx, dy = initialise(nx, ny, xmax, xmin, ymax, ymin)
#plot_3D(x, y, p_i)

omega, psi = calculate(omega_i.copy(), psi_i.copy(), dx, dy, l1_target)

# plot_3D(x_psi, y_psi, psi)
pyplot.contourf(x_psi,y_psi,psi,20, cmap=cm.viridis)

maxPsi = numpy.absolute(psi).max()

print('max psi: {0:f}'.format(maxPsi))

maxOmega = numpy.absolute(omega).max()

print('max omega: {0:f}'.format(maxOmega))

pts = numpy.round(psi[32,::8], 4)
print(pts)

pyplot.show()
