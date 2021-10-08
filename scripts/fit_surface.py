import numpy as np
import scipy as sp
import scipy.interpolate as interp
import pylab
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import argparse, sys


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:", \
                                       ["help", "input=", "output="])
    except getopt.GetoptError, err:
        print(str(err))
        usage()
        sys.exit(2)
    output = None
    input = None
    

if __name__ == '__main__':

    # Declare and parse command line arguments
    parser = argparse.ArgumentParser(\
        description= 'Fit a bivariate spline to the force grid')
    parser.add_argument(\
        'input', metavar='input', help='The force grid data')
    parser.add_argument(\
        'output', metavar='output', help='The output file')
    args = parser.parse_args()

    # Get the input filename
    infile = args.input 
    
    # Load the data
    Di = np.loadtxt(infile)
    dtype = [('y',float),('z',float),('fy',float),('fz',float)]
    values = [ (Di[i,0],Di[i,1],Di[i,2],Di[i,3]) for i in xrange(Di.shape[0]) ]
    D = np.array(values, dtype=dtype)
    D = np.sort(D, order=['y','z'])
    D = np.array( [ [D[i][0],D[i][1],D[i][2],D[i][3]] for i in xrange(D.shape[0]) ])
    print(D)
    D = D.reshape((81,113,4))
    print(D)

    fig = pylab.figure()
    ax = Axes3D(fig)
    fig2 = pylab.figure()
    ax2 = Axes3D(fig2)

    y, z, fy, fz = D[:,:,0], D[:,:,1],  D[:,:,2],  D[:,:,3]


    #bs = interp.SmoothBivariateSpline(z.T.flatten(), y.T.flatten(), fz.T.flatten())
    bs = interp.Rbf(y.T.flatten(), z.T.flatten(), fz.T.flatten(), epsilon=2.0, smooth=0.99)

    # print("{0} Bivariate Spline Knots (X): \n".format(len(bs.get_knots()[0])))
    # print(bs.get_knots()[0])

    # print("{0} Bivariate Spline Knots (Y): \n".format(len(bs.get_knots()[1])))
    # print(bs.get_knots()[1])

    # print("{0} Bivariate Spline Coeffs: \n".format(len(bs.get_coeffs().flatten())))
    # print(bs.get_coeffs())

    # print("Residual Error = {0}".format(bs.get_residual()))

    ax.plot_surface(y.T, z.T, fz.T, cmap=cm.jet)
    fig.add_axes(ax)
        
    fya = np.zeros(y.flatten().shape[0])
    for i, di in enumerate(zip(z.flatten(),y.flatten())):

        fya[i] = bs(di[0],di[1])#fy.flatten()[i]#bs.ev(di[0],di[1])
#        print("{0} {1} {2}\n".format(di[0], di[1], fya[i]))
        if( np.fabs(fya[i]) > 1000):
            print("{0} {1} {2}\n".format(di[0], di[1], fya[i]))
                                     
    print("Z is {0}\n".format(z.flatten().shape))
    print("Y is {0}\n".format(y.flatten().shape))
    print("FYA is {0}\n".format(fya.shape))
    
    ax2.plot_surface(y.T, z.T, fya.reshape(z.shape).T, cmap=cm.jet)
    fig2.add_axes(ax2)
    pylab.show()

    print(bs)
    print(y.T)
    print(z.T)
    print(fy.T)
