from optparse import OptionParser

import pylab
import numpy

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-x", "--xlabel", help="label for x axis")
    parser.add_option("-y", "--ylabel", help="label for y axis")
    parser.add_option("-t", "--title", help="title for plot")
    parser.add_option("-f", "--file", help="timing file")
 
    (options, args) = parser.parse_args()
    timings = []
    with open(options.file) as f:
        for l in f:
            lc = l.rstrip().split(' ')
            timings.append([float(i) for i in lc])

    timings = numpy.array(timings)
    pylab.figure()
    pylab.plot(timings[:,0], timings[:,1]) 
    pylab.xlabel(options.xlabel)
    pylab.ylabel(options.ylabel)
    pylab.title(options.title)
    ax = pylab.axes()
    ax.set_axis_bgcolor((1.0, 1.0, 1.0))
    pylab.show()
