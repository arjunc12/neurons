'''
Module for making bash files to run scripts on supercomputer.
Automates writing and naming of files, which especially helps with naming the
files correctly 
'''

import argparse
import pylab
import os

def run_simulations_snider(rmin, rmax, lmin, lmax, rstep, lstep, num_iters):
    for r in pylab.arange(rmin, rmax, rstep):
        for l in pylab.arange(lmin, lmax, lstep):
            print 'r = %f, l = %f' % (r, l)
            command = 'python neuron_builder.py -a snider -rp %f -l %f' % (r, l)
            print command
            os.system(command)

def main():
    parser = argparse.ArgumentParser()

    # ranges
    parser.add_argument('-rmin', type=float, default=0.05)
    parser.add_argument('-rmax', type=float, default=2)
    parser.add_argument('-lmin', type=float, default=0.05)
    parser.add_argument('-lmax', type=float, default=2)

    # increments
    parser.add_argument('-rstep', type=float, default=0.05)
    parser.add_argument('-lstep', type=float, default=0.05)

    # number of different processes to use, so number of trials to run for each
    # parameter combination
    parser.add_argument('-x', '--num_iters', type=int, default=50, dest='num_iters')

    args = parser.parse_args()
    print vars(args)

    rmin = args.rmin
    rmax = args.rmax
    lmin = args.lmin
    lmax = args.lmax

    rstep = args.rstep
    lstep = args.lstep

    num_iters = args.num_iters

    fname = 'neuron_builder_snider.sh'
    print fname
    f = open(fname, 'w')

    # create script variables for all the parameters
    f.write('algorithm=\'snider\'\n')
    f.write('\n')

    f.write('rmin=%f\n' % rmin)
    f.write('rmax=%f\n' % rmax)
    f.write('\n')

    f.write('lmin=%f\n' % lmin)
    f.write('lmax=%f\n' % lmax)
    f.write('\n')

    f.write('rstep=%f\n' % lstep)
    f.write('lstep=%f\n' % rstep)
    f.write('\n')

    f.write('x=%d\n' % num_iters)
    f.write('\n')

    pyscript = 'neuron_builder.py'

    py_command = 'python %s -r $r -l $l' % pyscript

    # logic for looping through parameter range and running the script
    f.write('for r in $(seq $rmin $rstep $rmax); do\n')
    f.write('    for l in $(seq $lmin $lstep $lmax); do\n')
    f.write('        for iter in $(seq 1 1 $x); do\n')
    f.write('           ' + py_command)
    f.write('        done\n')
    f.write('        wait\n')
    f.write('    done\n')
    f.write('done\n')

    f.close()
    run_simulations_snider(rmin, rmax, lmin, lmax, rstep, lstep, num_iters)

if __name__ == '__main__':
    main()
