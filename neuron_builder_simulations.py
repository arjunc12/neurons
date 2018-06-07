import argparse
import pylab
import os

def run_simulations_snider(rpmin, rpmax,\
                           rrmin, rrmax,\
                           lmin, lmax,\
                           proot_min, proot_max,\
                           rpstep, rrstep, lstep, proot_step,\
                           num_iters):
    for rp in pylab.arange(rpmin, rpmax, rpstep):
        for rr in pylab.arange(rrmin, rrmax, rrstep):
            for l in pylab.arange(lmin, lmax, lstep):
                for proot in pylab.arange(proot_min, proot_max, proot_step):
                    print 'rp = %f, rr = %f, l = %f, proot = %f' % (rp, rr, l, proot)
                    command = 'python neuron_builder.py -a snider -rp %f -rr %f -l %f -proot %f' % (rp, rr, l, proot)
                    print command
                    os.system(command)

def main():
    parser = argparse.ArgumentParser()

    # ranges
    parser.add_argument('-rpmin', type=float, default=1)
    parser.add_argument('-rpmax', type=float, default=2)
    
    parser.add_argument('-rrmin', type=float, default=1)
    parser.add_argument('-rrmax', type=float, default=2)
    
    parser.add_argument('-proot_min', type=float, default=0.05)
    parser.add_argument('-proot_max', type=float, default=0.95)
    
    parser.add_argument('-lmin', type=float, default=1)
    parser.add_argument('-lmax', type=float, default=2)

    # increments
    parser.add_argument('-rpstep', type=float, default=0.01)
    parser.add_argument('-rrstep', type=float, default=0.01)
    parser.add_argument('-lstep', type=float, default=0.01)
    parser.add_argument('-proot_step', type=float, default=0.05)

    # number of different processes to use, so number of trials to run for each
    # parameter combination
    parser.add_argument('-x', '--num_iters', type=int, default=50, dest='num_iters')

    args = parser.parse_args()

    rpmin = args.rpmin
    rpmax = args.rpmax

    rrmin = args.rpmin
    rrmax = args.rrmax

    lmin = args.lmin
    lmax = args.lmax

    proot_min = args.proot_min
    proot_max = args.proot_max

    rpstep = args.rpstep
    rrstep = args.rrstep
    lstep = args.lstep
    proot_step = args.proot_step

    num_iters = args.num_iters

    fname = 'neuron_builder_snider.sh'
    f = open(fname, 'w')

    # create script variables for all the parameters
    f.write('algorithm=\'snider\'\n')
    f.write('\n')

    f.write('rpmin=%f\n' % rpmin)
    f.write('rpmax=%f\n' % rpmax)
    f.write('\n')

    f.write('rrmin=%f\n' % rrmin)
    f.write('rrmax=%f\n' % rrmax)
    f.write('\n')
    
    f.write('lmin=%f\n' % lmin)
    f.write('lmax=%f\n' % lmax)
    f.write('\n')

    f.write('prootmin=%f\n' % proot_min)
    f.write('prootmax=%f\n' % proot_max)
    f.write('\n')

    f.write('rpstep=%f\n' % rpstep)
    f.write('rrstep=%f\n' % rrstep)
    f.write('lstep=%f\n' % lstep)
    f.write('prootstep=%f\n' % proot_step)
    f.write('\n')

    f.write('x=%d\n' % num_iters)
    f.write('\n')

    pyscript = 'neuron_builder.py'

    py_command = 'python %s -rp $rp -rr $rr -l $l' % pyscript

    # logic for looping through parameter range and running the script
    f.write('for rp in $(seq $rpmin $rpstep $rpmax); do\n')
    f.write('    for rp in $(seq $rpmin $rpstep $rpmax); do\n')
    f.write('        for l in $(seq $lmin $lstep $lmax); do\n')
    f.write('            for l in $(seq $lmin $lstep $lmax); do\n')
    f.write('                for iter in $(seq 1 1 $x); do\n')
    f.write('                    ' + py_command + '\n')
    f.write('                done\n')
    f.write('                wait\n')
    f.write('            done\n')
    f.write('        done\n')
    f.write('    done\n')
    f.write('done')

    f.close()
    run_simulations_snider(rpmin, rpmax,\
                           rrmin, rrmax,\
                           lmin, lmax,\
                           proot_min, proot_max,\
                           rpstep, rrstep, lstep, proot_step,\
                           num_iters)

if __name__ == '__main__':
    main()
