import os
from sys import argv

cell_type, lab, species, loc = argv[1:]

old_dir = 'datasets/neuron_nmo/%s/CNG\\ version' % lab 
new_dir = 'datasets/%s/%s/%s/%s' % (cell_type, species, loc, lab)

mkdir_command = 'mkdir -p %s' % new_dir
print mkdir_command
os.system(mkdir_command)

mv_command = 'mv %s/*.swc %s/' % (old_dir, new_dir)
print mv_command
os.system(mv_command)
