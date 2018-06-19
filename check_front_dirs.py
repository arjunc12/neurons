import os

no_front = open('no_front.csv', 'w')
no_trees = open('no_trees.csv', 'w')
fronts_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic'
good_dirs = 0
for neuron_name in os.listdir(fronts_dir):
    name_dir = '%s/%s' % (fronts_dir, neuron_name)
    for neuron_type in os.listdir(name_dir):
        type_dir = '%s/%s' % (name_dir, neuron_type)
        files = os.listdir(type_dir)
        tree_file = 'tree_costs.csv'
        front_file = 'pareto_front.csv'

        if (tree_file in files) and (front_file not in files):
            no_front.write('%s, %s\n' % (neuron_name, neuron_type))
        elif (tree_file not in files) and (front_file in files):
            no_trees.write('%s, %s\n' % (neuron_name, neuron_type))
        elif (tree_file in files) and (front_file in files):
            good_dirs += 1

no_front.close()
no_trees.close()
print good_dirs
