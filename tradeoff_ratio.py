import pandas as pd
import os

FRONTS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic'

def tradeoff_ratio(mcost, opt_mcost, scost, opt_scost):
    return (mcost / opt_mcost) / (scost / opt_scost)

def write_ratios():
    bad_fronts = open('bad_fronts.csv', 'w')
    with open('/iblsn/data/Arjun/neurons/pareto_steiner_output/tradeoff_ratio.csv', 'w') as tradeoff_file:
        tradeoff_file.write('neuron name, neuron type, tradeoff ratio\n')
        for neuron_name in os.listdir(FRONTS_DIR):
            for neuron_type in os.listdir(FRONTS_DIR + '/' + neuron_name):
                front_dir = FRONTS_DIR + '/' + neuron_name + '/' + neuron_type
                if 'pareto_front.csv' not in os.listdir(front_dir):
                    continue
                if 'tree_costs.csv' not in os.listdir(front_dir):
                    continue
                front_file = '%s/pareto_front.csv' % front_dir
                costs_file = '%s/tree_costs.csv' % front_dir
                pareto_front = pd.read_csv(front_file, skipinitialspace=True)
                mcosts, scosts = pareto_front['mcost'], pareto_front['scost']
                opt_mcost, opt_scost = min(mcosts), min(scosts)

                costs = pd.read_csv(costs_file, skipinitialspace=True)
                try:
                    trees = costs['tree']
                except:
                    bad_fronts.write('%s, %s\n' % (neuron_name, neuron_type))
                    print "bad front", neuron_name, neuron_type
                    continue
                costs = costs[costs['tree'] == 'neural']
                mcost = list(costs['mcost'])[0]
                scost = list(costs['scost'])[0]

                ratio = tradeoff_ratio(mcost, opt_mcost, scost, opt_scost)
                neuron_type = neuron_type.replace('_', ' ')
                tradeoff_file.write('%s, %s, %f\n' % (neuron_name, neuron_type, ratio))

    bad_fronts.close()

if __name__ == '__main__':
    write_ratios()
