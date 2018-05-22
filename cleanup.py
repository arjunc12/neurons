import pandas as pd
from dist_functions import pareto_dist_scale

missing = set([('knock-in-18', 'basal dendrite'), ('RatS1-15-80', 'basal dendrite'), ('5-HT1B-F-500017', 'truncated axon'), ('6-S18-3', 'apical dendrite'), ('d1s2-cell-1', 'basal dendrite'), ('MJ021506SL2', 'truncated axon'), ('45dd_11', 'basal dendrite'), ('Devel47-neuronH', 'apical dendrite'), ('6R-C3-S4', 'basal dendrite'), ('OT-IN-9028', 'basal dendrite'), ('18_2_1', 'apical dendrite'), ('2004-06-22-A-2', 'truncated axon'), ('TS072007s1', 'truncated axon'), ('AM61-1-2', 'axon'), ('Culture-60', 'apical dendrite'), ('C10-22', 'apical dendrite'), ('11-1-1-LN', 'basal dendrite'), ('Trh-M-000169', 'truncated axon'), ('CA2-P-P8-1', 'basal dendrite'), ('2-7-1-KD', 'apical dendrite'), ('107ML', 'basal dendrite'), ('13_10', 'basal dendrite'), ('13A2-Cell12', 'basal dendrite'), ('fru-F-200092', 'axon'), ('m2s2s4t-d-deep', 'basal dendrite'), ('veh15', 'basal dendrite'), ('F7-6', 'basal dendrite'), ('470364', 'basal dendrite'), ('MMP9KO-1V1-LayerIII-5', 'basal dendrite'), ('TML2', 'axon'), ('GFP-GT02-N02-all-gray', 'basal dendrite'), ('C88-27-03-13-C', 'truncated axon'), ('Pvalb-IRES-Cre_Ai14-199001-04-02-01_491119245_m', 'axon'), ('fru-M-300308', 'truncated axon'), ('10-2918-TTX48h-section5-left-cell4', 'basal dendrite'), ('185-2-4th', 'apical dendrite'), ('pcs3_2_3', 'apical dendrite'), ('O-2', 'basal dendrite'), ('MMP9KO-P7-15', 'basal dendrite'), ('inv-7-16-08-cell-1', 'basal dendrite'), ('L5SC-x140526c-cell-3', 'axon'), ('2303113_Disc1', 'axon'), ('C164-15-10-13-D', 'basal dendrite'), ('P7-Control-for-Mfn1-08-021', 'apical dendrite'), ('10-2918-TTX48h-section6-left-cell6', 'basal dendrite'), ('194-4-14nj', 'basal dendrite'), ('R2P20topneuN-1', 'basal dendrite'), ('fru-M-200345', 'truncated axon'), ('pcd16_3_3', 'basal dendrite'), ('13May2010', 'basal dendrite'), ('L23BPC-j131114c-cell-1', 'basal dendrite'), ('10-2912-TTX8h-slide3-section3-right-cell8', 'basal dendrite'), ('15_0', 'basal dendrite'), ('Z-2', 'basal dendrite'), ('DAS303F', 'truncated axon'), ('207-1-36LL', 'basal dendrite'), ('Badea2012Fig6A-C-L', 'axon'), ('Pvalb-IRES-Cre_Ai14-177452-03-02-01_491119336_m', 'axon'), ('KO-mPFC-B-20X-2-A', 'apical dendrite'), ('skin-A12-22', 'axon'), ('006-red1', 'axon')])

size_df = pd.read_csv('/iblsn/data/Arjun/neurons/neuron_size/neuron_size_synthetic.csv', skipinitialspace=True)

with open('/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner_synthetic2.csv', 'a') as f:
    for neuron_name, neuron_type in missing:
        neuron_type = neuron_type.replace(' ', '_')
        front_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic/%s/%s' % (neuron_name, neuron_type)
        pareto_front = pd.read_csv('%s/pareto_front.csv' % front_dir, skipinitialspace=True)
        tree_costs = pd.read_csv('%s/tree_costs.csv' % front_dir, skipinitialspace=True)
        mcosts, scosts, alphas = pareto_front['mcost'], pareto_front['scost'], pareto_front['alpha']
        tree_costs = tree_costs[tree_costs['tree'] == 'neural']
        neural_mcost, neural_scost = list(tree_costs['mcost'])[0], list(tree_costs['scost'])[0]
        
        neural_dist, neural_index = pareto_dist_scale(mcosts, scosts, neural_mcost, neural_scost)
        neural_alpha = alphas[neural_index]

        neuron_type = neuron_type.replace('_', ' ')
        size = size_df[(size_df['neuron name'] == neuron_name) & (size_df['neuron type'] == neuron_type)]
        points = int(list(size['points'])[0])

        f.write('%s, %s, %d, %f\n' % (neuron_name, neuron_type, points, neural_alpha))
