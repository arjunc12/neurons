library(ICC)
library(stringr)

getICCs <- function(df)
{
    print('cell type')
    print(ICCbare(x=cell_type, y=alpha, data=df))
    print('species')
    print(ICCbare(x=species, y=alpha, data=df))
    print('region')
    print(ICCbare(x=region, y=alpha, data=df))  
}

colnames = c('name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',
             'neural_dist', 'centroid_dist', 'random_dist', 'trials', 
             'successes', 'comparisons', 'dominates')
df = read.csv('pareto_mst.csv', col.names=colnames)

df$neuron_type = str_sub(df$name, -1, -1)

print("axons")
getICCs(df[df$neuron_type == '0',])

print("basal dendrite")
getICCs(df[df$neuron_type == '1',])

print("apical dendrite")
getICCs(df[df$neuron_type == '2',])
