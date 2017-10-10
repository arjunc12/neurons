library(ICC)
library(stringr)

getICCs <- function(df)
{
    for (cat in c("cell_type", "species", "region"))
    {
        df2 = unique(df, by=c("name", cat))
        df2$count = ave(df2$alpha, df2[,cat], FUN=length)
        df2 = df2[df2$count > 25,]
        print(cat)
        print(ICCbare(x=cat, y=alpha, data=df2))
    }
}

colnames = c('name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',
             'neural_dist', 'centroid_dist', 'random_dist', 'trials', 
             'successes', 'comparisons', 'dominates')
df = read.csv('pareto_mst.csv', col.names=colnames)
#print(df$cell_type == df$cell_type[3541])
#df = df[df$cell_type != df$cell_type[3541],]

df$neuron_type = str_sub(df$name, -1, -1)

print("-----all-----")
getICCs(df)

print("-----axons-----")
getICCs(df[df$neuron_type == '0',])

print("-----basal dendrite-----")
getICCs(df[df$neuron_type == '1',])

print("-----apical dendrite-----")
getICCs(df[df$neuron_type == '2',])

print("-----dendrite-----")
getICCs(df[df$neuron_type != '0',])

print("-----neuron type-----")
print(ICCbare(x=neuron_type, y=alpha, data=df))

print("-----axon vs dendrite-----")
df$axon = as.numeric(df$neuron_type == '0')
print(ICCbare(x=axon, y=alpha, data=df))
