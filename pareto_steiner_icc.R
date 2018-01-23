library(ICC)
library(stringr)

MIN_COUNT = 10
MIN_POINTS = 50

getICCs <- function(df)
{
    for (cat in c("cell_type", "species", "region"))
    {
        df2 = unique(df, by=c("name", cat))
        df2$count = ave(df2$alpha, df2[,cat], FUN=length)
        df2 = df2[df2$count > MIN_COUNT,]
        print(cat)
        print(ICCbare(x=cat, y=alpha, data=df2))
    }
}


#COLUMNS = ['name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',\
#           'norm_alpha', 'neural_dist', 'centroid_dist', 'random_dist',\
#           'norm_neural_dist', 'norm_centroid_dist', 'norm_random_dist',\
#           'trials', 'successes', 'norm_successes']

#cnames = c('name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',
#             'neural_dist', 'centroid_dist', 'random_dist', 'trials', 
#             'successes')

cnames = c('name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',
           'norm_alpha', 'neural_dist', 'centroid_dist', 'random_dist',
           'norm_neural_dist', 'norm_centroid_dist', 'norm_random_dist',
           'trials', 'successes', 'norm_successes')
df = read.csv('pareto_steiner.csv', col.names=cnames, header=FALSE)
df = unique(df, by=c("name", "cell_type", "species", "region"))
#print(df$cell_type == df$cell_type[3541])
#df = df[df$cell_type != df$cell_type[3541],]
df = df[df$points >= MIN_POINTS,]

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
