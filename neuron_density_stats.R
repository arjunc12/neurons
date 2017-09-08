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

colnames1 = c('name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',
             'neural_dist', 'centroid_dist', 'random_dist', 'trials', 
             'successes', 'comparisons', 'dominates')
df1 = read.csv('pareto_mst.csv', col.names=colnames1)
print(df1$cell_type == df1$cell_type[3541])
df = df1[df1$cell_type != df1$cell_type[3541],]

df1$neuron_type = str_sub(df1$name, -1, -1)

colnames2 = c('name', 'points', 'volume', 'mcost')
df2 = read.csv('neuron_density.csv', col.names=colnames2)

df = merge(df1, df2)

model = lm(log10(mcost) ~ log10(volume), data=df)
mcost_resid = model$residuals
df$mcost_resid = mcost_resid
mcost_hat = model$fitted.values
df$mcost_hat = mcost_hat
print(cor(df$alpha, df$mcost_resid))

df$density = df$mcost / df$volume
print(cor(df$alpha, df$volume))

model2 = lm(alpha~cell_type * species * region, data=df)
print(anova(model2))
