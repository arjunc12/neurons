library(ICC)
library(stringr)

OUTPUT_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner.csv'
SYNTHETIC_OUTPUT_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner_synthetic.csv'
CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES = c('cell.type', 'species', 'region', 'neuron.type', 'lab')

MIN_COUNT = 50
MIN_POINTS = 50

getICCs <- function(df)
{
    for (cat in CATEGORIES)
    {
        df2 = unique(df, by=c("neuron.name", "neuron.type", cat))
        df2$count = ave(df2$alpha, df2[,cat], FUN=length)
        df2 = df2[df2$count > MIN_COUNT,]
        print(cat)
        print(ICCbare(x=cat, y=alpha, data=df2))
    }
}

categories_df = read.csv(CATEGORIES_FILE, strip.white=TRUE)
output_df = read.csv(SYNTHETIC_OUTPUT_FILE, strip.white=TRUE)
df = merge(x=categories_df, y=output_df, by="neuron.name")
df = df[df$points >= MIN_POINTS,]

print("-----all-----")
getICCs(df)

print("-----axons-----")
getICCs(df[df$neuron.type == 'axon',])

print("-----basal dendrite-----")
getICCs(df[df$neuron.type == 'basal dendrite',])

print("-----apical dendrite-----")
getICCs(df[df$neuron.type == 'apical dendrite',])

print("-----dendrite-----")
getICCs(df[df$neuron.type != 'axon' & df$neuron.type != 'truncated axon',])

print("-----neuron type-----")
print(ICCbare(x=neuron.type, y=alpha, data=df))

print("-----axon vs dendrite-----")
df$axon = as.numeric(df$neuron.type == "axon" | df$neuron.type == 'truncated axon')
print(ICCbare(x=axon, y=alpha, data=df))
