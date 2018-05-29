library(ICC)
library(stringr)

OUTPUT_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner.csv'
SYNTHETIC_OUTPUT_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner_synthetic.csv'
CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES_FILE_FILTERED = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories_filtered.csv'
CATEGORIES = c('cell.type', 'species', 'region', 'neuron.type', 'lab')
RATIO_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/tradeoff_ratio_synthetic.csv'

MIN_COUNT = 50
MIN_POINTS = 50

getICCs <- function(df, yval)
{
    for (cat in CATEGORIES)
    {
        df2 = unique(df, by=c("neuron.name", "neuron.type", cat))
        df2$count = ave(df2$alpha, df2[,cat], FUN=length)
        df2 = df2[df2$count > MIN_COUNT,]
        print(cat)
        print(ICCbare(x=cat, y=yval, data=df2))
    }
}

alphaICCs <- function(df)
{
    print("alpha")
    getICCs(df, yval='alpha')
}

tradeoffICCs <- function(df)
{
    print("tradeoff ratio")
    getICCs(df, yval='tradeoff.ratio')
}
#alphaICCs <- function(df)
#{
#    for (cat in CATEGORIES)
#    {
#        df2 = unique(df, by=c("neuron.name", "neuron.type", cat))
#        df2$count = ave(df2$alpha, df2[,cat], FUN=length)
#        df2 = df2[df2$count > MIN_COUNT,]
#        print(cat)
#        print(ICCbare(x=cat, y=alpha, data=df2))
#    }
#}

categories_df = read.csv(CATEGORIES_FILE_FILTERED, strip.white=TRUE)
output_df = read.csv(SYNTHETIC_OUTPUT_FILE, strip.white=TRUE)
ratio_df = read.csv(RATIO_FILE, strip.white=TRUE)
print(colnames(categories_df))
print(colnames(output_df))
print(colnames(ratio_df))
df = merge(x=categories_df, y=output_df, by="neuron.name")
df = merge(x=df, y=ratio_df, by=c("neuron.name", "neuron.type"))
print(colnames(df))
df = df[df$points >= MIN_POINTS,]

print("-----all-----")
alphaICCs(df)
tradeoffICCs(df)

print("-----axons-----")
alphaICCs(df[df$neuron.type == 'axon',])
tradeoffICCs(df[df$neuron.type == 'axon',])

print("-----basal dendrite-----")
alphaICCs(df[df$neuron.type == 'basal dendrite',])
tradeoffICCs(df[df$neuron.type == 'basal dendrite',])

print("-----apical dendrite-----")
alphaICCs(df[df$neuron.type == 'apical dendrite',])
tradeoffICCs(df[df$neuron.type == 'apical dendrite',])

print("-----dendrite-----")
alphaICCs(df[df$neuron.type != 'axon' & df$neuron.type != 'truncated axon',])
tradeoffICCs(df[df$neuron.type != 'axon' & df$neuron.type != 'truncated axon',])

print("-----neuron type-----")
print(ICCbare(x=neuron.type, y=alpha, data=df))
print(ICCbare(x=neuron.type, y="tradeoff.ratio", data=df))

print("-----axon vs dendrite-----")
df$axon = as.numeric(df$neuron.type == "axon" | df$neuron.type == 'truncated axon')
print(ICCbare(x=axon, y=alpha, data=df))
print(ICCbare(x=axon, y="tradeoff.ratio", data=df))
