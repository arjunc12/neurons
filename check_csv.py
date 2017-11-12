FNAME = 'pareto_steiner.csv'
COLUMNS = ['name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',\
           'neural_dist', 'centroid_dist', 'random_dist', 'trials', 'successes']
NCOLS = len(COLUMNS)
print NCOLS

f = open(FNAME)
min_cols = float("inf")
max_cols = 0
for line in f:
    line = line.strip('\n')
    line = line.split(',')
    if len(line) != NCOLS:
        print line
    assert len(line) == NCOLS
    min_cols = min(min_cols, len(line))
    max_cols = max(max_cols, len(line))
f.close()
print min_cols
print max_cols
