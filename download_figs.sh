echo "figs"
scp -r achandrasekhar@doritos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs ./
echo "synthetic figs"
scp -r achandrasekhar@doritos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs_synthetic ./
echo "plots"
scp -r achandrasekhar@cheetos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_front_plots0.2 ./
echo "log plots"
scp -r achandrasekhar@cheetos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_front_log_plots0.2 ./
echo "synthetic plots"
scp -r achandrasekhar@cheetos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_front_plots_synthetic0.2 ./
echo "synthetic log plots"
scp -r achandrasekhar@cheetos.snl.salk.edu:/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_front_log_plots_synthetic0.2 ./
echo "drawings"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/drawings ./
echo "stats"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/steiner_stats ./
echo "synthetic stats"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/steiner_stats_synthetic ./
echo "imaris"
scp -r achandrasekhar@doritos.snl.salk.edu:/iblsn/data/Arjun/neurons/imaris/figs/. steiner_imaris
echo "neuron builder"
scp achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/neuron_builder/*.pdf neuron_builder/
echo "test runtimes"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/test_runtimes ./
echo "greedy eval"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/greedy_eval/*cost.pdf greedy_eval
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/greedy_eval/*hist.pdf greedy_eval
echo "neuron density"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/neuron_density ./
echo "cost bounds"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/cost_bounds ./
echo "truncation"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/truncation ./
echo "boutons"
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/boutons/drawings/ ./boutons/
scp -r achandrasekhar@doritos.snl.salk.edu:/home/achandrasekhar/neurons/boutons/histograms/ ./boutons/
echo "compare algorithms"
scp -r achandrasekhar@cheetos.snl.salk.edu:/home/achandrasekhar/neurons/compare_algorithms/compare_algorithms*.pdf compare_algorithms/
echo "triplets"
scp achandrasekhar@cheetos.snl.salk.edu:/home/achandrasekhar/neurons/triplet*.csv triplets/
echo "branch angles"
scp -r achandrasekhar@cheetos.snl.salk.edu:/home/achandrasekhar/neurons/branch_angles ./
echo "sandbox"
scp -r achandrasekhar@cheetos.snl.salk.edu:/home/achandrasekhar/neurons/sandbox ./
