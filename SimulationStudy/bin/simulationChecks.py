import pyslim, sys
import msprime, tskit
import matplotlib.pyplot as pyplot

## This script uses code and functions from the PySlim and TsKit tutorials to perform sanity checks and exploratory analyses on tree sequences that come from SLiM

def remove_singletons(ts):
	tables = ts.dump_tables()
	tables.sites.clear()
	tables.mutations.clear()
	for tree in ts.trees():
		for site in tree.sites():
			assert len(site.mutations) == 1
			mut = site.mutations[0]
			if tree.num_samples(mut.node) > 1:
				site_id = tables.sites.add_row(
					position = site.position, ancestral_state = site.ancestral_state
					)
				tables.mutations.add_row(
					site=site_id, node=mut.node, derived_state=mut.derived_state
					)
	return tables.tree_sequence()

def ld_matrix(ts, figure_name):
	ld_calc = tskit.LdCalculator(ts)
	A = ld_calc.r2_matrix()
	x = A.shape[0] / pyplot.rcParams["figure.dpi"]
	x = max(x, pyplot.rcParams["figure.figsize"][0])
	fig, ax = pyplot.subplots( figsize = (x,x) )
	fig.tight_layout(pad = 0)
	im = ax.imshow(A, interpolation = "none", vmin = 0, vmax =1, cmap ="Blues")
	ax.set_xticks([])
	ax.set_yticks([])
	for s in "top","bottom","left","right":
		ax.spines[s].set_visible(False)
	pyplot.gcf().colorbar(im, shrink = 0.5, pad = 0)
	pyplot.savefig(figure_name)
	
	
	
orig_ts = pyslim.load(sys.argv[1])


orig_max_roots = max(t.num_roots for t  in orig_ts.trees())

print("number of roots =",orig_max_roots)
if orig_max_roots != 1:
	print("\t\tTree did not Coalesce!!")

print("number of samples (individuals in the SLiM simulation) =", len(orig_ts.individuals_alive_at(0)))

variant_positions = [v.position for v in orig_ts.variants()]

print("number of funtional variants in the population =", len(variant_positions))

for v in orig_ts.variants():
	print( v.genotypes.sum() / len(v.genotypes) , v.position)

## Now let's mutate the tree and look at LD.. this will take a fair amount of time and memory...
mut_ts = msprime.mutate( orig_ts, rate = 1e-10, keep = True) 
all_variant_positions = [v.position for v in mut_ts.variants()]

print("number of functional + neutral variants in the population =", len(all_variant_positions))

ld_matrix(orig_ts, sys.argv[1]+'.LD.svg')



