#include <nauty.h>
#include <nautinv.h>

optionblk defaultoptions_graph() {
	DEFAULTOPTIONS_GRAPH(options);
	return options;
}

optionblk defaultoptions_digraph() {
	DEFAULTOPTIONS_DIGRAPH(options);
	return options;
}


long wordsize() {
	return WORDSIZE;
}

void baked_options(graph *g, int *canonical_labelling, int *partition, int *orbits,
		   statsblk *stats, int num_setwords, int num_vertices, graph *canonical_graph) {
	/* statsblk stats; */
	DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = 1;
	options.digraph = 1;

	return densenauty(g, canonical_labelling, partition, orbits, &options, stats, num_setwords, num_vertices, canonical_graph);
}

void baked_options_color(graph *g, int *canonical_labelling, int *partition, int *orbits,
		   statsblk *stats, int num_setwords, int num_vertices, graph *canonical_graph) {
	/* statsblk stats; */
	DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = 1;
	options.digraph = 1;
	options.defaultptn = 0;

	return densenauty(g, canonical_labelling, partition, orbits, &options, stats, num_setwords, num_vertices, canonical_graph);
}

void baked_options_and_stats(graph *g, int *canonical_labelling, int *partition, int *orbits,
		   int num_setwords, int num_vertices, graph *canonical_graph) {
	statsblk stats;
	DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = 1;
	options.digraph = 1;

	return densenauty(g, canonical_labelling, partition, orbits, &options, &stats, num_setwords, num_vertices, canonical_graph);
}

/*
void selftest() {
	graph g[4] = {0, 0, 0, 0};
	graph outg[4];
	int labelling[4];
	int partition[4];
	int orbits[4];
	canonical_form(g, 1, 4, labelling, partition, orbits, outg);
}

// Proof of concept C function that receives and performs some work on a graph.
setword graph_receiver(graph *g, int len) {
	setword acc = 0;
	for (int i=0; i<len; i++) {
		acc += g[i] >> 32;
	}
	return acc;
}

 nauty(g->matrix, g->lab, g->ptn, NULL, g->orbits,
         g->options, g->stats,  g->workspace, g->worksize,
         g->no_setwords, g->no_vertices, NULL); */
