# Phylogeny and Metadata Network Database (PMeND)

A graph database which enables storage of phylogenies as Tree Aligned Graphs (TAGs) and integrates these data with sample metadata. The database is then used to generate  Work in progress, there will be a refactor over the next few weeks to speed up a couple of steps in the pipeline, operationalize code meant to run in a particular file structure, and add a NextFlow workflow.

## Current Structure

`data_load`: contains a quick python script to extract an edge list for NWK trees in a directory and cypher statements to set up indexes, and add data to the Neo4j database. 

`phylogeny_generation`: contains bash scripts to generate an alignment on sequences and build ML trees for each.

`analysis`: Includes a subfolder for cytoscape styles used, a quick python script to generate lat/long coordinates for each zipcode, and a series of cypher statements used to generate a patristic distance network; build a forest of minimum spanning trees; and summarize the forest of minimum spanning trees.

`nextflow`: Contains a frame for generating a nextflow workflow to populate and maintain the database.

## Future Work

Issues show the current workstreams in progress. Patristic distance calculation updates and nextflow workflow updates will be coming up first.
