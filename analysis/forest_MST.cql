// Create integer field for subgraph generation
MATCH (n:sample)-[r:calculated_distance {source: "GISAID_combined.0801-0814.aligned.nwk"}]->(m:sample) 
SET r.source_id = 0

// Create graph before paring it down
CALL gds.graph.create(
  'distance_graph',                    
  'sample',                             
  {calculated_distance: {orientation: 'UNDIRECTED', properties: ['source_id', 'distance']}}  
)
YIELD
  graphName AS graph,
  relationshipProjection AS knowsProjection,
  nodeCount AS nodes,
  relationshipCount AS rels

// Create Subgraph for analysis
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
CALL gds.beta.graph.create.subgraph(dt_range,
  "distance_graph", '*',
  'r.source_id = ' + CASE WHEN dt_range = '0801-0814' THEN '0'
        WHEN dt_range = '0801-0828' THEN '1' 
        WHEN dt_range = '0801-0911' THEN '2'
        WHEN dt_range = '0801-0925' THEN '3'
        WHEN dt_range = '0801-0930' THEN '4' END
) YIELD graphName, fromGraphName, nodeFilter, relationshipFilter, 
    nodeCount, relationshipCount, createMillis
RETURN dt_range, graphName, fromGraphName, nodeFilter, relationshipFilter, 
    nodeCount, relationshipCount, createMillis


// Generate forest of Minimum Spanning Trees
MATCH (n:sample {name: 'hCoV-19/USA/TX-CDC-ASC210186970/2021'})
WITH n
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
UNWIND RANGE(0,100,1) AS iteration_number
CALL gds.alpha.spanningTree.minimum.write(dt_range, {
  startNodeId: id(n),
  relationshipWeightProperty: 'distance',
  writeProperty: 'MST_' + iteration_number + '_' + dt_range,
  weightWriteProperty: 'distance_combined'
})
YIELD createMillis, computeMillis, writeMillis, effectiveNodeCount
RETURN dt_range, iteration_number, createMillis, computeMillis, writeMillis, effectiveNodeCount

// Compute joint Minimum Spanning Tree
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
MATCH (n:sample)-[r]->(n2:sample)
WHERE type(r) ENDS WITH dt_range
RETURN n, n2, COUNT(DISTINCT r) AS edge_width, dt_range
