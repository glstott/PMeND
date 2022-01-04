// Create integer field for subgraph generation. 
//    Subgraphs are only available with numeric values on the relationships extracted.
MATCH (n:sample)-[r:calculated_distance {source: "GISAID_combined.0801-0814.aligned.nwk"}]->(m:sample) 
SET r.source_id = 0

// Calculate distance from each LICA. 
//    If we measure all distances to LICA, then sum incoming pairs, we will have patristic distances 
//    without necessitating slower algorithms. Batched to keep it from increasing memory consumption.
CALL apoc.periodic.iterate('MATCH (l:LICA)<-[:child_of {source: "GISAID_combined.0801-0925.aligned.nwk"}]-(:phylogeny) ' + 
'OPTIONAL MATCH (l)-[:child_of {source: "GISAID_combined.0801-0925.aligned.nwk"}]->(l2:LICA) ' + 
'WITH l, l2 ' + 
'WHERE l2 is null ' + 
'WITH l ' + 
'MATCH p=(s1:sample)-[:child_of * {source: "GISAID_combined.0801-0925.aligned.nwk"}]->(l) RETURN p, s1 ' ,
'FOREACH (t in nodes(p) | CREATE (s1)-[:lica_dist {source: "0801-0925", source_id: 3, distance: ' + 
'reduce(total=0.0, r IN relationships(p) | r.distance+total)}]->(t))', {batchSize:100});

// Sum distances to generate the patristic distance network.
CALL apoc.periodic.iterate('MATCH (s1:sample)-[r1:lica_dist]->(:LICA)<-[r2:lica_dist]-(s2:sample)' +
'WHERE id(s1) < id(s2) AND r1.source_id = r2.source_id' +
'RETURN s1, s2, r1, r2',
'CREATE (s1)-[:calculated_distance {source: r1.source, source_id: r1.source_id, distance: r1.distance + r2.distance}]->(s2)',
{batchSize: 100});

// Create complete graph before creating subgraphs for each data subset. 
CALL gds.graph.drop(["distance_graph",'0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930']);
CALL gds.graph.create(
  'distance_graph',                    
  {sample: {properties: ['source_id']}},                             
  {calculated_distance: {orientation: 'UNDIRECTED', properties: ['source_id', 'distance']}}  
)
YIELD
  graphName AS graph,
  relationshipProjection AS knowsProjection,
  nodeCount AS nodes,
  relationshipCount AS rels
  ;

// Create Subgraph for analysis. Only the relevant nodes and relationships for each date range are selected.
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
CALL gds.beta.graph.create.subgraph(dt_range,
  "distance_graph", 'n.source_id <= ' + CASE WHEN dt_range = '0801-0814' THEN '0'
        WHEN dt_range = '0801-0828' THEN '1' 
        WHEN dt_range = '0801-0911' THEN '2'
        WHEN dt_range = '0801-0925' THEN '3'
        WHEN dt_range = '0801-0930' THEN '4' END,
  'r.source_id = ' + CASE WHEN dt_range = '0801-0814' THEN '0'
        WHEN dt_range = '0801-0828' THEN '1' 
        WHEN dt_range = '0801-0911' THEN '2'
        WHEN dt_range = '0801-0925' THEN '3'
        WHEN dt_range = '0801-0930' THEN '4' END
) YIELD graphName, fromGraphName, nodeFilter, relationshipFilter, 
    nodeCount, relationshipCount, createMillis
RETURN dt_range, graphName, fromGraphName, nodeFilter, relationshipFilter, 
    nodeCount, relationshipCount, createMillis
;

// Generate forest of Minimum Spanning Trees, starting with each of the day 0 subjects as an origin. 
//    This should decrease the bias 
MATCH (n:sample {collection_date: '2021-08-01'})
WITH n
WITH collect(id(n)) AS ids
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
UNWIND ids AS iteration_number
CALL gds.alpha.spanningTree.minimum.write(dt_range, {
  startNodeId: iteration_number,
  relationshipWeightProperty: 'distance',
  writeProperty: 'MST_' + iteration_number + '_' + dt_range,
  weightWriteProperty: 'distance_combined'
})
YIELD createMillis, computeMillis, writeMillis, effectiveNodeCount
RETURN dt_range, iteration_number, createMillis, computeMillis, writeMillis, effectiveNodeCount
;

// Compute joint Minimum Spanning Tree
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
MATCH (n:sample)-[r]->(n2:sample)
WHERE type(r) ENDS WITH dt_range
WITH n, n2, COUNT(DISTINCT r) AS edge_width, AVG(r.distance_combined) AS mean_dist, dt_range
CREATE (n)-[:joint_MST {source:dt_range, edge_width:edge_width, mean_dist: mean_dist}]->(n2)
;


/////////////////// Sketch of Clustering script. Needs work. ////////////////////////////

call apoc.custom.asFunction('getLastCluster',
"MATCH (n) where id(n) = id($n)
WITH n 
MATCH clustPath = (n)-[:IN_CLUSTER*0..]->(c)
WITH clustPath order by length(clustPath) desc
LIMIT 1
RETURN nodes(clustPath)[-1] as lastCluster",
"NODE",
[["n","NODE"]],
false,
"Starting from a given node, get the node at the end of the longest path of IN_CLUSTER relationships")
;

MATCH (a)-[m:joint_MST]->(b)
WHERE m.evaluated is null
WITH a, b, m
ORDER BY m.distance
LIMIT 1
WITH custom.getLastCluster(a).lastCluster as c1, custom.getLastCluster(b).lastCluster as c2, m
CREATE (c1)-[r1:IN_CLUSTER]->(c:Cluster)<-[r2:IN_CLUSTER]-(c2)
SET r1.splitDistance = m.distance,
r2.splitDistance = m.distance,
c.size = coalesce(c1.size, 1) + coalesce(c2.size, 1),
m.evaluated = True
RETURN *
;

call apoc.periodic.commit(
"MATCH (a)-[m:joint_MST]->(b)
WHERE m.evaluated is null
WITH a, b, m
ORDER BY m.distance
LIMIT 1 WITH 
custom.getLastCluster(a).lastCluster as cluster1, custom.getLastCluster(b).lastCluster as cluster2, m
CREATE (cluster1)-[r1:IN_CLUSTER]->(c:Cluster)
                <-[r2:IN_CLUSTER]- (cluster2)
SET r1.splitDistance = m.distance,
r2.splitDistance = m.distance,
c.size = coalesce(cluster1.size, 1) + coalesce(cluster2.size, 1),
m.evaluated = True
RETURN count(*)", {})
;