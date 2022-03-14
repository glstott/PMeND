// Calculate distance from each LICA. 
//    If we measure all distances to LICA, then sum incoming pairs, we will have patristic distances 
//    without necessitating slower algorithms. 
UNWIND [0] AS idval
MATCH p=(s1:sample)-[:child_of * {source_id: idval}]->(:LICA {type: "root"}) 
WHERE s1.source_id <= idval 
WITH s1, idval, p
UNWIND range(1,size(nodes(p))-1) AS t_index 
WITH nodes(p)[t_index] AS target, s1, reduce(total=0.0, r IN relationships(p)[..t_index] | r.distance+total) AS calc_dist, idval
WHERE calc_dist <= .5  
WITH target,s1, calc_dist, idval 
CREATE (s1)-[:lica_dist {source_id: idval, distance: calc_dist}]->(target);

// Sum distances to generate the patristic distance network.
CALL apoc.periodic.iterate('MATCH (l:LICA) ' +
'WITH l ' +
'MATCH (s1:sample)-[r1:lica_dist]->(l:LICA)<-[r2:lica_dist]-(s2:sample) ' +
'WHERE id(s1) < id(s2) AND r1.source_id = r2.source_id  ' +
'RETURN s1, s2, r1.source_id AS srcid, r1.source AS src, r1.distance + r2.distance AS dist',
'CREATE (s1)-[:calculated_distance {source_id: srcid, distance: dist}]->(s2) ',
{batchSize: 1000});
CALL apoc.periodic.iterate("MATCH (:phylogeny)-[r:lica_dist]->(:phylogeny) RETURN r", "DELETE r", {batchSize: 1000});

// Create complete graph before creating subgraphs for each data subset. 
//CALL gds.graph.drop(["distance_graph",'new']);
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


// Generate forest of Minimum Spanning Trees, starting with each of the day 0 subjects as an origin. 
//    This should decrease the bias 
UNWIND ['delta'] AS dt_range
// CALL {
//   WITH dt_range
//   CALL gds.graph.drop('new', false) YIELD graphName AS gn
//   WITH dt_range, gn
//   CALL gds.beta.graph.create.subgraph("new",
//     "distance_graph", 'n.source_id <= 0',
//     'r.source_id = 0'
//   ) YIELD graphName, fromGraphName, nodeFilter, relationshipFilter, 
//       nodeCount, relationshipCount, createMillis
//       RETURN graphName, fromGraphName, nodeFilter, relationshipFilter, 
//       nodeCount, relationshipCount, createMillis
// }
CALL {
  WITH dt_range
  MATCH (n:sample)
  WITH dt_range, n LIMIT 2000
  WITH dt_range, collect(id(n)) AS ids
  UNWIND ids AS iteration_number
  CALL gds.alpha.spanningTree.minimum.write("new", {
    startNodeId: iteration_number,
    relationshipWeightProperty: 'distance',
    writeProperty: 'MST_' + iteration_number + '_' + dt_range,
    weightWriteProperty: 'distance_combined'
  })
  YIELD createMillis, computeMillis, writeMillis, effectiveNodeCount
  RETURN createMillis AS cM, computeMillis AS coM, writeMillis AS wM, effectiveNodeCount AS enc, iteration_number
}
RETURN dt_range, iteration_number, cM, coM, wM, enc
;

// Compute joint Minimum Spanning Tree
UNWIND ['delta'] AS dt_range
CALL apoc.periodic.iterate(

    'MATCH (:sample)-[r]->(:sample) '+
    'WHERE type(r) ENDS WITH "' + dt_range + '" '+
    'WITH COUNT(DISTINCT type(r)) AS max_width '+
    'MATCH (n:sample)-[r]->(n2:sample) '+
    'WHERE type(r) ENDS WITH "' + dt_range + '"  AND id(n) < id(n2)'+
    'RETURN max_width, n, n2, COUNT(DISTINCT r) AS edge_width, AVG(r.distance_combined) AS mean_dist',

    'CREATE (n)-[:joint_MST {source:"' + dt_range + '", edge_width:edge_width / toFloat(max_width), mean_dist: mean_dist}]->(n2)', {batchSize: 10000}

) YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;


// Delete each single MST
UNWIND ['0801-0814', '0801-0828', '0801-0911', '0801-0925', '0801-0930'] AS dt_range
CALL apoc.periodic.iterate("MATCH (:sample)-[r]->(:sample)
WHERE type(r) ENDS WITH '"+ dt_range +"' RETURN r", "DELETE r", {batchSize: 10000}) 
YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;