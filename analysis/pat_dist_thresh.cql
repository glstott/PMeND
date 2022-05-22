// Loop through thresholds to plot average size and number of clusters
UNWIND range(4,180,1) AS threshold
CALL gds.graph.drop('calc_dist') YIELD graphName as graph1
WITH graph1, threshold
CALL gds.graph.project.cypher(
  'calc_dist',
  'MATCH (n:sample) RETURN id(n) AS id',
  'MATCH (n:sample)-[r:calculated_distance]->(m:sample) WHERE r.distance <= ' + threshold / 365.0 + ' RETURN id(n) AS source, id(m) AS target')
YIELD
  graphName AS graph, nodeQuery, nodeCount AS nodes, relationshipQuery, relationshipCount AS rels
WITH graph, threshold
CALL gds.wcc.stream(graph)
YIELD nodeId, componentId
WITH threshold, componentId, COUNT(DISTINCT nodeId) AS clusterSize
WHERE clusterSize >= 2
RETURN threshold, COUNT(DISTINCT componentId) AS clusterCount, AVG(clusterSize) AS avg_ClusterSize, SUM(clusterSize) AS totalRepresented, 
  MAX(clusterSize) AS max_cluster;
