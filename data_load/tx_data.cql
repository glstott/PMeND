LOAD CSV WITH HEADERS FROM 'file:///Texas_cumulative.csv' AS line
CREATE (s:sample {name: line.GISAID_name, age: toInteger(line.age), GISAID_accession: line.GISAID_accession, collection_date: line.collection_date, receive_date: line.receive_date, date_received_by_CDC: line.Date_received_by_CDC})
MATCH (sex:characteristic:sex {id:line.sex})
MATCH (race:characteristic:race {id:line.race})
MATCH (ethnicity:characteristic:ethnicity {id:line.ethnicity})
CREATE (s)-[:has_characteristic]->(sex)
CREATE (s)-[:has_characteristic]->(race)
CREATE (s)-[:has_characteristic]->(ethnicity)