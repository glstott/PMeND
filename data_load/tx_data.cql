// Index relationships to speed up downstream queries.
CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);
CREATE INDEX source_relationship2 FOR ()-[r:lica_dist]->() ON (r.source_id);
CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);
CREATE INDEX phylogeny_id FOR (n:phylogeny) ON (n.id);
CREATE INDEX sample_name FOR (n:sample) ON (n.name);

// Install APOC 4.4.0.2 and graph data science 1.8.3 plugins prior to running each script 
// Load data from a CSV file, Texas_cumulative.csv, in the import subdirectory of the database. 
//  CSV with column headers: vendor_accession	vendor	ID	GISAID_accession	GISAID_name	NCBI_accession	NCBI_name	country	state	zip	host	collection_date	receive_date	Date_received_by_CDC	sex	age	race	ethnicity	targeted_sequencing	sequence_machine	assembly_method	specimen_source	original_PANGO_lineage	pango_learn_version	pango_designation_version	pango_class_date	lineage_PANGO_lineage	clade_Nextclade_clade	nextclade_version
//  Note: collection dates are used to generate source numbers for future filtering by tree. 
LOAD CSV WITH HEADERS FROM 'file:///Texas_cumulative.csv' AS line
WITH line, SPLIT(SPLIT(line.zip, "-")[0], ".")[0] AS simple_zip, CASE WHEN toInteger(SPLIT(line.collection_date, "-")[1]) = 8 AND 
                                                                        toInteger(SPLIT(line.collection_date, "-")[2]) <= 14 THEN 0
                    WHEN toInteger(SPLIT(line.collection_date, "-")[1]) = 8 AND toInteger(SPLIT(line.collection_date, "-")[2]) > 14 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) <= 28  THEN 1
                    WHEN (toInteger(SPLIT(line.collection_date, "-")[1]) = 9 AND toInteger(SPLIT(line.collection_date, "-")[2]) > 0 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) <= 11 ) OR (toInteger(SPLIT(line.collection_date, "-")[1]) = 8 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) > 28 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) <= 31) THEN 2
                    WHEN toInteger(SPLIT(line.collection_date, "-")[1]) = 9 AND toInteger(SPLIT(line.collection_date, "-")[2]) > 11 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) <= 25 THEN 3
                    WHEN toInteger(SPLIT(line.collection_date, "-")[1]) = 9 AND toInteger(SPLIT(line.collection_date, "-")[2]) > 25 
                        AND toInteger(SPLIT(line.collection_date, "-")[2]) <= 30 THEN 4 END AS source_num
WHERE source_num IS NOT NULL
// note: we have a manually entered case statement to identify a set of zip codes of note, those in the Houston area.
CREATE (s:sample:phylogeny {name: line.GISAID_name, age: toInteger(line.age), 
    GISAID_accession: line.GISAID_accession, 
    collection_date: line.collection_date, 
    receive_date: line.receive_date, 
    date_received_by_CDC: line.Date_received_by_CDC, clade: line.clade_Nextclade_clade, 
    pango_lineage: line.lineage_PANGO_lineage, sex: line.sex, race: line.race, ethnicity: line.ethnicity,
    city: CASE WHEN simple_zip IN ['77001', '77002', '77003', '77004', '77005', '77006', 
            '77007', '77008', '77009', '77010', '77011', '77012', '77013', '77014', 
            '77015', '77016', '77017', '77018', '77019', '77020', '77021', '77022', 
            '77023', '77024', '77025', '77026', '77027', '77028', '77029', '77030', 
            '77031', '77032', '77033', '77034', '77035', '77036', '77037', '77038', 
            '77039', '77040', '77041', '77042', '77043', '77044', '77045', '77046', 
            '77047', '77048', '77049', '77050', '77051', '77052', '77053', '77054', 
            '77055', '77056', '77057', '77058', '77059', '77060', '77061', '77062', 
            '77063', '77064', '77065', '77066', '77067', '77068', '77069', '77070', 
            '77071', '77072', '77073', '77074', '77075', '77076', '77077', '77078', 
            '77079', '77080', '77081', '77082', '77083', '77084', '77085', '77086', 
            '77087', '77088', '77089', '77090', '77091', '77092', '77093', '77094', 
            '77095', '77096', '77098', '77099', '77201', '77202', '77203', '77204', 
            '77205', '77206', '77207', '77208', '77209', '77210', '77212', '77213', 
            '77215', '77216', '77217', '77218', '77219', '77220', '77221', '77222', 
            '77223', '77224', '77225', '77226', '77227', '77228', '77229', '77230', 
            '77231', '77233', '77234', '77235', '77236', '77237', '77238', '77240', 
            '77241', '77242', '77243', '77244', '77245', '77248', '77249', '77251', 
            '77252', '77253', '77254', '77255', '77256', '77257', '77258', '77259', 
            '77261', '77262', '77263', '77265', '77266', '77267', '77268', '77269', 
            '77270', '77271', '77272', '77273', '77274', '77275', '77277', '77279', 
            '77280', '77282', '77284', '77287', '77288', '77289', '77290', '77291', 
            '77292', '77293', '77297', '77299'] THEN "HOUSTON" ELSE "OTHER" END , 
    simple_zip: simple_zip, zip:line.zip, source_id: source_num, 
    source: CASE WHEN source_num = 0 THEN '0801-0814' 
                WHEN source_num = 1 THEN '0801-0828' 
                WHEN source_num = 2 THEN '0801-0911' 
                WHEN source_num = 3 THEN '0801-0925' 
                WHEN source_num = 4 THEN '0801-0930' END});


// Load phylogenies from neo_ready.csv, which is generated using tree_to_neo4j.ipynb
CALL apoc.periodic.iterate( 
    "LOAD CSV WITH HEADERS FROM 'file:///neo_ready.csv' AS line WITH line WHERE line.type = 'root' OR line.type = 'node' RETURN DISTINCT line.unique_id AS id, line.type AS type", 
    "CREATE (child:LICA:phylogeny {id: id, type: type, source_id: -1})", 
    {batchSize: 1000} );

// Load phylogenies from neo_ready.csv, which is generated using tree_to_neo4j.ipynb
CALL apoc.periodic.iterate(" LOAD CSV WITH HEADERS FROM 'file:///neo_ready.csv' AS line " + "WITH line WHERE line.type = 'leaf' " +
"MATCH (child:sample:phylogeny {name: line.long_id}) " +
"MATCH (parent:LICA:phylogeny {id: line.parent_unique_id}) " +
"RETURN child, parent, line",
"CREATE (child)-[:child_of {source: line.source, distance: toFloat(line.distance), "+
    'source_id: CASE WHEN line.source = "GISAID_combined.0801-0814.aligned.nwk" THEN 0 '+
                    'WHEN line.source = "GISAID_combined.0801-0828.aligned.nwk" THEN 1 '+
                    'WHEN line.source = "GISAID_combined.0801-0911.aligned.nwk" THEN 2 '+
                    'WHEN line.source = "GISAID_combined.0801-0925.aligned.nwk" THEN 3 '+
                    'WHEN line.source = "GISAID_combined.0801-0930.aligned.nwk" THEN 4 END}]->(parent)',
                    {batchSize: 1000})
;

// Load phylogenies from neo_ready.csv, which is generated using tree_to_neo4j.ipynb
CALL apoc.periodic.iterate("LOAD CSV WITH HEADERS FROM 'file:///neo_ready.csv' AS line " +
"WITH line " +
'WHERE line.type <> "leaf" AND line.type <> "root" ' +
'MATCH (child:LICA:phylogeny {id: line.unique_id}) ' +
"MATCH (parent:LICA:phylogeny {id: line.parent_unique_id}) " +
"RETURN child, parent, line",
'CREATE (child)-[:child_of {source: line.source, distance: toFloat(line.distance), ' +
    'source_id: CASE WHEN line.source = "GISAID_combined.0801-0814.aligned.nwk" THEN 0 ' +
                   ' WHEN line.source = "GISAID_combined.0801-0828.aligned.nwk" THEN 1 ' +
                   ' WHEN line.source = "GISAID_combined.0801-0911.aligned.nwk" THEN 2 ' +
                   ' WHEN line.source = "GISAID_combined.0801-0925.aligned.nwk" THEN 3 ' +
                   ' WHEN line.source = "GISAID_combined.0801-0930.aligned.nwk" THEN 4 END}]->(parent) ', 
                   {batchSize: 1000})
;