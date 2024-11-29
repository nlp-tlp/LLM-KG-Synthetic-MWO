import os
import csv
import json
from dotenv import load_dotenv
from neo4j import GraphDatabase

# Function to create label names
def create_label_name(tokens, start, end):
    """ Create label name from tokens. """
    label = " ".join(tokens[start:end])
    label = label.lower().strip()
    return label

# Function to create nodes
def create_nodes(tx, entities, unique_entities, tokens, entry_id):
    """ Create nodes for entities. """
    current_entities = []
    for entity in entities:
        entity_text = create_label_name(tokens, entity["start"], entity["end"])
        entity_type = entity["type"].split('/')[0]
        sub_types = entity["type"].split('/')[1:]
        unique_entity_key = (entity_text, entity_type)
        current_entities.append(unique_entity_key)

        if unique_entity_key not in unique_entities:
            entity_id = len(unique_entities)
            unique_entities[unique_entity_key] = entity_id
            properties = {"id": entity_id, "text": entity_text, "type": entity_type, "entry_id": [entry_id]}
            query = f"MERGE (n:{entity_type} {{id: $id, text: $text, type: $type, entry_id: $entry_id}}) "
            for i, subtype in enumerate(sub_types):
                properties[f"subtype{i}"] = subtype
                query += f"SET n.subtype{i} = $subtype{i} "
            tx.run(query, properties)
        else:
            entity_id = unique_entities[unique_entity_key]
            tx.run(
                "MATCH (n) WHERE n.id = $id "
                "SET n.entry_id = n.entry_id + $entry_id",
                id=entity_id, entry_id=[entry_id]
            )
    return current_entities, unique_entities

# Function to create relations
def create_relations(tx, relations, unique_entities, current_entities):
    """ Create relations between entities. """
    for relation in relations:
        head = relation["head"]
        tail = relation["tail"]
        rel_type = relation["type"].replace("/", "_")
        head_id = unique_entities[current_entities[head]]
        tail_id = unique_entities[current_entities[tail]]
        tx.run(
            f"MATCH (a {{id: $head_id}}), (b {{id: $tail_id}}) "
            "WHERE NOT a:Entry AND NOT b:Entry "
            f"MERGE (a)-[:{rel_type}]->(b)", 
            head_id=head_id, tail_id=tail_id, type=rel_type
        )

# Function to create entry nodes and connect to its entities
def create_entry(tx, entry_text, entry_id):
    """ Create nodes for entries and connect to its entities. """
    # Create entry node
    tx.run(
        "MERGE (n:Entry {id: $id, text: $text})",
        id=entry_id, text=entry_text
    )
    # Connect entry node to its entities
    tx.run(
        "MATCH (a:Entry {id: $id}), (b) "
        "WHERE ANY(eid IN b.entry_id WHERE eid = $entry_id) "
        "MERGE (b)-[:comesFrom]->(a)",
        id=entry_id, entry_id=entry_id
    )

# Function to read MaintIE manual mapping of failure mode codes
def read_failure_mode_mapping(filepath):
    with open(filepath, 'r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)
        failure_mode_mapping = {}
        for row in reader:
            if row[3] != "":
                failure_mode_mapping[row[0]] = row[3]
    return failure_mode_mapping # entry_text: failure_mode

# Function to label entry to failure mode codes
def entry_failure_mode(tx, failure_mode_mapping):
    for entry_text, failure_mode in failure_mode_mapping.items():
        tx.run(
            "MATCH (a:Entry {text: $entry_text}) "
            "SET a.failure_mode = $failure_mode",
            entry_text=entry_text, failure_mode=failure_mode
        )

# Function to create graph
def create_graph(tx, data):
    """ Create graph from MaintIE dataset. """
    unique_entities = {}
    for entry_id, entry in enumerate(data):
        entities = entry["entities"]
        relations = entry["relations"]
        tokens = entry["tokens"]
        current_entities, unique_entities = create_nodes(tx, entities, unique_entities, tokens, entry_id)
        create_relations(tx, relations, unique_entities, current_entities)
        create_entry(tx, entry['text'], entry_id) 
        entry_failure_mode(tx, read_failure_mode_mapping('../data/MaintIE/gold_undesirable_mapped.csv'))

    print(f"Created {len(unique_entities)} entities.")
    print(f"Created {len(data)} entry entities.")
    print(f"Total number of nodes: {len(unique_entities) + len(data)}")

if __name__ == "__main__":
    # Connect to Neo4j
    load_dotenv()
    URI = os.getenv("NEO4J_URI")
    USERNAME = os.getenv("NEO4J_USER")
    PASSWORD = os.getenv("NEO4J_PASSWORD")
    driver = GraphDatabase.driver(URI, auth=(USERNAME, PASSWORD))

    # Get data from MaintIE dataset
    with open('../data/MaintIE/gold_release.json', 'r', encoding='utf-8') as file:
        data = json.load(file)

    # Process data
    with driver.session() as session:
        session.execute_write(create_graph, data)

    driver.close()
