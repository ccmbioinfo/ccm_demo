import requests
import json
import os

class Protein:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.json = self.gather()
        self.process_uniprot()
        self.get_description()
        if "citationCrossReferences" in self.json:
            self.references = self.extract_references()
        else:
            self.references = []

        self.cross_references = self.get_xrefs()

    def gather(self):
        url = "https://rest.uniprot.org/uniprotkb/{}?format=json".format(self.uniprot_id)
        response = requests.get(url)
        response.raise_for_status()

        content = response.content.decode().strip()
        content = json.loads(content)
        return content

    def process_uniprot(self):
        self.sequence = self.json["sequence"]["value"]
        self.organism = {"name": self.json["organism"]["scientificName"],
                         "taxid": self.json["organism"]["taxonId"]}
        self.name = self.json["proteinDescription"]["recommendedName"]["fullName"]["value"]
        self.genes = [item["geneName"]["value"] for item in self.json["genes"]]
        self.comments = self.json["extraAttributes"]["countByCommentType"]
        self.features = self.json["extraAttributes"]["countByFeatureType"]

    def extract_features(self, feature_types=None):
        if feature_types is not None:
            features = [feat for feat in self.json["features"] if feat["type"] in feature_types]
        else:
            features = [feat for feat in self.json["features"]]

        return features

    def get_comments(self, types=None):
        if types is not None:
            if type(types) == str:
                types = [types]
            comments = [feat for feat in self.json["comments"] if feat["commentType"] in types]
        else:
            comments = [feat for feat in self.json["comments"]]
        return comments

    def extract_references(self):
        refs = []
        for reference in self.json["references"]:
            ref = {"title": reference["citation"]["title"],
                   "sources": reference["citation"]["citationCrossReferences"]}
            refs.append(ref)
        self.references = refs
        return self

    def get_xrefs(self):
        xrefs = self.json["uniProtKBCrossReferences"]
        self.xref_types = list(set([item["database"] for item in xrefs]))
        self.xrefs = xrefs
        return self

    def get_cross_references(self, types=None):
        if type(types) is str:
            types = [types]

        if type is None:
            return self.xrefs
        else:
            return [item for item in self.xrefs if item["database"] in types]


    def get_description(self):
        desc = []
        for comment in self.json["comments"]:
            if "texts" in comment.keys():
                desc.append("\n".join([item["value"] for item in comment["texts"]]))

        self.description = "\n".join(desc)
        return self

    def __str__(self):
        return self.description

    def __repr__(self):
        return "Protein instance with name {} and {} aa long".format(self.name, len(self.sequence))

class ProteinNetwork:
    def __init__(self, name, species=9606): 
        # name parameter supports use of uniprot id, uniprot accession #, gene name, or gene name synonyms
        # default organism is homo sapiens 

        self.species = species
        self.string_id, self.common_name = self.get_identifiers(name)
        self.interactions = self.get_interactions(common_name)
       

    def get_identifiers(self, name):
        response = requests.get(f"https://string-db.org/api/json/get_string_ids?identifiers={self.name}&species={self.species}")
        response.raise_for_status()

        content = response.content.decode().strip()
        content = json.loads(content)

        string_id = content[0]["stringId"]
        common_name = content[0]["preferredName"]
        return string_id, common_name
    
    def get_interactions(self, common_name): 
        response = requests.get(f"https://string-db.org/api/json/network?identifiers={self.common_name}&species={self.species}")
        response.raise_for_status()

        content = response.content.decode().strip()
        content = json.loads(content)

        # len(content) = number_of_interactors

        # for item in content: 
        #     string_id_A = item["stringId_A"]
        #     string_id_B = item["stringId_B"] 
        #     common_name_A = item["preferredName_A"]
        #     common_name_B = item["preferredName_B"]
        #     score = item["score"]
        #     nscore = item["nscore"]
        #     fscore = item["fscore"]
        #     pscore = item["pscore"]
        #     ascore = item["ascore"]
        #     escore = item["escore"]
        #     dscore = item["dscore"]
        #     tscore = item["tscore"]
        
        return content

        
    def get_network(self):




        
       

        
