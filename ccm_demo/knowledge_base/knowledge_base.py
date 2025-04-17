import pandas as pd

from sqlalchemy import MetaData, select
from sqlalchemy.orm import Session

from ccm_demo.literature.literature import Paper

class KnowledgeBase:
    def __init__(self, engine):
        self.engine = engine
        self.meta = MetaData(bind=self.engine)
        self.meta.reflect(bind=self.engine)
        self.session = Session(self.engine)
        self.db_tables = self.meta.tables

    def add_paper(self, paper):
        pass

    def remove_paper(self, paper_id):
        pass

    def query(self, **kwargs):
        pass

    def question(self, question):
        pass

    def RAG(self, query):
        pass

    #def chat not sure, depends on the overall model

