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
        new_paper = Papers(
            source_id = paper.id,
            source = paper.source,
            title = paper.title,
            abstract = paper.abstract,
            pdf_url = paper.download_link,
            text = paper.text,
            figures = paper.figures,
            tables = paper.tables
        )
        self.session.add(new_paper)
        self.session.commit()

    def remove_paper(self, source_id):
        paper_to_delete = session.query(Papers).filter_by(source_id).one_or_none
        if paper_to_delete:
            session.delete(paper_to_delete)
        else:
            pass

    def query(self, **kwargs):
        pass

    def question(self, question):
        pass

    def RAG(self, query):
        pass

    #def chat not sure, depends on the overall model

