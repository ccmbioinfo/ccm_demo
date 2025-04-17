from bs4 import BeautifulSoup as bs
import requests
import tarfile

from ccm_demo.literature.utils import *

class NoPapersError(Exception):
    pass

class LitSearch:
    def __init__(self, pubmed_api_key=None):
        """
        create the ncessary framework for searching
        :param pubmed_url:
        :param arxiv_url:
        :param pubmed_api_key:
        """
        self.pubmed_key = pubmed_api_key

    #TODO advanced search
    def search(self, query, database="pubmed", results="id", max_results=1000):
        """
        search pubmed and arxiv for a query, this is just keyword search no other params are implemented at the moment
        :param query:
        :param database:
        :param results: what to return, default is paper id PMID and arxiv id
        :param max_results:
        :return:
        """
        #TODO implement pubmed api key for non-free papers
        if database == "pubmed":
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={}&retmax={}".format(query, max_results)
            search_response = requests.get(search_url)
            search_response.raise_for_status()

            soup = bs(search_response.text, "xml")
            ids = [item.text for item in soup.find_all("Id")]

            if results == "doi":
                dois = []
                for paperid in ids:
                    response = requests.get(
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}".format(paperid))
                    response.raise_for_status()
                    soup = bs(response.text, "xml")
                    dois.append([item.text for item in soup.find_all("ArticleId") if item.attrs["IdType"] == "doi"])
                to_ret=dois
            else:
                to_ret=ids

        elif database == "arxiv":
            search_url="http://export.arxiv.org/api/search_query?{}&max_results={}".format(query, str(max_results))
            search_response = requests.get(search_url)
            search_response.raise_for_status()
            soup = bs(search_response.text, "xml")
            ids=[item.text.split("/").pop() for item in soup.find_all("id")][1:] #first one is the search id
            to_ret= ids
        return to_ret

class Paper:
    def __init__(self, id):
        self.id = id
        self.title = None
        self.abstract = None
        self.text = None
        self.figures = None
        self.tables = None

    def get_meta(self, source):
        if source =="pubmed":
            response=requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}".format(self.id))
            response.raise_for_status()
            soup=bs(response.text, "xml")
            self.title = soup.find("ArticleTitle").text
            self.abstract = soup.find("AbstractText").text
            self.source=source
            ids=[item for item in soup.find_all("ArticleIdList") if item.name=="PubmedData"]
            pmc=[item.text for item in list(ids[0].children) if item["IdType"] == "pmc"][0]
            pmc_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={}".format(pmc)
            response = requests.get(pmc_url)
            response.raise_for_status()
            self.download_link=bs(response.text, "xml").find_all("link")[0]["href"]

        elif source == "arxiv":
            response = requests.get("http://export.arxiv.org/api/query?search_query=id:{}".format(self.id))
            response.raise_for_status()
            soup=bs(response.text, "xml")
            self.abstract = soup.find("summary").text
            if self.abstract is None:
                raise NoPapersError("could not find a paper with id {}".format(self.id))
            else:
                self.title = soup.find("title").text
                self.source=source
                pdf=[item["href"] for item in soup.find_all("link") if item.attrs["type"] == "application/pdf"]
                self.download_link = pdf
        else:
            raise NotImplementedError("source must be pubmed or arxiv other sources are not implemented")

        return self

    def download(self, destination):
        if self.source == "pubmed":
            download = requests.get(self.download_link, stream=True)
            download.raise_for_status()
            with open("{}/{}.tar.gz".format(destination, self.id), "wb") as f:
                f.write(download.content)
            #TODO extract pdf and cleanup
        elif self.source == "arxiv":
            download = requests.get(self.download_link, stream=True)
            download.raise_for_status()
            with open("{}/{}.pdf".format(destination, self.id), "wb") as f:
                f.write(download.content)

    def process(self, pdf):
        text, figures, captions= process_pdf(pdf)
        text=[cleanup_text(item) for item in text]
        self.text=text
        self.figures=figures






