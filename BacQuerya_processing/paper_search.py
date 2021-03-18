#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download information of interest available through NCBI.
"""
from Bio import Entrez
from tqdm import tqdm
from urllib.parse import quote

def search_pubmed(searchTerm,
                  email,
                  number):
    """Download the attribute for accessions of interest using Biopython Entrez"""
    Entrez.email = email
    handle = Entrez.read(Entrez.esearch(db="pubmed", term=searchTerm, retmax = number))
    idList = handle["IdList"]
    searchResult = []
    for idTerm in tqdm(idList):
        try:
            esummary_handle = Entrez.esummary(db="pubmed", id=idTerm, report="full")
            esummary_record = Entrez.read(esummary_handle, validate = False)
            encodedDOI = quote(esummary_record[0]["DOI"], safe='')
            esummary_record[0]["encodedDOI"] = encodedDOI
            searchResult += esummary_record
        except:
            pass
    return searchResult

if __name__ == '__main__':
    search_pubmed("", "email", 100)