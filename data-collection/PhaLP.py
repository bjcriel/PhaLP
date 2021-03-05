#############
###IMPORTS###
#############
import requests
from Bio import SeqIO, Entrez, SwissProt, Medline, pairwise2
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import joblib
from http.client import IncompleteRead
from urllib.error import URLError, HTTPError
from requests.exceptions import ConnectionError, ChunkedEncodingError, Timeout
from itertools import zip_longest
from difflib import SequenceMatcher
import re
import mysql.connector
import pymysqlpool
import pandas
from io import StringIO
import numpy
import collections
import csv
from pathlib import Path
import os
import sys
from setuptools.dist import sequence
from copy import deepcopy
import subprocess
from Bio.SeqIO import FastaIO
from datetime import datetime
from inspect import currentframe, getframeinfo
from time import sleep
import json
import ijson
from ftplib import FTP
import zipfile
from Bio.ExPASy import Enzyme
from bio_embeddings.embed import SeqVecEmbedder
import gc
import shelve
import concurrent.futures

###############
###FUNCTIONS###
###############
class Logger(object):

    def __init__(self,):
        self.terminal = sys.stdout
        starttime = datetime.now().strftime("%Y%m%d-%H:%M:%S")
        self.log = open("../log/%s.log"%starttime, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass  

def get_url():
    ## GENERATE QUERY
    url_base = 'https://www.uniprot.org/uniprot/?query='
    
    
    tax_list = ['Caudovirales [28883]',
                'unclassified bacterial viruses [12333]',
                'Tectiviridae [10656]',
                'Cystoviridae [10877]',
                'Microviridae [10841]',
                'Corticoviridae [10659]', 
                'Sphaerolipoviridae [1714267]',
                'Picobirnaviridae [585893]',
                'Leviviridae [11989]'
                ]
                
    GO_list = ['catabolism by organism of cell wall peptidoglycan in other organism [51672]',
               'lytic transglycosylase activity [8933]',
               'peptidoglycan beta-N-acetylmuramidase activity [33922]',
               'N-acetylmuramoyl-L-alanine amidase activity [8745]',
               'lysozyme activity [3796]',
               'peptidoglycan binding [0042834]',
               'peptidoglycan catabolic process [9253]',
               'peptidoglycan N-acetylglucosaminidase activity [61784]',
               'peptidoglycan endopeptidase activity [61785]',
               'disruption by virus of host cell wall peptidoglycan during virus entry [98932]',
               'peptidoglycan muralytic activity [61783]',
               'peptidoglycan stem peptide endopeptidase activity [61786]',
               'peptidoglycan cross-bridge peptide endopeptidase activity [61787]'
               ]    
    
    IPR_list = ['IPR000064', # NLP_P60_dom #745
                'IPR000726', # Glyco_hydro_19_cat #278 many chitinases
                'IPR001165', # T4-type lysozyme #17
                'IPR002053', # Glyco_hydro_25 # 15
                'IPR002196', # Glyco_hydro_24 #65
                'IPR002502', # N-acetylmuramoyl-L-alanine amidase domain # Ami_2 # 79
                'IPR002508', # N-acetylmuramoyl-L-alanine amidase, catalytic domain #5
                'IPR002901', # Mannosyl-glycoprotein endo-beta-N-acetylglucosamidase-like domain # 481
                'IPR003709', # Peptidase M15B (Pfam VanY) #55
                'IPR005490', # YkuD #4
                'IPR006619', # Peptidoglycan recognition protein family domain, metazoa/bacteria #26
                'IPR007921', # CHAP_dom #781
                'IPR008044', # Phage_lysin #237
                'IPR008258', # Transglycosylase_SLT_dom_1 #941
                'IPR008565', # DUF847 (Pfam GH_108) #245
                'IPR010618', # Transglycosylas # 203
                'IPR016284', # Bacteriophage PRD1, P15, lysozyme #6
                'IPR013207', # LGFP #21
                'IPR013230', # Peptidase_M15A_C #245
                'IPR015020', # DUF1906; overlaps glycohydro family #9
                'IPR015510', # PGRP #28
                'IPR016047', # Peptidase_M23 #437
                'IPR018077', # Glycoside hydrolase, family 25 subgroup #5
                'IPR018537', # PGB_3
                'IPR019505', # Peptidase_U40 #37
                'IPR021976', # Amidase02_C # 2
                'IPR022016', # DUF3597 #1
                'IPR024408', # N-acetylmuramidase #164
                'IPR033907', # Endolysin/autolysin #34
                'IPR034689', # Endolysin T7 type #2
                'IPR034690', # Endolysin T4 type #1
                'IPR034691', # Endolysin lambda type #1
                'IPR036505', # PGRP domain superfamily #110
                'IPR038994', # Internal virion protein Gp16 #2
                'IPR039561', # Peptidase_M15C #724
                'IPR039564', # Peptidase_C39-like #251
                'IPR041219', # Phage tail lysozyme #
                'IPR023346' # Lysozyme-like 1009
               ] 

    query_tax = '+OR+'.join([ 'taxonomy:%22' + i.replace(' ','+') + '%22' for i in tax_list])
    query_GO = '+OR+'.join([ 'goa:(%22' + i.replace(' ','+') + '%22)' for i in GO_list])
    query_IPR = '+OR+'.join([ 'database:(type:interpro+' + i + ')' for i in IPR_list])
    url_html = url_base + '(' + query_tax + ')' + '+AND+' + '((' + query_GO + ')+OR+(' + query_IPR+'))'
    return url_html

def read_joblib(name):
    with open(name, 'rb') as f:
        d = joblib.load(f)
        return d

def save_joblib(d, name):
    with open(name, 'wb') as f:
        joblib.dump(d, f)

def save_json(d, name):
    with open(name, 'w') as f:
        json.dump(d, f, indent=4)

def update_json(key, value, name):
    with open(name, 'r+') as f:
        d = json.load(f)
        d.update({key:value})
        json.dump(d, f, indent=4)
        
def read_json(name):
    with open(name, 'r') as f:
        d = json.load(f)
        return d

def download_VirusHostDB():
    if not os.path.exists('downloads'):
        os.makedirs('downloads')
        q
    ftp = FTP('ftp.genome.jp')
    ftp.login()
    ftp.cwd('/pub/db/virushostdb/')
    ftp.retrbinary('RETR virushostdb.tsv', open('downloads/virushostdb.tsv', 'wb').write)
    ftp.quit()
    with open('downloads/virushostdb.tsv', 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        VH_dict={}

        for line in reader:
            if 'Bacteria' in line[9].split('; '):
                host = (line[7], line[8]) #tuple (host_id, host_name)
                VH_dict.setdefault(line[0],[]).append(host)
        with shelve.open("picklejar/VH_shelf.db", writeback = True) as VH_shelf:
            for key, value in VH_dict.items():
                VH_shelf[key]=value
            VH_shelf.sync()
            save_json(list(VH_shelf.keys()), 'picklejar/VH_shelf_keys.json')
        save_json(list(VH_dict.keys()), 'picklejar/VH_dict.json')
    return

def download_NCBI_taxonomy():
    if not os.path.exists('downloads'):
        os.makedirs('downloads')
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/pub/taxonomy/new_taxdump/')
    with open('downloads/new_taxdump.zip', 'wb') as f:
        ftp.retrbinary('RETR new_taxdump.zip', f.write)
    ftp.quit()
    with zipfile.ZipFile('downloads/new_taxdump.zip','r') as z:
        z.extractall('downloads/new_taxdump/')
    # update_taxonomy_dictionaries
    rank={}
    ranks = []
    with open('downloads//new_taxdump/nodes.dmp', 'rt', encoding="utf-8") as f:
        lines = csv.reader(f, delimiter='|')
        for line in lines:
            line = [i.strip('\t') for i in line]
            rank[line[0]]=line[2]
            if line[2] not in ranks:
                ranks.append(line[2])              
    tax_lineage_dict = {}
    with open('downloads//new_taxdump/rankedlineage.dmp', 'rt', encoding="utf-8") as f:
        lines = csv.reader(f, delimiter='|')
        for line in lines:
            line = [i.strip('\t') for i in line]
            species, genus, family, order, classs, phylum, kingdom, superkingdom = line[2:-1]
            tax_lineage_dict[line[0]]={'species':species , 'genus':genus , 'family':family, 'order':order, 'class':classs, 'phylum':phylum, 'kingdom':kingdom , 'superkingdom':superkingdom}
            if rank[line[0]] in tax_lineage_dict[line[0]].keys():
                tax_lineage_dict[line[0]][rank[line[0]]]=line[1]
    #make dictionaries for all entries in the NCBI taxonomy DB
    #from tax to name(scientific only)
    tax_name_dict = {} 
    #from name (all, incl. synonyms) to taxID
    name_tax_dict = {} 
    with open('downloads/new_taxdump/names.dmp', 'rt', encoding="utf-8") as f:
        lines = csv.reader(f, delimiter='|')
        for line in lines:
            line = [i.strip('\t') for i in line]
            
            if tax_lineage_dict[line[0]]['superkingdom'] == 'Bacteria':
                bacterial = True
            else:
                bacterial = False
            #if there are multiple organisms (different taxID) with the same name, append the tax_IDs in a string
            if line[1] in name_tax_dict.keys():               
                name_tax_dict[line[1]]['tax_ID'].append(line[0])
                name_tax_dict[line[1]]['bacterial'].append(bacterial)
            else:
                name_tax_dict[line[1]]= {'tax_ID':[line[0]],'bacterial':[bacterial]}
            if line[3]=='scientific name':
                if line[0]in tax_name_dict.keys():
                    raise ValueError('%s already in use for other organism'%line[0]) 
                tax_name_dict[line[0]]=line[1]            
    save_joblib(tax_lineage_dict,'picklejar/taxID_to_lineage.joblib')
    save_joblib(tax_name_dict,'picklejar/taxID_to_name.joblib')
    save_joblib(name_tax_dict,'picklejar/name_to_taxID.joblib')
    return

def download_Expasy_enzyme():
    if not os.path.exists('downloads'):
        os.makedirs('downloads')
        
    ftp = FTP('ftp.expasy.org')
    ftp.login()
    ftp.cwd('/databases/enzyme')
    with open('downloads/enzyme.dat', 'wb') as f:
        ftp.retrbinary('RETR enzyme.dat', f.write)
    with open('downloads/enzclass.txt', 'wb') as f:
        ftp.retrbinary('RETR enzclass.txt', f.write)
    ftp.quit()

    # get enzyme classes
    classes_dict={}
    handle = open('downloads/enzclass.txt')
    for line in handle:
        if re.match(r'\d\.', line):
            ID = line[:9].replace(' ','')
            description = re.compile(r'\s\s+').split(line)[-1][:-1]
            classes_dict[ID]=description
            
    # update EC dict
    EC_dict = {}
    handle = open('downloads/enzyme.dat')
    
    for rec in Enzyme.parse(handle):
        EC_ID = rec['ID']
        EC_dict[EC_ID]={'name':rec['DE'], 'alt_names':'; '.join(rec['AN']), 'reaction':rec['CA']}
    
        while True:
            class_ID = '.'.join([EC_ID.split('.')[0],'-','-','-'])
            EC_dict[EC_ID]['class_ID']=class_ID
            EC_dict[EC_ID]['class_name']=classes_dict[class_ID]
            
            if EC_ID.split('.')[1]=='-':
                EC_dict[EC_ID]['subclass_ID']=None
                EC_dict[EC_ID]['subclass_name']=None
                EC_dict[EC_ID]['subsubclass_ID']=None
                EC_dict[EC_ID]['subsubclass_name']=None
                break
            subclass_ID = '.'.join([EC_ID.split('.')[0],EC_ID.split('.')[1],'-','-'])
            EC_dict[EC_ID]['subclass_ID']=subclass_ID
            EC_dict[EC_ID]['subclass_name']=classes_dict[subclass_ID]
            
            if EC_ID.split('.')[2]=='-':
                EC_dict[EC_ID]['subsubclass_ID']=None
                EC_dict[EC_ID]['subsubclass_name']=None
                break
            subsubclass_ID = '.'.join([EC_ID.split('.')[0],EC_ID.split('.')[1],EC_ID.split('.')[2],'-'])
            EC_dict[EC_ID]['subsubclass_ID']=subsubclass_ID
            EC_dict[EC_ID]['subsubclass_name']=classes_dict[subsubclass_ID]
            break
                
    for EC_ID, value in classes_dict.items():
        
        EC_dict[EC_ID]={'name':None, 'alt_names':None, 'reaction':None}
        while True:
            class_ID = '.'.join([EC_ID.split('.')[0],'-','-','-'])
            EC_dict[EC_ID]['class_ID']=class_ID
            EC_dict[EC_ID]['class_name']=classes_dict[class_ID]
            
            if EC_ID.split('.')[1]=='-':
                EC_dict[EC_ID]['subclass_ID']=None
                EC_dict[EC_ID]['subclass_name']=None
                EC_dict[EC_ID]['subsubclass_ID']=None
                EC_dict[EC_ID]['subsubclass_name']=None
                break
            subclass_ID = '.'.join([EC_ID.split('.')[0],EC_ID.split('.')[1],'-','-'])
            EC_dict[EC_ID]['subclass_ID']=subclass_ID
            EC_dict[EC_ID]['subclass_name']=classes_dict[subclass_ID]
            
            if EC_ID.split('.')[2]=='-':
                EC_dict[EC_ID]['subsubclass_ID']=None
                EC_dict[EC_ID]['subsubclass_name']=None
                break
            subsubclass_ID = '.'.join([EC_ID.split('.')[0],EC_ID.split('.')[1],EC_ID.split('.')[2],'-'])
            EC_dict[EC_ID]['subsubclass_ID']=subsubclass_ID
            EC_dict[EC_ID]['subsubclass_name']=classes_dict[subsubclass_ID]
            break        
    with shelve.open("picklejar/EC_shelf.db", writeback = True) as EC_shelf:
        for key, value in EC_dict.items():
            EC_shelf[key]=value
        EC_shelf.sync()
    save_json(EC_dict, 'picklejar/EC_dict.json')

def update( last_PhaLP_version, PhaLP_version,url_html):
    global analyzed, VH_shelf_keys, tax_lineage_dict,tax_name_dict,name_tax_dict,name_tax_dict_lowercase, EC_dict, canonical2db_name, dbname2canonical, canonical2version_date, embedder, clf, SeqVecEmbeddings_index, KfoldCVpredictions
    connect_mysql('user','password','host','port','database') # fill your own mysqlconnection parameters
    proteins_in_PhaLP = retrieve_protein_acc_in_PhaLP()
    temporary_UniProt = 'temporary_uniprot'
    Entrez.email = 'fill.your@email.com' #fill your own email
    retrieve_UniProt(url_html, temporary_UniProt,'txt')
    num_analyzed = 0
    canonical2db_name, dbname2canonical, canonical2version_date = retrieve_InterPro_version_info()
    # retrieve required dictionaries from pickle file
    VH_shelf_keys = read_json('picklejar/VH_shelf_keys.json')
    tax_lineage_dict = read_joblib('picklejar/taxID_to_lineage.joblib')
    #Taxonomy ID:sientific name
    tax_name_dict = read_joblib('picklejar/taxID_to_name.joblib')
    #all names(including synonyms):Taxonomy ID
    name_tax_dict = read_joblib('picklejar/name_to_taxID.joblib')
    name_tax_dict_lowercase = {k.lower():v for k,v in name_tax_dict.items()}
    EC_dict = read_json('picklejar/EC_dict.json')
    embedder = SeqVecEmbedder()
    clf = read_joblib('picklejar/RF_clf_SeqVecEmbeddings_trained.pickle')
    SeqVecEmbeddings_index = list(pandas.read_csv('picklejar/SeqVecEmbeddings.csv', index_col = 0).index)
    KfoldCVpredictions=pandas.read_csv('picklejar/10foldCVpredictions_SeqVecEmbeddings.csv', index_col = 0)

    if Path('picklejar/analyzed.json').exists():
        a = read_json('picklejar/analyzed.json')
        if PhaLP_version not in a.keys():
            analyzed={'phalps':[],'phages':[],'UPIs':[]}
        else:
            analyzed = a[PhaLP_version]
    else:
        analyzed = {'phalps':[],'phages':[],'UPIs':[]}
        save_json(analyzed, 'picklejar/analyzed.json')

    with open(temporary_UniProt+'.txt', 'rt') as f:
        num_total = sum([1 for i in SwissProt.parse(f)])

    while len(analyzed['phalps'])< num_total:
        try:
            with open(temporary_UniProt+'.txt', 'rt') as f, open(temporary_UniProt+'.txt', 'rt') as f2:
                with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                    swiss_recs = []
                    seq_recs = []
                    for swiss_rec, seq_rec in zip_longest(SwissProt.parse(f), SeqIO.parse(f2, 'swiss')):
                        if seq_rec.id not in analyzed['phalps']:
                            swiss_recs.append(swiss_rec)
                            seq_recs.append(seq_rec)
                        else:
                            num_analyzed += 1
                    print('%s PhaLPs out of %s already analyzed'%(num_analyzed, num_total))
                    results = executor.map(analyze_PhaLP, swiss_recs, seq_recs)
                    for ACC, PHAGE_ID, UPI, submit_check in results: 
                        if ACC not in analyzed['phalps']:
                            print('\tend ', ACC)
                            #add accession to the list of analyzed phalps
                            if False not in [value for key, value in submit_check.items()]:
                                analyzed['phalps'].append(ACC)
                                if UPI not in  analyzed['UPIs']:
                                    analyzed['UPIs'].append(UPI)
                                if PHAGE_ID not in  analyzed['phages']:
                                    analyzed['phages'].append(PHAGE_ID)
                            else:
                                for key, value in submit_check.items():
                                    if value == False:
                                        print('Submission of PhalP %s not completed because of problem with %s submit'%(ACC, key))
                        else:
                            print('\tfixed ', ACC)
                        gc.collect()
                        num_analyzed += 1
                        if num_analyzed%10 == 0:
                            print('%s PhaLPs out of %s analyzed'%(num_analyzed, num_total))
                            update_json(PhaLP_version, analyzed,'picklejar/analyzed.json')
                    
                    print('%s PhaLPs out of %s analyzed'%(num_analyzed, num_total))
                    update_json(PhaLP_version, analyzed,'picklejar/analyzed.json')
               
               #Remove UniProt entries that are not in query anymore
                for acc in proteins_in_PhaLP:
                    if acc not in analyzed['phalps']:
                        print (acc, 'must be removed')
                        delete_UniProt_entry(acc)     
                # Remove entries from parent tables that have no children anymore
                clean_up_parent_tables()
        except Exception as e:
            print('Starting new attempt, error occured: ', e)
            
    return analyzed
    
def analyze_PhaLP(swiss_rec, seq_rec):
    assert swiss_rec.accessions[0]== seq_rec.id  
    ACC = seq_rec.id
    last_annotation_update = datetime.strptime(seq_rec.annotations['date_last_annotation_update'],'%d-%b-%Y')
    # if the annotation of the gb-file has been changed since the last PhaLP update, REPLACE existing entries from the protein with new ones.                       
    submit_check = {}
    print('start ', ACC)
    DATE_CREATED = datetime.strptime(seq_rec.annotations['date'], '%d-%b-%Y') #format .strftime("%Y-%m-%d")
    DATE_LAST_UPDATED = datetime.strptime(seq_rec.annotations['date_last_annotation_update'], '%d-%b-%Y') #format .strftime("%Y-%m-%d")
    #find all the phage info   
    PHAGE_NAME, PHAGE_ID, PHAGE_LINEAGE = extract_phage_taxonomy(seq_rec)
    PROTEIN_NAME = extract_protein_name(seq_rec.description)
    # Try to submit the phage. If it already exists in the database, the hosts are not analyzed. Otherwise, find the hosts and put them in the db.
    if PHAGE_ID not in analyzed['phages']:
        submit_check['phage'] = submit_phage(PHAGE_ID, PHAGE_NAME, PHAGE_LINEAGE)
        #Remove all previous phage-host links
        HOSTS = extract_host_taxonomy(seq_rec, PHAGE_ID, PHAGE_NAME)
        clean_link_phage_host(PHAGE_ID)
        submit_check['host'] = True
        for host in HOSTS:
            submit_check_host_temp = submit_host(host, PHAGE_ID)
            if submit_check_host_temp == False:
                submit_check['host'] = False
    #Find the stable UniParc identifier (unique protein sequence) for the protein sequence
    UPI = ID_mapping(ACC,'ID','UPARC')
    # if the protein sequence (UniParcID) is not yet in the database, analyze it and put it in
    if UPI not in analyzed['UPIs']:
        #analyze protein sequence from representative
        SEQ, LEN, Mw, pI, Aromaticity, Hydropathy  = analyze_protein_sequence(seq_rec)
        submit_check['UniParc'] = submit_UniParc(UPI, SEQ, LEN, Mw, pI, Aromaticity, Hydropathy )
        DOMAINS = extract_domains(ACC, seq_rec)
        clean_link_UniParc_domains(UPI)
        submit_check['domains'] = submit_domains(DOMAINS, UPI)
        domain_accs = DOMAINS.keys()
        TYPE, TYPE_EVIDENCE, TYPE_PROBABILITY = determine_type(UPI, domain_accs, SEQ)
    else:
        select_accs = ("SELECT UniProt_ID, type, type_evidence, type_probability FROM UniProt where UniRef_ID = '%s' ;"%UPI)
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(select_accs)
        data = cursor.fetchall()
        cursor.close()
        cnx.close()
        print('start analyzing %s'%ACC)
        for row in data:
            print(ACC, row)
            if row[0] != ACC and row[0] in analyzed['phalps']:
                TYPE = row[1]
                TYPE_EVIDENCE = row[2]
                TYPE_PROBABILITY = row[3]
                break
        submit_check['UniProt']= submit_UniProt(ACC, PROTEIN_NAME, TYPE, TYPE_EVIDENCE, TYPE_PROBABILITY, UPI, PHAGE_ID, DATE_CREATED, DATE_LAST_UPDATED)
         
        REF_EXP_EVID = extract_references_experimental_evidence(swiss_rec,ACC)
        clean_UniProt_child(ACC,'experimental_evidence')
        submit_check['exp_evid'] = submit_experimental_evidence(ACC, REF_EXP_EVID)
     
        GOs = extract_GOs(swiss_rec)
        clean_UniProt_child(ACC,'gene_ontologies')
        submit_check['GO'] = submit_GOs(ACC,GOs)
         
        ECs = extract_EC(seq_rec.description)
        clean_UniProt_child(ACC,'link_EC_UniProt')
        submit_check['EC'] = submit_EC(ECs, ACC)

        PDBs = extract_3D(ACC,swiss_rec)
        clean_UniProt_child(ACC,'tertiary_structures')
        submit_3D_structures(ACC, PDBs)

        CDSs = extract_CDSs(ACC,swiss_rec)
        clean_UniProt_child(ACC,'CDSs')
        submit_check['CDS'] = submit_CDSs(CDSs, ACC)
        return ACC, PHAGE_ID, UPI, submit_check
         
        
def connect_mysql(user,password,host,port,db):
    config = {'user': user, 'password': password, 'host': host, 'port': port, 'database': db, 'connect_timeout':1}
    #create a Pool of 4 MySQL connection so each Thread has its own connection 
    global cnxpool
    cnxpool = pymysqlpool.ConnectionPool(size=4, name='cnxpool', **config)

def retrieve_protein_acc_in_PhaLP(): 
    select_accs = ("SELECT UniProt_ID FROM UniProt;")
    cnx = cnxpool.get_connection()
    cursor = cnx.cursor()
    cursor.execute(select_accs)
    data = cursor.fetchall()
    cursor.close()
    cnx.close()
    proteins_PhaLP = []
    for row in data:
        proteins_PhaLP.append(row[0])
    return proteins_PhaLP

def retrieve_UniProt(url_html, temporary_UniProt, ext):
    ### Download XML file of all UniProt records that comply with the query (encoded in URL) and save ad temporary file
    # Build url with REST API: https://www.ebi.ac.uk/proteins/api/doc/
    # url = url from uniprot 'share' + '&format=xml'
    url = url_html+'&format=%s'%ext
    while True:
        try:
            r = requests.get(url)
            file = open(temporary_UniProt+'.%s'%ext, 'wt')
            file.write(r.text)
            file.close()
            break
        except (ConnectionError, ChunkedEncodingError) as err:
            #frameinfo= getframeinfo(currentframe())
            print('ERROR:\t',err)#, 'line %s, file'%(frameinfo.lineno()),frameinfo.filename()) 
    return
    
def retrieve_InterPro_version_info():
    canonical2db_name={}
    dbname2canonical={}
    canonical2version_date ={}

    dbname2canonical['GENE3D']='cathgene3d'
    canonical2db_name ['cathgene3d'] ='GENE3D'
    url = "https://www.ebi.ac.uk/interpro/api/"
    while True:
        try:
            r = requests.get(url)
            if r.status_code == 200:
                break
            if r.status_code == 408 or r.status_code ==404 or r.status_code == 500:
                print(r.status_code, 'wait 1 min', url)
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL                continue
            

            print('ERROR Interpro problem %s %s'%(r.status_code,r.text ))
           
        except (ConnectionError, ChunkedEncodingError, Timeout, HTTPError) as err:
            if err.code == 404:
                sleep(61)
                continue
            else:
                raise err
                #print('ERROR:\t%s\t'%ACC,err) 
    js = r.json()
    for value in js['databases'].values():
        canonical2db_name[value['canonical']]=value['name']
        dbname2canonical[value['name']]=value['canonical']
        if value['releaseDate']:
            datetimeobj = datetime.strptime(value['releaseDate'], '%Y-%m-%dT%H:%M:%SZ') #format .strftime("%Y-%m-%d")
        else:
            datetimeobj = None
        canonical2version_date[value['canonical']]=(value['version'],datetimeobj)

    return canonical2db_name, dbname2canonical, canonical2version_date

def extract_phage_taxonomy(seq_rec): 
    # extract phage name from UniProt file
    phage_name= seq_rec.annotations['organism'].split(' (')[0]
    phage_ID = seq_rec.annotations['ncbi_taxid'][0]
    # if tax_id is nonexistant in current ncbi taxonomy, lookup new tax_id for phage name
    if phage_ID not in tax_name_dict.keys():
        phage_ID = name_tax_dict[phage_name]['tax_ID'][name_tax_dict[phage_name]['bacterial'].index(False)]
    phage_lineage = tax_lineage_dict[phage_ID]
    return phage_name, phage_ID, phage_lineage

def extract_protein_name(description):
    names=description.split('; ')
    names[-1]=names[-1][:-1]            
    for name in names:
        if re.match(r'RecName:|SubName:', name):
            NAME = re.search(r'Full=(.*)', name).group(1).split(' {')[0]
            return NAME

def submit_phage(phage_ID, phage_name, phage_lineage):
    submit_check = True
    phage_already_analyzed = False
    for key, value in phage_lineage.items():
        if value == '':
            phage_lineage[key]=None
    submit_phage = ("INSERT INTO phages "
                   "(phage_name, phages_ID, lineage_order, lineage_family, lineage_genus, lineage_species) "
                   "VALUES (%s, %s, %s, %s, %s, %s) "
                   "ON DUPLICATE KEY UPDATE phage_name=VALUES(phage_name), phages_ID=VALUES(phages_ID), lineage_order=VALUES(lineage_order), lineage_family=VALUES(lineage_family), lineage_genus=VALUES(lineage_genus), lineage_species=VALUES(lineage_species)")
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(submit_phage, (phage_name, phage_ID, phage_lineage['order'], phage_lineage['family'], phage_lineage['genus'], phage_lineage['species']))
        cnx.commit()
        cursor.close()
        cnx.close()
    except mysql.connector.Error as err:
        if err.errno == 1062:
            #print('\n')
            #print('%s already exists in phages' %phage_name)
            None
        else:
            print('\n')
            print('ERROR',err)
            submit_check=False
    return submit_check

def extract_host_taxonomy(seq_rec, phage_ID, phage_name):
    hosts=[]
    #look in VirusHostDB
    if phage_ID in VH_shelf_keys:
        hosts = []
        with shelve.open("picklejar/VH_shelf.db", flag='r') as VH_shelf:
            for host in VH_shelf[phage_ID]:
                hosts.append({'host_ID':host[0],'host_name':host[1],'lineage':tax_lineage_dict[host[0]]})
    #look in UniProt file
    if 'host_ncbi_taxid' in seq_rec.annotations.keys():
        for i in range(len(seq_rec.annotations['host_ncbi_taxid'])):
            if seq_rec.annotations['host_ncbi_taxid'][i] not in [host['host_ID'] for host in hosts]:
                host_ID = seq_rec.annotations['host_ncbi_taxid'][i]
                host_name = seq_rec.annotations['organism_host'][0].split(' (')[0]
                if True not in name_tax_dict[host_name]['bacterial']:
                    break
                host_lineage = tax_lineage_dict[host_ID]
                hosts.append({'host_ID':host_ID,'host_name':host_name,'lineage':host_lineage})
    if not hosts:
        hosts.extend(extract_host_from_gb(seq_rec))
        if not hosts:
            extract_host_from_phage_name(phage_ID, phage_name)
    return hosts

def clean_link_phage_host(PHAGE_ID):
    clean_link_phage_host = ("DELETE FROM link_phage_host  WHERE link_phage_host.phages_ID =%s; ")
    
    try:
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_link_phage_host,(PHAGE_ID,))  
        cnx.commit()
        cursor.close()
        cnx.close()
        #print('hosts for phage', PHAGE_ID, 'cleaned',cursor.rowcount,'rows deleted')
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
            #print('\n')
            #print('%s already exists in UniProt' %ACC)
        else:
            print('\n')
            print('ERROR',err)
    return 

def extract_host_from_gb(seq_rec):    
    #find gb accession numbers (crosref) in uniprot file
    EMBL = []
    for ref in seq_rec.dbxrefs:   
        if ref.split(':')[0]=='EMBL':
            EMBL.append(ref.split(':')[1])
    while True:
        try:
            handle = Entrez.efetch(db='nucleotide', id=','.join(EMBL),rettype="gb", retmode="text")
            break
        except (URLError, IncompleteRead, ValueError, RuntimeError, TimeoutError) as err:
            print(err)
    host_names=[]
    # retrieve host names from gb files 
    for gb in SeqIO.parse(handle,'gb'):
        for feat in gb.features:
            if feat.type == 'source':
                if 'host' in feat.qualifiers and feat.qualifiers['host'][0] not in host_names:
                    #retrieve host information
                    host_names.append(str(feat.qualifiers['host'][0]))
                if 'lab_host' in feat.qualifiers and feat.qualifiers['lab_host'][0] not in host_names:
                    host_names.append(str(feat.qualifiers['lab_host'][0]))
    # if host name(s) is/are found, check wether it is a known name in the ncbi taxonomy and whether it is a bacterium or not. 
    if len(host_names)>0:
        for host_name in list(host_names):
            host_name_new = host_name
            
            while host_name_new.lower() not in [x.lower() for x in name_tax_dict.keys()]:
                
                if host_name_new.replace(',','') in name_tax_dict.keys():
                    host_name_new = host_name_new.replace(',','')
                    break
                host_name_new = cut_last_word(host_name_new)
                #reset is an empty string is reached
                if host_name_new == '':
                    host_name_new = host_name
                    break
                
            while host_name_new.lower() not in [x.lower() for x in name_tax_dict.keys()]:
                host_name_new = cut_first_word(host_name_new)
                # check if the name was not misspelled
                if host_name_new == '':
                    host_name_new = correct_typo_annotation(host_name)                      
                    break
    
            
            # check if name is known in ncbi taxonomy (either as scientific name or as synonym
            if host_name_new and host_name_new.lower() in [x.lower() for x in name_tax_dict.keys()]:
                #check whether the host is a bacterium, otherwise remove
                if True not in name_tax_dict_lowercase[host_name_new.lower()]['bacterial']:
                    host_names.remove(host_name)
                    print(host_name,'removed. Not bacterial')
                #if the name is not the scientific name, but a synonym, replace by the scientific name
                elif host_name not in tax_name_dict.values():
                    host_names[host_names.index(host_name)] = tax_name_dict[name_tax_dict_lowercase[host_name_new.lower()]['tax_ID'][0]]
                else:
                    host_names[host_names.index(host_name)] = host_name_new
            else:
                host_names.remove(host_name)
    if len(host_names)>0:
        hosts = []
        for host_name in host_names:
            hosts.append({'host_ID':name_tax_dict[host_name]['tax_ID'][0],'host_name':host_name,'lineage':tax_lineage_dict[name_tax_dict[host_name]['tax_ID'][0]]}) 
        return hosts
    else:
        return []

def cut_last_word(host_name):

    return ' '.join( host_name.split(' ')[:-1])

def cut_first_word(host_name):

    return ' '.join( host_name.split(' ')[1:])

def correct_typo_annotation(host_name):
    host_name_new = host_name
    print("Problem with '%s'. Checking annotation for typos" %host_name_new)
    while True:
        options = []
        for name in name_tax_dict.keys():
            if True in name_tax_dict[name]['bacterial']:
                similarity = string_similarity(name, host_name_new)
                if similarity > 0.9 and (len(host_name_new.split(' '))==len(name.split(' '))):
                    #print(name, similarity, host_name_new)
                    options.append({'name':name,'tax_ID':name_tax_dict[name]['tax_ID'][0],'similarity':similarity})
        if options:
            if len(list(set([x['tax_ID'] for x in options]))) == 1:
                return tax_name_dict[options[0]['tax_ID']]
            else:
                print('Multiple options found:\t',options)
                return None
        else:
            host_name_new = cut_last_word(host_name_new)
        if host_name_new == '':
            print('No alternative names found')
            return None

def string_similarity(a, b):

    return SequenceMatcher(None, a, b).ratio()

def extract_host_from_phage_name(phage_ID, phage_name):
    # find host in phage name
    phage_synonyms = []
    for key, value in name_tax_dict.items():
        if value['tax_ID'][0] == phage_ID:
            phage_synonyms.append(key)
    host_ID_options = []
    for phage_synonym in phage_synonyms:
        for x in ['phage','bacteriophage', 'virus']:
            if x in phage_synonym.split(' '):
                host_option = phage_synonym.split(' '+x)[0]
                if host_option.lower() in name_tax_dict_lowercase.keys():
                    d = name_tax_dict_lowercase[host_option.lower()]
                    if True in d['bacterial']:
                        host_ID = d['tax_ID'][d['bacterial'].index(True)]
                        host_ID_options.append(host_ID)    
    if len(host_ID_options)>1:
        rank = []
        ranks = ['species', 'genus', 'family', 'order', 'classs', 'phylum', 'kingdom', 'superkingdom']
        for host_ID in host_ID_options:
            for r in ranks:
                if tax_lineage_dict[host_ID][r]!='':
                    rank.append(ranks.index(r))
                    break
            
        host_ID =   host_ID_options[rank.index(min(rank))]
        hosts = [{'ID':host_ID,'name':tax_name_dict[host_ID],'lineage':tax_lineage_dict[host_ID]}]
    elif len(host_ID_options)==1:
        hosts = [{'ID':host_ID_options[0],'name':tax_name_dict[host_ID_options[0]],'lineage':tax_lineage_dict[host_ID_options[0]]}]

    else:
        print('no host:\t', phage_name)
        hosts = []
    return hosts

def submit_host(host, phage_ID):   
    submit_check_host = True 
    for key, value in host['lineage'].items():
        if value == '':
            host['lineage'][key]=None
    submit_host = ("INSERT INTO hosts"
                   "(host_name, hosts_ID, lineage_phylum, lineage_class, lineage_order, lineage_family, lineage_genus, lineage_species)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
                    "ON DUPLICATE KEY UPDATE host_name=VALUES(host_name), hosts_ID=VALUES(hosts_ID), lineage_phylum=VALUES(lineage_phylum), lineage_class=VALUES(lineage_class), lineage_order=VALUES(lineage_order), lineage_family=VALUES(lineage_family), lineage_genus=VALUES(lineage_genus), lineage_species=VALUES(lineage_species)")
    try:
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(submit_host,(host['host_name'], host['host_ID'], host['lineage']['phylum'], host['lineage']['class'], host['lineage']['order'], host['lineage']['family'], host['lineage']['genus'], host['lineage']['species']))
        cnx.commit()
        cursor.close()
        cnx.close()
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
            print('\n')
            print('%s already exists in hosts' %host['host_name'])
        else:
            print('\n')
            print('ERROR',err)
            submit_check_host =  False
    
    
    submit_phage_host = ("INSERT INTO `link_phage_host` "
                         "(phages_ID, hosts_ID) "
                         "VALUES (%s, %s) "
                         "ON DUPLICATE KEY UPDATE phages_ID=VALUES(phages_ID), hosts_ID=VALUES(hosts_ID)")
    try:
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(submit_phage_host,(phage_ID, host['host_ID']))
        cnx.commit()
        cursor.close()
        cnx.close()
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
            print('\n')
            print('%s - %s already exists in link_phage_host' %(phage_ID, host['host_ID']))
        else:
            print('\n')
            print('ERROR',err)
            submit_check_host = False
    return submit_check_host           

def ID_mapping(ID, db_from, db_to):
    mapping_url = 'https://www.uniprot.org/uploadlists/'
    mapping_params = {
        'from':'%s'%db_from,
        'to':'%s'%db_to,
        'format':'tab',
        'query':'%s'%ID
        }
    while True:
        try:
            response = requests.get(mapping_url, params=mapping_params)
            return response.text.split('\n')[1].split('\t')[1]
        except (ConnectionError, IndexError, requests.exceptions.TooManyRedirects) as err:
            print('ERROR:\t%s\t'%ID, err) 

def analyze_protein_sequence(seq_rec):
    SEQ = str(seq_rec.seq)
    LEN = len(SEQ)
    analysed_seq = ProteinAnalysis(SEQ)
    #if there is an ambiguous AA in the sequence remove it for the Mw calculation
    try:
        Mw = analysed_seq.molecular_weight()
        Hydropathy = analysed_seq.gravy()
    except ValueError as err:
        m = re.search(r"'(.)' is not a valid unambiguous letter for protein",str(err))
        ambiguous_AA = m.groups()[0]
        cut_seq = SEQ.replace(ambiguous_AA,'')
        analysed_cut_seq = ProteinAnalysis(cut_seq)
        Mw = analysed_cut_seq.molecular_weight() 
        Hydropathy = analysed_cut_seq.gravy()       
    pI = analysed_seq.isoelectric_point()
    Aromaticity = analysed_seq.aromaticity()
    return SEQ, LEN, Mw, pI, Aromaticity, Hydropathy

def submit_UniParc(UPI, SEQ, LEN, Mw, pI, Aromaticity, Hydropathy):
    submit_check = True
    submit_UniParc = ("INSERT INTO UniRef"
                   "(UniRef_ID, representative_accession, representative_name, protein_sequence, length, Mw, pI, aromaticity, hydropathy)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"
                   #Every time PhaLP is updated, the UniRef database is renewed. 
                   #The UPI remains static, but the representative might change, therefore this is updated at every run.
                   "ON DUPLICATE KEY UPDATE representative_accession=VALUES(representative_accession), representative_name=VALUES(representative_name), protein_sequence=VALUES(protein_sequence), length=VALUES(length), Mw=VALUES(Mw), pI=VALUES(pI), aromaticity=VALUES(aromaticity), hydropathy=VALUES(hydropathy)")
    try:
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(submit_UniParc,(UPI,None,None, SEQ, LEN, Mw, pI, Aromaticity, Hydropathy))
        cnx.commit()
        cursor.close()
        cnx.close()
    
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
            #print('\n')
            #print('%s already exists in UniParc' %UPI)

        else:
            print('\n')
            print('ERROR',err)
            submit_check = False
    return submit_check

def extract_domains(ACC, seq_rec):
    # get the response from the interpro API
    url = 'https://www.ebi.ac.uk/interpro/api/entry/all/protein/UniProt/%s/?page_size=100'%ACC 
    while True:
        try: 
            r = requests.get(url)
            if r.status_code == 408 or r.status_code ==404 or r.status_code == 500:
                print(r.status_code, 'wait 1 min', url)
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            if r.status_code == 200:
                interpro=True
                break
            elif r.status_code == 204:
                interpro = False
                break
            print('ERROR Interpro problem %s %s'%(r.status_code,r.text ))
        except (ConnectionError, ChunkedEncodingError, Timeout, HTTPError) as err:
            if err.code == 404:
                sleep(61)
                continue
            else:
                raise err
                #print('ERROR:\t%s\t'%ACC,err) 

    DOMAINS = {}
    # if the response is OK, read data.
    if interpro == True:
        js =r.json()
        if js['next']:
            raise ValueError('%s has more than 100 IPRs'%ACC)
        ipr_version = dir(r.headers['InterPro-Version'])       
        # get a dict with to get the interpro name from the interpo ACC
        ipr_list = []
        for element in js['results']:
            if element['metadata']['source_database'] != 'interpro':                
                accession = element['metadata']['accession']
                source_db = element['metadata']['source_database']
                source_db_version, source_db_release_date = canonical2version_date[source_db]
                name, short_name, integrated, found = fetch_InterPro(accession, source_db)
                if found:
                    if integrated:
                        ipr_ID = element['metadata']['integrated']
                        ipr_name, ipr_short_name = fetch_InterPro(ipr_ID, 'interpro')
                        ipr_list.append(ipr_ID)
                        ipr_version, ipr_release_date = canonical2version_date['interpro']
                    else:
                        ipr_ID = None
                        ipr_name = None
                        ipr_short_name = None
                        ipr_version = None
                        ipr_release_date = None
                    if accession not in DOMAINS.keys():
                            DOMAINS[accession]=[]
                    for match in element['proteins'][0]['entry_protein_locations']:                              
                        # for some matches, the match is not continious. This is ignored and the match is considered over the complete stretch. 
                        if match['fragments'][0]['dc-status']!='CONTINUOUS' or len(match['fragments']) > 1:
                            starts = []
                            ends = []
                            for fragment in match['fragments']:
                                starts.append(fragment['start'])
                                ends.append(fragment['end'])
                            start = min(starts)
                            end = max(ends)
                        else:
                            start = match['fragments'][0]['start']
                            end = match['fragments'][0]['end']

                        DOMAINS[accession].append({'source_db':source_db,
                                                    'source_db_version':source_db_version,
                                                    'source_db_release_date':source_db_release_date,
                                                    'name':name,
                                                    'short_name':short_name,
                                                    'ipr_ID':ipr_ID,
                                                    'ipr_name':ipr_name,
                                                    'ipr_short_name':ipr_short_name,
                                                    'ipr_version':ipr_version,
                                                    'ipr_release_date':ipr_release_date,
                                                    'start':start,
                                                    'end':end
                                                    })                
        for element in js['results']:
            if element['metadata']['source_database'] == 'interpro':
                if element['metadata']['accession'] not in ipr_list:
                    print('ERRORIPR: accessions %s has an IPR entry %s that is not yet in domains'%(ACC,element['metadata']['accession']) )
    else:
        print('%s not avaialbe in InterPro'%ACC)
    return DOMAINS
  
def fetch_InterPro(accession, source_db):
    url = "http://www.ebi.ac.uk/interpro/api/entry/%s/%s/"%(source_db, accession)
    while True:
        try:
            r = requests.get(url)
            if r.status_code == 408 or r.status_code ==404 or r.status_code == 500:
                js = r.json()
                print(js, type(js))
                print(js['Error'])
                if r.status_code ==404 and js['Error'] == "the level '%s' is not a valid DB member level"%accession:
                    print('Domain profile %s skipped, not in InterPro'%accession)
                    return None, None ,None , None
                print(r.status_code, 'wait 1 min', url)
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            if r.status_code == 200:
                break
            print('ERROR Interpro problem %s %s'%(r.status_code,r.text ))
           
        except (ConnectionError, ChunkedEncodingError, Timeout, HTTPError) as err:
            if err.code == 404:
                sleep(61)
                continue
            else:
                raise err
    js = r.json()
    if source_db == 'interpro':
        return js['metadata']['name']['name'], js['metadata']['name']['short']
    else:
        return js['metadata']['name']['name'], js['metadata']['name']['short'], js['metadata']['integrated'], True

def get_ipr_source_db(accession):
    url = "https://www.ebi.ac.uk/interpro/api/utils/accession/%s/"%( accession)
    while True:
        try:
            r = requests.get(url)
            if r.status_code == 408 or r.status_code == 500:
                print(r.status_code, 'wait 1 min', url)
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            if r.status_code == 200:
                js = r.json()
                return js['source_database']

            elif r.status_code == 404:
                return None

            print('ERROR Interpro problem %s %s'%(r.status_code,r.text ))
           
        except (ConnectionError, ChunkedEncodingError, Timeout, HTTPError) as err:
            if err.code == 404:
                sleep(61)
                continue
            else:
                raise err

def clean_link_UniParc_domains(UPI):
    clean_link_UniParc_domains = ("DELETE FROM link_UniRef_domains  WHERE link_UniRef_domains.UniRef_ID =%s; ")
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_link_UniParc_domains,(UPI,))
        cnx.commit()
        cursor.close()
        cnx.close()
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)
    return 
      
def submit_domains(DOMAINS, UPI):
    submit_domain = ("INSERT INTO domains"
                   "(domains_ID, domain_name, domain_name_short, source_database, InterPro_ID, InterPro_name, InterPro_name_short)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE domains_ID=VALUES(domains_ID), domain_name=VALUES(domain_name), domain_name_short=VALUES(domain_name_short), source_database=VALUES(source_database), InterPro_ID=VALUES(InterPro_ID), InterPro_name=VALUES(InterPro_name), InterPro_name_short=VALUES(InterPro_name_short)")    
    submit_link_domain = ("INSERT INTO link_UniRef_domains"
                   "(domains_ID, UniRef_ID, start, end, InterPro_version, InterPro_release_date)"
                   "VALUES (%s, %s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE domains_ID=VALUES(domains_ID), UniRef_ID=VALUES(UniRef_ID), start=VALUES(start), end=VALUES(end), InterPro_version=VALUES(InterPro_version), InterPro_release_date=VALUES(InterPro_release_date)")
    submit_check = True

    for dom_acc, value in DOMAINS.items():
        try:
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(submit_domain,(dom_acc, value[0]['name'], value[0]['short_name'], canonical2db_name[value[0]['source_db']], value[0]['ipr_ID'], value[0]['ipr_name'], value[0]['ipr_short_name']))
            cnx.commit()
            cursor.close()
            cnx.close() 
        except mysql.connector.Error as err:
                
            if err.errno == 1062:
            else:
                print('\n')
                print('ERROR',err)
                submit_check = False
        for value2 in value:
            try:
                cnx = cnxpool.get_connection()
                cursor = cnx.cursor()
                cursor.execute(submit_link_domain,(dom_acc, UPI, int(value2['start']), int(value2['end']), value[0]['ipr_version'], value[0]['ipr_release_date']) )
                cnx.commit()
                cursor.close()
                cnx.close() 
            except mysql.connector.Error as err:
                        
                if err.errno == 1062:
                else:
                    print('\n')
                    print('ERROR',err)
                    submit_check = False
    return submit_check

def determine_type(UPI, doms, SEQ):
    ip_locations, ip_functions, ip_names, ip_GO_ids, ip_GOs, ip_ECs = extract_identical_protein_annotations(UPI)
    TYPE = None 
    TYPE_EVIDENCE = None
    match = None        
    if ip_locations:
        for loc in ip_locations:
            if  loc == 'Host cell inner membrane' or loc == 'Host cell wall' or loc == 'Secreted':
                TYPE='endolysin'
                TYPE_EVIDENCE = 'UniProt location annotation'
                match = loc
                break
            elif loc == 'Virion' or loc == 'Virion membrane':
                TYPE_EVIDENCE = 'UniProt location annotation'
                match = loc
                break           
    
    if not TYPE and ip_functions:
        terms_endolysin = ['Essential for lysis of bacterial cell wall, by showing cell wall hydrolyzing activity.',
                            'The bacteriophage lysin is a powerful lytic agent',
                            'participates in the liberation of progeny bacteriophage',
                            'Helps to release the mature phage particles from the cell wall by breaking down the peptidoglycan',
                            'Endolysin that degrades host peptidoglycans'                   
                            ]
        for fun in ip_functions:
            m = re.search(r'|'.join(terms_endolysin),fun)
            if m:
                TYPE='endolysin'
                TYPE_EVIDENCE = 'UniProt function annotation'
                match = fun
                break                     
                                     
    if not TYPE and ip_GO_ids:
        for GO_id, GO in zip(ip_GO_ids, ip_GOs):
            if GO_id in ['GO:0044659', 'GO:0019076', 'GO:0051715','GO:0001897']:
                TYPE='endolysin'
                TYPE_EVIDENCE = 'GO annotation'
                match = GO
                break
            elif GO_id in ['GO:0098025', 'GO:0019062','GO:0007155','GO:0019012','GO:0055036','GO:0039641','GO:0098003','GO:0098015','GO:0039678','GO:0019867','GO:0099001','GO:0019069','GO:0043493','GO:0019073','GO:0098023','GO:0098004','GO:0044694']:
                TYPE='VAL'
                TYPE_EVIDENCE = 'GO annotation'
                match = GO
                break
        
    if not TYPE and doms:
        for dom in doms:
            if dom in ['cd16901','cd06415']:
                TYPE='endolysin'
                TYPE_EVIDENCE = 'InterPro domain annotation'
                match = dom
                break
            elif dom in ['PF06715','SSF69349','SSF69255','G3DSA:2.40.50.260','PF06714','PF10145','TIGR02675','TIGR01760','PF06605','TIGR01665','SSF69279','PF06791','PF16838']:
                TYPE='VAL'
                TYPE_EVIDENCE = 'InterPro domain annotation'
                match = dom

    if not TYPE and ip_names:
        terms_VAPGH = ['tail','virion','capsid','core','head','injection','plate','tapemeasure','tape measure']
        terms_endolysin = ['endolysin','endolysiin','Lysin A']#,'lysin']#ply|lysA|lytic|lysis]
        for name in ip_names:
            if re.search(r'|'.join(terms_endolysin),name, re.IGNORECASE ) and not re.search(r'putative',name, re.IGNORECASE ):
                TYPE='endolysin'
                TYPE_EVIDENCE = 'Protein name annotation'
                match = name
                break
            elif re.search(r'|'.join(terms_VAPGH),name, re.IGNORECASE ) and not re.search(r'putative',name, re.IGNORECASE ):
                TYPE='VAL'
                TYPE_EVIDENCE = 'Protein name annotation'
                match = name
                break
    # if the protein was used to train the final RF classification model, use the predicted value from the 10-fold CV model where the protein was witheld from the training dataset. 
    if UPI in list(KfoldCVpredictions.index):
        # if the current proteion annotations doesn't allow it to be classified anymore, use ML prediction as evidence.
        if not TYPE: 
            TYPE_EVIDENCE = 'ML prediction'
        row = KfoldCVpredictions.loc[UPI]
        if row['true_val'] != row['pred_val']:
            TYPE = row['true_val']
            TYPE_PROBABILITY = round(1.0 - row['pred_prob'], 2)*100
        else:
            TYPE = row['true_val']
            TYPE_PROBABILITY = row['pred_prob']*100
    # if not, but it's type was classified based on annotation, keep the type and get the predicted probability according to the RF classification model
    elif TYPE:
        TYPE_predict, TYPE_PROBABILITY_predict = predict_type(UPI, SEQ)
        if TYPE == TYPE_predict:
            TYPE_PROBABILITY = TYPE_PROBABILITY_predict*100
        else:
            TYPE_PROBABILITY = round(1.0 - TYPE_PROBABILITY_predict, 2)*100
    else:
        TYPE, TYPE_PROBABILITY = predict_type(UPI, SEQ)
        TYPE_PROBABILITY = TYPE_PROBABILITY*100
        TYPE_EVIDENCE = 'ML prediction'
    
    return TYPE, TYPE_EVIDENCE, TYPE_PROBABILITY

def predict_type(UPI, SEQ):
    if UPI in SeqVecEmbeddings_index:
        embedding = pandas.read_csv('picklejar/SeqVecEmbeddings.csv', index_col = 0).loc[UPI].to_frame().T 
    else:
        embedding = embedder.embed(SEQ)
        embedding = embedder.reduce_per_protein(embedding)
        embedding = pandas.DataFrame(embedding, columns = [UPI]).T
        with open('picklejar/SeqVecEmbeddings.csv', 'a') as f:
            embedding.to_csv(f, header=False)
            SeqVecEmbeddings_index.append(UPI)
    TYPE = clf.predict(embedding)[0]
    TYPE_PROBABILITY = max(clf.predict_proba(embedding)[0])
    return TYPE, TYPE_PROBABILITY

def extract_identical_protein_annotations(upi):
    uniprot_url = 'http://www.uniprot.org/uniprot/'
    search_params = {
        'query':'uniparc:(%s)'%upi,
        'format': 'tab',
        'columns': 'comment(SUBCELLULAR LOCATION),comment(FUNCTION),protein names,go-id,go,id',
        }
    while True:
        try:
            response = requests.get(uniprot_url, search_params)
            df = pandas.read_csv(StringIO(response.text), sep="\t")
            break
        except (ConnectionError, pandas.errors.ParserError, pandas.errors.EmptyDataError) as err:
            print('ERROR:\t%s\t'%upi,err) 
    
    all_loc = []
    all_func = []
    all_nam = []
    all_EC = []
    all_GO_id = []
    all_GO = []
    for loc in filter(lambda x : str(x)!='nan',df['Subcellular location [CC]'].tolist()):
        all_loc.extend([j.split(' {')[0].split(': ')[-1] for i in loc.split('SUBCELLULAR LOCATION: ')[1:] for j in i.split('.')[0].split('; ') ])
    for func in filter(lambda x : str(x)!='nan', df['Function [CC]'].tolist()):
        all_func.extend([i.split(' {')[0] for i in func[10:].split('; FUNCTION: ')])
    for GO_id in filter(lambda x : str(x)!='nan',df['Gene ontology IDs'].tolist()):
        all_GO_id.extend(GO_id.split('; '))
    for GO in filter(lambda x : str(x)!='nan',df['Gene ontology (GO)'].tolist()):
        all_GO.extend(GO.split('; '))
    for nam in [str(i) for i in filter(lambda x : str(x)!='nan',df['Protein names'].tolist())]:
        all_EC.extend(re.findall(r'\((EC.*?)\)',nam))
        all_nam.append(nam.split(' (')[0])
        all_nam.extend(list(filter(lambda x: x not in all_EC, re.findall(r'\((.*?)\)', nam))))
    locations = None
    functions = None
    names = None
    ECs = None
    GOs = None
    GO_ids = None
    
    # remove all duplicate entries in the list and order on frequency (most frequent first)
    if all_loc:
        locations = [i[0] for i in collections.Counter(all_loc).most_common()]
    if all_func:
        functions = [i[0] for i in collections.Counter(all_func).most_common()]
    if all_EC:
        ECs = [i[0] for i in collections.Counter(all_EC).most_common()]
    if all_nam:
        names = [i[0] for i in collections.Counter(all_nam).most_common()]
    if all_GO_id:
        GO_ids = [i[0] for i in collections.Counter(all_GO_id).most_common()]    
    if all_GO:
        GOs = [i[0] for i in collections.Counter(all_GO).most_common()]   
    return locations, functions, names, GO_ids, GOs, ECs

def submit_UniProt(ACC,PROTEIN_NAME, TYPE, TYPE_EVIDENCE, TYPE_PROBABILITY, UPI, PHAGE_ID,DATE_CREATED, DATE_LAST_UPDATED):
    submit_check = True
    if TYPE == '':
        TYPE = None
    submit_UniProt = ("INSERT INTO UniProt"
                   "(UniProt_ID, name, type, type_evidence, type_probability, UniRef_ID, phages_ID, date_created, date_last_updated)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE UniProt_ID=VALUES(UniProt_ID), name=VALUES(name), type=VALUES(type),type_evidence=VALUES(type_evidence), type_probability=VALUES(type_probability), UniRef_ID=VALUES(UniRef_ID), phages_ID=VALUES(phages_ID), date_created =VALUES(date_created), date_last_updated =VALUES(date_last_updated) ")
    try:
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(submit_UniProt,(ACC, PROTEIN_NAME, TYPE, TYPE_EVIDENCE, str(TYPE_PROBABILITY), UPI, PHAGE_ID, DATE_CREATED, DATE_LAST_UPDATED))## 
        cnx.commit()
        cursor.close()
        cnx.close()     
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
            #print('\n')
            #print('%s already exists in UniProt' %ACC)
        else:
            print('\n')
            print('ERROR',err)
            submit_check = False
    return submit_check

def extract_references_experimental_evidence(swiss_rec,ACC):   
    pmid_list=[]
    for comment in swiss_rec.comments:
        if re.search(r'ECO:0000269', comment):
            if re.match(r'FUNCTION|CATALYTIC ACTIVITY',comment) :
                if re.search('exolysin|endolysin|lysozyme', comment, re.IGNORECASE) or (re.search('hydrolyzes|hydrolysis|cleavage', comment, re.IGNORECASE) and re.search('peptidoglycan|cell-wall glycopeptides', comment, re.IGNORECASE)):
                                
                    for evidences in re.findall(r'{(.*?)}', comment):
                        for evidence in evidences.split(', '):
                            m = re.match('ECO:0000269\|(.*)',evidence)
                            if m:
                                for ref in m.groups():
                                    pmid_list.append(ref.replace('PubMed:',''))
                        
    GO_exp_list = ['GO:0003796','GO:0051672','GO:0033922','GO:0085027','GO:0008933','GO:0044659','GO:0008745','GO:0051715','GO:0098932','GO:0061783','GO:0061784','GO:0061785','GO:0009253','GO:0061786','GO:0061787']
    GO_tsv_downloaded = False
    for xref in swiss_rec.cross_references:
        if xref[0] == 'GO' and xref[-1][:3] in ['EXP','IDA', 'IPI','IMP','IGI','IEP']:
            GO_ID =xref[1] 
            if GO_ID in GO_exp_list:
                
                if GO_tsv_downloaded is False:
                    url_xml =   'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=%s'%ACC  
                    while True:
                        try:
                            
                            r = requests.get(url_xml,headers={ "Accept" : "text/tsv"})
                            with open('downloads/temporary_GO.tsv', 'wt') as f:
                                f.write(r.text)
                            break
                        except (ConnectionError, ChunkedEncodingError) as err:
                            print('ERROR:\t%s\t'%ACC,err) 
                    
                    GO_tsv_downloaded = True
                with open('downloads/temporary_GO.tsv', 'rt') as GO_handle:
                    for line in csv.reader(GO_handle, delimiter='\t'):
                        if line[4]==GO_ID and line[7] in ['EXP','IDA', 'IPI','IMP','IGI','IEP']:
                            pmid_list.append(line[8].replace('PMID:',''))
    pmid_list = list(set(pmid_list))
    ref_dict={}
    if pmid_list:
        while True:
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode="text")
                break
            except (URLError, IncompleteRead, ValueError, RuntimeError, TimeoutError) as err:
                print(err)
       
        for rec in Medline.parse(handle):
            authors = ', '.join(rec['AU'])
            title = rec['TI']
            pmid = rec['PMID']
            date_publ = rec['DP']
            ref_dict[pmid]=[title, authors, date_publ]
    return((ref_dict))

def clean_UniProt_child(ACC,child_table):
    if child_table in ['experimental_evidence','gene_ontologies','link_EC_UniProt','tertiary_structures','CDSs']:
        clean_child = ("DELETE FROM {}  WHERE {}.UniProt_ID =%s; ".format(child_table, child_table))
        try:
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(clean_child,( ACC,))
            cnx.commit()
            cursor.close()
            cnx.close()  
        except mysql.connector.Error as err:
            if err.errno == 1062:
                None
            else:
                print('\n')
                print('ERROR',err)
        return
    else:
        raise ValueError('Table is no child of Uniprot table. Cannot delete entries.')

def submit_experimental_evidence(ACC,REF_EXP_EVID):
    submit_check = True
    submit_exp = ("INSERT INTO experimental_evidence"
                   "(UniProt_ID, pmid, title, authors, date_of_publication)"
                   "VALUES (%s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE UniProt_ID=VALUES(UniProt_ID), pmid=VALUES(pmid), title=VALUES(title), authors=VALUES(authors), date_of_publication=VALUES(date_of_publication)")
    for pmid, value in REF_EXP_EVID.items():
        try: 
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(submit_exp,(ACC, pmid, value[0], value[1], value[2]))
            cnx.commit()
            cursor.close()
            cnx.close()      
        except mysql.connector.Error as err:
            
            if err.errno == 1062:
                None
            else:
                print('\n')
                print('ERROR',err)
                submit_ckeck = False

def extract_GOs(swiss_rec):
    GOs={}
    category_dict = {'P':'Biological process','C':'Cellular component','F':'Molecular function'}
    evidence_dict = {   'EXP': 'Inferred from Experiment',
                        'IBA': 'Inferred from Biological Aspect of Ancestor',
                        'IC': 'Inferred by Curator',
                        'IDA': 'Inferred from Direct Assay',
                        'HDA': 'High-throughput Direct Assay',
                        'IEA': 'Inferred from Electronic Annotation',
                        'IEP': 'Inferred from Expression Pattern',
                        'HEP': 'High-throughput Expression Pattern',
                        'IGC': 'Inferred from Genomic Context',
                        'IGI': 'Inferred from Genetic Interaction',
                        'HGI': 'High-throughput Genetic Interaction',
                        'IMP': 'Inferred from Mutant Phenotype',
                        'HMP':' High-throughput Mutant Phenotype',
                        'IPI': 'Inferred from Physical Interaction',
                        'ISA': 'Inferred from Sequence Alignment',
                        'ISM': 'Inferred from Sequence Model',
                        'ISO': 'Inferred from Sequence Orthology',
                        'ISS': 'Inferred from Sequence or Structural Similarity',
                        'NAS': 'Non-traceable Author Statement',
                        'TAS': 'Traceable Author Statement'
                        }
    for xref in swiss_rec.cross_references:
        if xref[0] == 'GO':
            GO_ID = xref[1]
            GO_name = xref[2][2:]
            GO_category = category_dict[xref[2][0]]
            GO_evidence = evidence_dict[xref[3].split(':')[0]]
            GO_source = xref[3].split(':')[1]
            if GO_ID in GOs.keys():
                print('DOUBLE')
            GOs[GO_ID]=[GO_name,GO_category,GO_evidence,GO_source]
            
    return GOs

def submit_GOs(ACC,GOs):
    submit_check = True
    submit_GO = ("INSERT IGNORE INTO gene_ontologies"
                   "(UniProt_ID, gene_ontologies_ID, category, name, evidence, source)"
                   "VALUES (%s, %s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE UniProt_ID=VALUES(UniProt_ID), gene_ontologies_ID=VALUES(gene_ontologies_ID), category=VALUES(category), name=VALUES(name), evidence=VALUES(evidence), source=VALUES(source)")
    for GO_ID, value in GOs.items():
        try:
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(submit_GO,(ACC, GO_ID, value[1], value[0],value[2], value[3]))
            cnx.commit()
            cursor.close()
            cnx.close()    
        except mysql.connector.Error as err:
            if err.errno == 1062:
                None
            else:
                print('\n')
                print('ERROR',err)
                submit_check = False
    return  submit_check

def extract_EC(description):
    ECs = {}
    for EC in re.findall(r'EC=(.*?)(?: {(.*?)})?;', description):
        EC_ID = EC[0]
        ECs[EC_ID] = {'evidence':[]}
        if EC[1]:
            for evidence in EC[1].split(', '):
                ECs[EC_ID]['evidence'].append(evidence.split('|'))
    return ECs

def submit_EC(ECs, ACC):
    submit_check = True
    submit_EC = ("INSERT IGNORE INTO EC"
                   "(EC_ID, entry_name, alternative_names, reaction_catalyzed, class_ID, class_name, subclass_ID, subclass_name, subsubclass_ID, subsubclass_name)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE EC_ID=VALUES(EC_ID), entry_name=VALUES(entry_name), alternative_names=VALUES(alternative_names), reaction_catalyzed=VALUES(reaction_catalyzed), class_ID=VALUES(class_ID), class_name=VALUES(class_name), subclass_ID=VALUES(subclass_ID), subclass_name=VALUES(subclass_name), subsubclass_ID=VALUES(subsubclass_ID), subsubclass_name=VALUES(subsubclass_name)")
    submit_evidence = ("INSERT INTO link_EC_UniProt"
                   "(EC_ID, UniProt_ID, source, evidence, ECO)"
                   "VALUES (%s, %s, %s, %s, %s)"
                   "ON DUPLICATE KEY UPDATE EC_ID=VALUES(EC_ID), UniProt_ID=VALUES(UniProt_ID), source=VALUES(source), evidence=VALUES(evidence), ECO=VALUES(ECO)")

    with shelve.open("picklejar/EC_shelf.db", flag='r') as EC_shelf:
        for EC_ID, value in ECs.items():

            for key, value2 in EC_shelf[EC_ID].items():
                if value2 == '':
                    EC_shelf[EC_ID][key]=None
            try:
                cnx = cnxpool.get_connection()
                cursor = cnx.cursor()
                cursor.execute(submit_EC,(EC_ID, EC_shelf[EC_ID]['name'], EC_shelf[EC_ID]['alt_names'], EC_shelf[EC_ID]['reaction'], EC_shelf[EC_ID]['class_ID'], EC_shelf[EC_ID]['class_name'], EC_shelf[EC_ID]['subclass_ID'], EC_shelf[EC_ID]['subclass_name'], EC_shelf[EC_ID]['subsubclass_ID'], EC_shelf[EC_ID]['subsubclass_name']))
                cnx.commit()
                cursor.close()
                cnx.close()       
            except mysql.connector.Error as err:
                    
                if err.errno == 1062:
                    None
                else:
                    print('\n')
                    print(EC_ID, ACC )
                    print('ERROR',err)
                    submit_check = False
            for ECO, source in value['evidence']:
                ECO_dict = {'ECO:0000269':'experimental evidence used in manual assertion',
                            'ECO:0000303':'non-traceable author statement used in manual assertion',
                            'ECO:0000305':'curator inference used in manual assertion',
                            'ECO:0000250':'sequence similarity evidence used in manual assertion',
                            'ECO:0000255':'match to sequence model evidence used in manual assertion',
                            'ECO:0000256':'match to sequence model evidence used in automatic assertion',
                            'ECO:0000259':'match to InterPro member signature evidence used in automatic assertion',
                            'ECO:0000312':'imported information used in manual assertion',
                            'ECO:0000313':'imported information used in automatic assertion',
                            'ECO:0000244':'combinatorial evidence used in manual assertion',
                            'ECO:0000213':'combinatorial evidence used in automatic assertion'
                    }
                try:
                    cnx = cnxpool.get_connection()
                    cursor = cnx.cursor()
                    cursor.execute(submit_evidence,(EC_ID, ACC, source, ECO_dict[ECO], ECO ))
                    cnx.commit()
                    cursor.close()
                    cnx.close()   
                except mysql.connector.Error as err:
                        
                    if err.errno == 1062:
                        None
                    else:
                        print('\n')
                        print(EC_ID, ACC, source, ECO_dict[ECO], ECO )
                        print('ERROR',err)
                        submit_check = False
    return submit_check

def extract_3D(ACC, swiss_rec):
    PDBs = {}
    for xref in swiss_rec.cross_references:
        if xref[0] == 'PDB':#=='RefSeq' or xref[0]=='EMBL':
            PDB_ID, method, resolution, chain_position =  xref[1:]
            PDBs[PDB_ID]={'method':method,'resolution':resolution,'chain_position':chain_position}
    return PDBs
        
def submit_3D_structures(ACC, PDBs):
    for PDB_ID in PDBs.keys():
        if PDBs[PDB_ID]['resolution'] == '-':
            PDBs[PDB_ID]['resolution'] = None
        submit_3D = ("INSERT INTO tertiary_structures"
                       "(UniProt_ID, PDB_ID, method, resolution, chain_position)"
                       "VALUES (%s, %s, %s, %s, %s)"
                       "ON DUPLICATE KEY UPDATE UniProt_ID=VALUES(UniProt_ID), PDB_ID=VALUES(PDB_ID), method=VALUES(method), resolution=VALUES(resolution), chain_position=VALUES(chain_position)")
        try: 
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(submit_3D,(ACC,PDB_ID, PDBs[PDB_ID]['method'], PDBs[PDB_ID]['resolution'], PDBs[PDB_ID]['chain_position']))
            cnx.commit()
            cursor.close()
            cnx.close()     
        except mysql.connector.Error as err:
            
            if err.errno == 1062:
                None
            else:
                print('\n')
                print('ERROR',err)

def extract_CDSs(ACC, swiss_rec):
    CDSs = {}
    sequence_conflicts = {}
    for feat in swiss_rec.features:
        if feat.type=='CONFLICT':
            print(feat)##TEMP##
            ##TODO## change feat indexes to names
            m = re.search(r'.*?;\s(.*?)\)', feat.qualifiers['note'])
            if m and m.groups()[0] != 'AA sequence':
                genpept_acc=m.groups()[0]
                residue_location = str(int(feat.location.split(':')[0][1:])+1)
                if feat.qualifiers['note'].split(' (')[0] == 'Missing':
                    genpept_residue='-'
                    uniprot_residue=swiss_rec.sequence[residue_location] 
                else:
                    uniprot_residue, genpept_residue= feat.qualifiers['note'].split(' (')[0].split(' -> ')
                if genpept_acc in sequence_conflicts.keys():
                    sequence_conflicts[genpept_acc].append({'res_num':residue_location, 'genpept_res':genpept_residue, 'uniprot_res':uniprot_residue})
                else:
                    sequence_conflicts[genpept_acc]=[{'res_num':residue_location, 'genpept_res':genpept_residue, 'uniprot_res':uniprot_residue}]
    
    for xref in swiss_rec.cross_references:
        if xref[0]=='EMBL' and xref[3] == '-':
            gp_rec = fetch_entrez(xref[2].split('\.'[0]),'protein')
                        
            gp_seq=deepcopy(gp_rec.seq)
            up_seq=deepcopy(swiss_rec.sequence)
                    
            gp_acc=gp_rec.id
            if gp_acc.split('.')[0] in sequence_conflicts.keys():
                for sequence_conflict in sequence_conflicts[gp_acc.split('.')[0]]:
                    assert up_seq[sequence_conflict['res_num']-1] == sequence_conflict['uniprot_res']
                    up_seq = up_seq[:sequence_conflict['res_num']-1]+sequence_conflict['genpept_res']+up_seq[sequence_conflict['res_num']:]

            if str(up_seq[1:]).replace('-','') == str(gp_seq)[1:]:
                #print('\tCORRECT')
                coded_by =None
                for feature in gp_rec.features:
                    if feature.type == 'CDS' and "coded_by" in feature.qualifiers:
                        coded_by = feature.qualifiers['coded_by'][0]
                else:
                    if '<' not in coded_by and '>' not in coded_by:
                        CDS , gb_acc, start, end, strand = fetch_nucleotide_from_coded_by_string(coded_by) 
                        if gp_acc in CDSs.keys():
                            print('DUPLICATE CDS')
                        CDSs[gp_acc]={'gb_acc':gb_acc,'CDS':CDS,'start':start,'end':end,'strand':strand,'sequence_conflicts':[]}
                        if gp_acc.split('.')[0] in sequence_conflicts.keys():
                            for sequence_conflict in sequence_conflicts[gp_acc.split('.')[0]]:
                                CDSs[gp_acc]['sequence_conflicts'].append(sequence_conflict)
    return CDSs

def fetch_entrez(acc, db, start=None, end=None):
    while True:
        try:
            handle = Entrez.efetch(db=db, id=acc, seq_start=start, seq_stop=end, rettype="gb", retmode="text")
            break
        except (URLError, IncompleteRead, ValueError, RuntimeError, TimeoutError) as err: #OSError: [Errno 65] No route to host
            print(err)
    seq_rec=SeqIO.read(handle,'gb')
      return seq_rec

def fetch_nucleotide_from_coded_by_string(coded_by, complement = False, join = False):      
    if coded_by.startswith("complement("):
        assert coded_by.endswith(")")
        return fetch_nucleotide_from_coded_by_string(coded_by[11:-1], complement=True, join=join)
    if coded_by.startswith("join(") :
        assert coded_by.endswith(")")
        return fetch_nucleotide_from_coded_by_string(coded_by[5:-1], complement = complement, join = True) 
    if "(" in coded_by or ")" in coded_by or ((coded_by.count(":") != 1 or coded_by.count("..") != 1) and join == False):
        raise ValueError("Don't understand %s" % repr(coded_by))
    if join  == False:
        name, loc = coded_by.split(":")
        start, end = [int(x.lstrip("<>")) for x in loc.split("..")]
        if complement == False:
            strand = 1
            CDS = fetch_fasta(name.strip(), 'nucleotide', start, end)
        else:
            strand = -1
            CDS = fetch_fasta(name.strip(), 'nucleotide', start, end).reverse_complement()
        return CDS, name, start, end, strand
    else:
        names = []
        starts = []
        ends = []
        seqs = []
        for part in coded_by.split(","):
            part = part.strip(' ')
            name_part, loc_part = part.split(':')
            start_part, end_part = [str(x.lstrip("<>")) for x in loc_part.split("..")]
            if complement == False:
                seq_part = str(fetch_nucleotide_from_coded_by_string(part)[0])
            else:
                seq_part = str(fetch_nucleotide_from_coded_by_string(part)[0].reverse_complement())
            names.append(name_part)
            starts.append(start_part)
            ends.append(end_part)
            seqs.append(seq_part)   
        name=";".join(names)
        end=";".join(ends)
        start=";".join(starts)
        if complement == True:
            seqs = seqs[::-1]       #reverse order of parts
        seq = SeqRecord(Seq("".join(seqs)))
        if complement == False:
            strand = 1
        else:
            strand = -1  
        return seq.seq , name, start, end, strand    

def fetch_fasta(acc, db, start=None, end=None):
    while True:
        try:
            handle = Entrez.efetch(db=db, id=acc, seq_start=start, seq_stop=end, rettype="fasta", retmode="text")
            break
        except (URLError, IncompleteRead, ValueError, RuntimeError, TimeoutError) as err:
            print(err)
    seq_rec=SeqIO.read(handle,'fasta')
    
    return seq_rec.seq

def submit_CDSs(CDSs, ACC):
    submit_check = True
    submit_CDS = ("INSERT INTO CDSs"
                   "(GenBank_accession_id, CDSs_ID, CDS, start, end, strand, UniProt_ID)"
                   "VALUES (%s, %s, %s, %s, %s, %s, %s)"
                    "ON DUPLICATE KEY UPDATE GenBank_accession_id=VALUES(GenBank_accession_id), CDSs_ID=VALUES(CDSs_ID), CDS=VALUES(CDS), start=VALUES(start), end=VALUES(end), strand=VALUES(strand), UniProt_ID=VALUES(UniProt_ID)")
    submit_sequence_conflict = ("INSERT INTO sequence_conflicts"
                                "(UniProt_ID, CDSs_ID, residue_number, UniProt_residue, GenPept_residue)"
                                "VALUES (%s, %s, %s, %s, %s)"
                                "ON DUPLICATE KEY UPDATE UniProt_ID=VALUES(UniProt_ID), CDSs_ID=VALUES(CDSs_ID), residue_number=VALUES(residue_number), UniProt_residue=VALUES(UniProt_residue), GenPept_residue=VALUES(GenPept_residue)")

    for gp_acc, value in CDSs.items():
        try:
            cnx = cnxpool.get_connection()
            cursor = cnx.cursor()
            cursor.execute(submit_CDS,(value['gb_acc'], gp_acc, str(value['CDS']), value['start'], value['end'], value['strand'], ACC))
            cnx.commit()
            cursor.close()
            cnx.close()
        except mysql.connector.Error as err:
            if err.errno == 1062:
                print('\n')
                print('%s already exists in CDS' %gp_acc)
            else:
                print('\n')
                print('ERROR',err)
                submit_check = False
        for sequence_conflict in value['sequence_conflicts']:
            try:
                cnx = cnxpool.get_connection()
                cursor = cnx.cursor()
                cursor.execute(submit_sequence_conflict,(ACC, gp_acc, sequence_conflict['res_num'], sequence_conflict['uniprot_res'], sequence_conflict['genpept_res'] ))
                cnx.commit()
                cursor.close()
                cnx.close()
            except mysql.connector.Error as err:
                if err.errno == 1062:
                    None
                else:
                    print('\n')
                    print('ERROR',err)
                    submit_check = False
    return submit_check

def delete_UniProt_entry(ACC):
    delete_UniProt = ("DELETE FROM UniProt "
                        "WHERE UniProt_ID =%s; ")
    
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(delete_UniProt,(ACC,))
        cnx.commit()
        cursor.close()
        cnx.close() 
        print(ACC, 'deleted from database.',cursor.rowcount,'rows deleted')
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)
    return 

def clean_up_parent_tables():
    clean_phages = ("DELETE FROM phages WHERE NOT EXISTS (SELECT phages_ID FROM UniProt WHERE phages.phages_ID =UniProt.phages_ID);")
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_phages)
        cnx.commit()
        cursor.close()
        cnx.close() 
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)

    
    clean_hosts = ("DELETE FROM hosts WHERE NOT EXISTS (SELECT hosts_ID FROM link_phage_host WHERE hosts.hosts_ID =link_phage_host.hosts_ID);")
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_hosts)
        cnx.commit()
        cursor.close()
        cnx.close() 
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)


    clean_UniRef = ("DELETE FROM UniRef "
                    "WHERE NOT EXISTS ("
                        'SELECT UniRef_ID FROM UniProt '
                        'WHERE UniRef.UniRef_ID = UniProt.UniRef_ID);')
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_UniRef)
        cnx.commit()
        cursor.close()
        cnx.close() 
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)

    clean_domains = ("DELETE FROM domains "
                    "WHERE NOT EXISTS ("
                        'SELECT domains_ID FROM link_UniRef_domains '
                        'WHERE link_UniRef_domains.domains_ID = domains.domains_ID);')
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_domains)
        cnx.commit()
        cursor.close()
        cnx.close() 
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)
     

    clean_EC = ("DELETE  FROM EC "
                "WHERE NOT EXISTS ("
                    'SELECT EC_ID FROM link_EC_UniProt '
                    'WHERE link_EC_UniProt.EC_ID = EC.EC_ID);')
    try: 
        cnx = cnxpool.get_connection()
        cursor = cnx.cursor()
        cursor.execute(clean_EC)
        cnx.commit()
        cursor.close()
        cnx.close() 
    except mysql.connector.Error as err:
        
        if err.errno == 1062:
            None
        else:
            print('\n')
            print('ERROR',err)
    return 
