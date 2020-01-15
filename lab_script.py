# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 17:01:52 2020

@author: carina afonso e laura duro
"""

from Bio import SeqIO
from Bio import Entrez as E
from Bio import Medline
import re


# ____________________________ Literature analysis ___________________________

E.email = 'pg40959@alunos.uminho.pt'

## Pubmed articles - general

# __________ Genes __________
names_genes = ['KLHL22', 'SCARF2', 'ZNF74']

dic_art = { names_genes[0]: {}, names_genes[1]: {}, names_genes[2]: {}}
for gene in dic_art.keys():
    term = 'Human[Orgn] AND ' + gene + '[Gene]'
    handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
    list_id_art =  E.read(handle)['IdList']
    file_name = 'Articles_' + gene + '.txt'
    file = open(file_name, 'w')
    handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
    records = list(Medline.parse(handle))
    for record in records: #for each pubmed id, store following parameters
        id_pubmed = str(record.get('PMID'))
        title = str(record.get('TI'))
        authors = str(record.get('AU'))
        abstract = str(record.get('AB'))
        source = str(record.get('SO'))
        file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
    dic_art[gene] = list_id_art
    file.close()

# __________ Cancer __________
## Adenocarcinoma


key = 'Adenocarcinoma'
dic_art[key] = {}
term = 'Human[Orgn] AND Adenocarcinoma'
handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
list_id_art =  E.read(handle)['IdList']
file_name = key + '.txt'
file = open(file_name, 'w')
handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
records = list(Medline.parse(handle))
for record in records: #for each pubmed id, store following parameters
    id_pubmed = str(record.get('PMID'))
    title = str(record.get('TI'))
    authors = str(record.get('AU'))
    abstract = str(record.get('AB'))
    source = str(record.get('SO'))
    file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
file.close()
dic_art[key] = list_id_art

## Adenocarcinoma caracterization
key = 'Ovary adenocarcinoma characterization'
dic_art[key] = {}
term = 'Human[Orgn] AND ' + key
handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
list_id_art =  E.read(handle)['IdList']
file_name = key + '.txt'
file = open(file_name, 'w')
handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
records = list(Medline.parse(handle))
for record in records: #for each pubmed id, store following parameters
    id_pubmed = str(record.get('PMID'))
    title = str(record.get('TI'))
    authors = str(record.get('AU'))
    abstract = str(record.get('AB'))
    source = str(record.get('SO'))
    file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
file.close()
dic_art[key] = list_id_art

## Ovary adenocarcinoma
key = 'Ovary Adenocarcinoma'
dic_art[key] = {}
term = 'Human[Orgn] AND Ovary adenocarcinoma'
handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
list_id_art =  E.read(handle)['IdList']
file_name = key + '.txt'
file = open(file_name, 'w')
handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
records = list(Medline.parse(handle))
for record in records: #for each pubmed id, store following parameters
    id_pubmed = str(record.get('PMID'))
    title = str(record.get('TI'))
    authors = str(record.get('AU'))
    abstract = str(record.get('AB'))
    source = str(record.get('SO'))
    file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
file.close()
dic_art[key] = list_id_art

## celular line OV90
key = 'OV90 characterization'
dic_art[key] = {}
term = 'Human[Orgn] AND ' + key
handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
list_id_art =  E.read(handle)['IdList']
file_name = key + '.txt'
file = open(file_name, 'w')
handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
records = list(Medline.parse(handle))
for record in records: #for each pubmed id, store following parameters
    id_pubmed = str(record.get('PMID'))
    title = str(record.get('TI'))
    authors = str(record.get('AU'))
    abstract = str(record.get('AB'))
    source = str(record.get('SO'))
    file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
file.close()
dic_art[key] = list_id_art


# __________ Cancer/Line + Genes __________
## Adenocarcionama and genes
for i in range(3):
    key = 'Adenocarcinoma AND ' + names_genes[i]
    dic_art[key] = {}
    term = 'Human[Orgn] AND ' + key + '[Gene]'
    handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
    list_id_art =  E.read(handle)['IdList']
    file_name = key.replace(' ', '_') + '.txt'
    file = open(file_name, 'w')
    handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
    records = list(Medline.parse(handle))
    for record in records: #for each pubmed id, store following parameters
        id_pubmed = str(record.get('PMID'))
        title = str(record.get('TI'))
        authors = str(record.get('AU'))
        abstract = str(record.get('AB'))
        source = str(record.get('SO'))
        file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
    file.close()
    dic_art[key] = list_id_art

##Ovary adenocarcionama and genes
for i in range(3):
    key = 'Ovary adenocarcinoma AND ' + names_genes[i]
    dic_art[key] = {}
    term = 'Human[Orgn] AND ' + key + '[Gene]'
    handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
    list_id_art =  E.read(handle)['IdList']
    file_name = key.replace(' ', '_') + '.txt'
    file = open(file_name, 'w')
    handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
    records = list(Medline.parse(handle))
    for record in records: #for each pubmed id, store following parameters
        id_pubmed = str(record.get('PMID'))
        title = str(record.get('TI'))
        authors = str(record.get('AU'))
        abstract = str(record.get('AB'))
        source = str(record.get('SO'))
        file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
    file.close()
    dic_art[key] = list_id_art

## line + genes
for i in range(3):
    key = 'OV90 AND ' + names_genes[i]
    dic_art[key] = {}
    term = 'Human[Orgn] AND ' + key + '[Gene]'
    handle = E.esearch(db = 'pubmed', term = term) #search in pubmed
    list_id_art =  E.read(handle)['IdList']
    file_name = key.replace(' ', '_') + '.txt'
    file = open(file_name, 'w')
    handle = E.efetch(db = 'pubmed', id = list_id_art, rettype = 'Medline', retmode = 'text')             
    records = list(Medline.parse(handle))
    for record in records: #for each pubmed id, store following parameters
        id_pubmed = str(record.get('PMID'))
        title = str(record.get('TI'))
        authors = str(record.get('AU'))
        abstract = str(record.get('AB'))
        source = str(record.get('SO'))
        file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
    file.close()
    dic_art[key] = list_id_art

# # ______________ Retrieve sequences and pubmed related articles _____________    

dic_genes = {names_genes[0]: {}, names_genes[1]: {}, names_genes[2]: {}}

for gene in dic_genes.keys():
    term = 'Human[Orgn] AND ' + gene + '[Gene]'
    handle = E.esearch(db = 'nucleotide', term = term) #search in NCBI gene ids
    list_id_gene =  E.read(handle)['IdList']
    if len(list_id_gene) != 0:
        dic_genes[gene]['IdList_nt'] = list_id_gene; 
        dic_genes[gene]['Files'] = {} 
        for ID_gene in list_id_gene: #for each id, retrieve GeneBank file and save it
            list_id_pubmed = []
            file_name = 'GB_' + gene + '_' + ID_gene + '.txt'
            file = open(file_name, 'w')
            record = E.efetch(db = 'nucleotide', id = ID_gene, rettype = 'gb', retmode = 'text')
            record_string = record.read()
            file.write(record_string)
            file.close()
            gene_data = SeqIO.read(E.efetch(db = 'nucleotide', id = ID_gene, rettype = 'gb', retmode = 'text'), 'genbank')
            ptr_id = 0; transl_seq = 0; seq = 0
            try:
                for data in gene_data.features:
                    if data.type == 'CDS':
                        ptr_id = data.qualifiers['protein_id']
                        transl_seq = data.qualifiers['translation']
                for i in range(len(gene_data.features)):
                    if gene_data.features[i].type == 'CDS':
                        seq = gene_data.features[i].extract(gene_data.seq)
            except:
                ptr_id = 0; transl_seq = 0; seq = 0
            list_id_pubmed = []
            pos_pub = [(m.start(0), m.end(0)) for m in re.finditer('PUBMED', record_string)]
            for pos in pos_pub: #retrieve pubmed ids in GeneBank file
                list_id_pubmed += ((record_string[pos[1]+1:pos[1]+11]).split())
            dic_genes[gene]['Files'][file_name] = {'ids_pubmed': list_id_pubmed, 'Protein_id': ptr_id, 'Protein_seq': transl_seq, 'Seq': seq}
            if len(list_id_pubmed) != 0:
                file_name = 'Articles_' + gene + '_' + ID_gene + '.txt'
                file = open(file_name, 'w')
                handle = E.efetch(db = 'pubmed', id = list_id_pubmed, rettype = 'Medline', retmode = 'text')
                records = list(Medline.parse(handle))
                for record in records: #for each pubmed id, store following parameters
                    id_pubmed = str(record.get('PMID'))
                    title = str(record.get('TI'))
                    authors = str(record.get('AU'))
                    abstract = str(record.get('AB'))
                    source = str(record.get('SO'))
                    file.writelines('IDPubmed: ' + id_pubmed + '\nTitle: ' + title + '\nAuthors: ' + authors + '\nAbstract: ' + abstract + '\nSource: ' + source + '\n\n')
                file.close()
    else:
        dic_genes[gene]['IdList_nt'] = 'No gene IDs available'
        
## IDProtein SwissProt
for gene in dic_genes.keys():
    term = 'Human[Orgn] AND ' + gene + '[Gene]'
    handle = E.esearch(db = 'gene', term = term) #search in NCBI gene ids
    id_gene =  E.read(handle)['IdList']
    dic_genes[gene]['IdList_gene'] = [id_gene[0]] 
    handle = E.efetch(db = 'gene', id = id_gene[0], retmode = 'xml')
    result =  E.read(handle)
    if gene == 'SCARF2':
        dic_genes[gene]['IdList_gene'].append(result[0]['Entrezgene_comments'][4]['Gene-commentary_comment'][0]['Gene-commentary_products'][0]['Gene-commentary_products'][0]['Gene-commentary_comment'][2]['Gene-commentary_comment'][0]['Gene-commentary_source'][1]['Other-source_src']['Dbtag']['Dbtag_tag']['Object-id']['Object-id_str']) 
    elif gene == 'ZNF74':  
        dic_genes[gene]['IdList_gene'].append(result[0]['Entrezgene_comments'][4]['Gene-commentary_products'][19]['Gene-commentary_products'][0]['Gene-commentary_source'][1]['Other-source_anchor'][-6:])
    else:
         dic_genes[gene]['IdList_gene'].append(result[0]['Entrezgene_comments'][5]['Gene-commentary_products'][10]['Gene-commentary_products'][0]['Gene-commentary_source'][1]['Other-source_anchor'][-6:])
       
# ____________________________________ Blast ____________________________________
        
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
    
## SCARF 2

gene = 'SCARF2'
ids_prot = ['1677556846', '1677530789']

for id_prot in ids_prot:
    name_file_gb = 'GB_' + gene + '_' + id_prot + '.txt'
    file_name = gene + '_' + id_prot + '.fasta'
    file = open(file_name, 'w+')
    file.writelines('> ' + 'id ' + id_prot + '\n' + str(dic_genes[gene]['Files'][name_file_gb]['Protein_seq'])[2:-2])
    file.close()

for id_prot in ids_prot:
    file_fasta = gene + '_' + id_prot + '.fasta'
    name_file = 'Blast_' + gene + '_' + id_prot + '.txt'
    record = SeqIO.read(open(file_fasta), 'fasta')
    result_handle = NCBIWWW.qblast( 'blastp', 'nr', record.format('fasta'))
    blast_records = NCBIXML.parse(result_handle)
    file = open(name_file, 'w+')
    e_value_thresh = 0.05
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:
                    file.writelines('****Alignment****\n')
                    file.writelines('\nsequence:' + str(alignment.title))
                    file.writelines('\nlength:' + str(alignment.length))
                    file.writelines('\ne value:' + str(hsp.expect))
                    file.writelines( hsp.query [0:75] + '\n...\n')
                    file.writelines( hsp.match [0:75] + '...\n')
                    file.writelines( hsp.sbjct [0:75] + '...\n\n')
    file.close()

## ZNF74
gene = 'ZNF74'
id_prot = '1036030677'

file_gb = 'GB_' + gene + '_' + id_prot + '.txt'
file = open(file_gb)
text = file.read()
file.close()
re.findall('translation=', text)
pos_int_protein = [(m.start(0), m.end(0)) for m in re.finditer('translation="', text)]
pos_end_protein = [(m.start(0), m.end(0)) for m in re.finditer('"\n', text)]
for i in range(4):
    pos_init = pos_int_protein[i][1]
    for j in range(len(pos_end_protein)):   
        pos_end = pos_end_protein[j][0]
        if pos_end > pos_init:
            break    
    file_fasta = gene + '_' + id_prot + '_' + str(i+1) + '.fasta'
    file = open(file_fasta, 'w+')
    file.write('> id ' + id_prot + 'isoforma ' + str(i+1) + '\n' + ((text[pos_init:pos_end]).replace('\n', '')).replace(' ', ''))
    file.close()

for i in range(4):
    file_fasta = gene + '_' + id_prot + '_' + str(i+1) + '.fasta'
    name_file = 'Blast_' + gene + '_' + id_prot + '_' + str(i+1) + '.txt'
    record = SeqIO.read(open(file_fasta), 'fasta')
    result_handle = NCBIWWW.qblast( 'blastp', 'nr', record.format('fasta'))
    blast_records = NCBIXML.parse(result_handle)
    file = open(name_file, 'w+')
    e_value_thresh = 0.05
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:
                    file.writelines('****Alignment****\n')
                    file.writelines('sequence:' + str(alignment.title))
                    file.writelines('\nlength:' + str(alignment.length))
                    file.writelines('\ne value:' + str(hsp.expect))
                    file.writelines( '\n' + hsp.query [0:75] + '...\n')
                    file.writelines( hsp.match [0:75] + '...\n')
                    file.writelines( hsp.sbjct [0:75] + '...\n')
    file.close()


## KLHL22
    
gene = 'KLHL22'
id_prot = '1519245277'

# blastp
name_file_gb = 'GB_' + gene + '_' + id_prot + '.txt'
file_name = gene + '_' + id_prot + '.fasta'
file = open(file_name, 'w+')
file.writelines('> ' + 'id ' + id_prot + '\n' + str(dic_genes[gene]['Files'][name_file_gb]['Protein_seq'])[2:-2])
file.close()

file_fasta = gene + '_' + id_prot + '.fasta'
name_file = 'Blast_' + gene + '_' + id_prot + '.txt'
record = SeqIO.read(open(file_fasta), 'fasta')
result_handle = NCBIWWW.qblast( 'blastp', 'nr', record.format('fasta'))
blast_records = NCBIXML.parse(result_handle)
file = open(name_file, 'w+')
e_value_thresh = 0.05
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_thresh:
                file.writelines('****Alignment****\n')
                file.writelines('\nsequence:' + str(alignment.title))
                file.writelines('\nlength:' + str(alignment.length))
                file.writelines('\ne value:' + str(hsp.expect))
                file.writelines( hsp.query [0:75] + '\n...\n')
                file.writelines( hsp.match [0:75] + '...\n')
                file.writelines( hsp.sbjct [0:75] + '...\n\n')
file.close()


# ___________________ Multiple aligment and phylogenetic tree ___________________

## SCARF2
#manual analysis - choose ids 
id_alinm_1 = ['XP_016795221.2', 'XP_030771343.1', 'XP_024095947.1', 'XP_023046686.1']
id_alinm_2 = ['XP_016795222.2', 'XP_030861352.1', 'XP_024095948.1', 'XP_010367300.2' ]
id_alinm= ['NP_699165.3', 'NP_878315.2'] + id_alinm_1 + id_alinm_2

file = open('SCARF2_multiple.fasta', 'w+')
for id_al in id_alinm:
    record = E.efetch(db = 'protein', id = id_al, rettype = 'fasta', retmode = 'text')
    file.writelines(record.read())
file.close()

## KLHL22
#manual analysis - choose ids 
id_alinm= ['NP_116164.2', 'XP_024095566.1', 'XP_002830921.1', 'XP_003807728.1', 'NP_001248517.1', 'XP_023046688.1', 'XP_031512345.1']

file = open('KLHL22_multiple.fasta', 'w+')
for id_al in id_alinm:
    record = E.efetch(db = 'protein', id = id_al, rettype = 'fasta', retmode = 'text')
    file.writelines(record.read())
file.close()

## ZNF74
#manual analysis - choose ids 
id_alinm_1 = ['XP_008958064.1', 'XP_009436111.2', 'XP_010367301.2']
id_alinm_2 = ['XP_009232403.1', 'XP_011738108.1', 'XP_011738107.1']
id_alinm_3 = ['XP_008958065.1', 'XP_023046697.1', 'XP_023046696.1']
id_alinm_4 = ['XP_009436114.2', 'XP_030771342.1', 'XP_005567989.1']
id_alinm= ['AAF21777.1', 'AAF21778.1', 'AAF21779.1', 'AAF21780.1'] + id_alinm_1 + id_alinm_2 + id_alinm_3 + id_alinm_4

file = open('ZNF74_multiple.fasta', 'w+')
for id_al in id_alinm:
    record = E.efetch(db = 'protein', id = id_al, rettype = 'fasta', retmode = 'text')
    file.writelines(record.read())
file.close()

## ficheiros obtidos a partir da ferramenta docker - imagem ClustalW

## Trees
from Bio import Phylo

scarf_tree = Phylo.read("tree_scarf", "newick")
klhl22_tree = Phylo.read("KLHL22_multiple.dnd", "newick")
znf74_tree = Phylo.read("ZNF74_multiple.dnd", "newick")
Phylo.draw_ascii(znf74_tree)
