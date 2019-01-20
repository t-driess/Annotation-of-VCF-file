# read in the VCF file and change the column names
# later on use the read_vcf function when calling the data frame

#!users/

# import pandas, io, re, codecs, json, requests, and numpy
import pandas as pd
import io
import re
import codecs, json
import requests
import numpy as np


# function to read the vcf file and edit columns
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_table(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})

# read vcf file as Pandas dataframe
vcf = pd.DataFrame(read_vcf('Challenge_data.vcf'))

# extract the location of variant and ref/alt allele, and rename dataframe as 'variant'
chrom_pos = pd.DataFrame(vcf.apply(lambda x: '%s-%s-%s-%s' %(x['CHROM'],x['POS'],x['REF'],x['ALT']), axis=1))
chrom_pos=chrom_pos.rename(columns={0:'chrom_position'})

# The type of variant is found in INFO TYPE. Extract the TYPE of variant and rename column as 'variant_type'.
variant_type = pd.DataFrame(vcf.INFO.apply(lambda x: re.findall("TYPE=(?=\w+)\w+",x)))
variant_type=variant_type.rename(columns={'INFO':'variant_type'})

# The depth of sequence coverage at the locus is found at INFO=<ID=DP. Make the sequence coverage it's own dataframe and rename it as "depth_coverage".
depth_coverage = pd.DataFrame(vcf.INFO.apply(lambda x: re.findall("DP=(?=\w+)\w+", x)))
depth_coverage=depth_coverage.rename(columns={'INFO':'depth_coverage'})

# The number of reads supporting the variant is INFO=<ID=AO. Make the number of reads supporting
# the variant it's own dataframe and rename it 'reads_for_variant'.
reads_for_variant = pd.DataFrame(vcf.INFO.apply(lambda x: re.findall(";AO=(?=\w+)\w+", x)))
reads_for_variant=reads_for_variant.rename(columns={'INFO':'reads_for_variant'})

# read from ExAC database. Convert variant chromosomal location and ref/alt alleles as a string, then call POST with a request body to retrieve a JSON array with variant strings from ExAC.
chrom = vcf.apply(lambda x: '%s-%s-%s-%s' %(x['CHROM'],x['POS'],x['REF'],x['ALT']), axis=1)
chrom_for_payload=chrom.values.tolist()
payload = json.dumps(chrom_for_payload)
r=requests.post("http://exac.hms.harvard.edu/rest/bulk/variant/variant", data=payload)
x_json=r.json()
y_json=str(x_json)

# make a dataframe for the allele frequencies. Split the dataframe into two columns, and keep the allele frequency values.
allele_freq = pd.DataFrame(re.findall('allele_freq\D\D\D(?=\w+)\w+\D(?=\w+)\w+',y_json))
allele_freq2 = allele_freq.rename(columns={0:'allele'})
new2 = allele_freq2['allele'].str.split(": ",n=1, expand=True)
allele_freq2['allele_f']=new2[0]
allele_freq2['allele_freq']=new2[1]
allele_freq2.drop(columns=['allele','allele_f'], inplace=True)

# make a dataframe for the variant
variant_id=pd.DataFrame(re.findall('.variant_id\D\D\D\D....(?=\w+)\w+....',y_json))
variant_id2 = variant_id.rename(columns={0:'variant_id'})
new=variant_id2['variant_id'].str.split(": '", n=1, expand=True)
variant_id2['variant']=new[0]
variant_id2['chrom_position']=new[1]
variant_id2.drop(columns=['variant_id','variant'], inplace=True)

variant_allele_concat=pd.concat([variant_id2, allele_freq2], axis=1)
variant_allele_concat['id']=variant_allele_concat.groupby('chrom_position').cumcount()

# assemble the dataframe and export as CSV
output = pd.concat([chrom_pos, variant_type, depth_coverage, reads_for_variant], axis=1)
output['id']=output.groupby('chrom_position').cumcount()
#output_final=output.merge(x, on=['actual_variant'],how='outer').drop('id',axis=1)

output_final=output.merge(variant_allele_concat, on=['id','chrom_position'],how='outer').drop('id', axis=1)
pd.DataFrame.to_csv(output_final, path_or_buf="Annotation_output.csv", sep=',')


#    if variant_id=y['variant_id']:
#        return 'allele_freq\D\D\D(?=\w+)\w+\D(?=\w+)\w+'
#    else:
#        return 'Nan'

#re.findall('allele_f(?=\w+)\w+\D\D\D\D(?=\w+)\w+',y) #gets allele_freq' : 0
#z=re.findall("'allele_f(?=w+)\w+",y) #returns allele_freq for all "DP=(?=\w+)\w+"
#return(json.dumps(x['consequence']))

#allele_freq': 0.7955390334572491
#'variant_id': '4-126240510-T-C'



#import requests
#import json
#mystring=["14-21853913-T-C","22-46615746-A-G"]
#chrom_pos2 = chrom_pos.apply(lambda x: ', '.join(x.astype(str)) if x.dtype=='int64' else ', '.join("\'"+x.astype(str)+"\'"), reduce=False)

#payload = json.dumps(mystring)
#payload=chrom_pos.to_json(orient="records" or split or index or table or columns or values) doesn't work
# values and columns got the closest to what I want though...
#payload=chrom_pos.to_json(orient="split" or "table" with index=False, gets closer..

#chrom.values returns an array that is damn close.


#def function(x, variant_id):
#    return [obj in obj in x if obj['variant_id']==variant_id[0]['allele_freq']]

#r = requests.post("http://exac.hms.harvard.edu/rest/bulk/variant/variant", data=payload)




#def make_URL(path):
#    with open(chrom2, 'r') as f:
#        lines = [l for l in f if not l.startswith('0')]
#    for lines in f:
#        return('http://exac'+lines)
#this one works, and chrom should be a string.
#def URL(x):
#    for x in chrom:
#        return('http://exac.hms.harvard.edu/rest/variant/variant/'+str(x))

#def do_itt(x):
#    r=requests.get(url=URL(chrom))
#    x=r.json()
#    return(json.dumps(x['consequence']))

#r=requests.get(url=URL(chrom))

#z=list(chrom)
#def URL3(z):
#    for i in z:
#        URL2='http://exac.hms.harvard.edu/rest/variant/variant/'
#        URL3=URL2+i
#        r=requests.get(url=URL3)
#        x=r.json()
#        return(json.dumps(x['allele_freq']))

#def URL(x):
#    URL2='http://exac.hms.harvard.edu/rest/bulk/variant/variant/'
#    for x in chrom_pos['chrome']:
#        URL3=URL2+x
#        r=requests.get(url=URL3)
#        x=r.json()
#        return(json.dumps(x['allele_freq']))

#chrom=str("14-21853913-T-C","14-21854139-A-G")
#x2 = chrom.to_json(orient='split')
#def URL3(x):
#    for i in URL(chrom):
#        URL2='http://exac.hms.harvard.edu/rest/variant/variant/'
#        URL3=URL2+i
#        r=requests.get(url=URL3)
#        x=r.json()
#        return(json.dumps(x['allele_freq']))
#def URL3(x):
#    for x in chrom2:
#        URL2='http://exac.hms.harvard.edu/rest/variant/variant/'
#        return(URL2+x)


#def URL(x):
#    for x in chrom:
#        return('http://exac.hms.harvard.edu/rest/variant/variant'+str(x))

#r = requests.get(url=URL(chrom))
#x = r.json()
#json_str = json.dumps(x['allele_freq'])


#def CHROME_POS(i):
#    return i in chrom_pos

#newDF = CHROME_POS(chrom_pos, '1-931393')


#    chrome = pd.DataFrame(vcf.CHROM)
#    POS = pd.DataFrame(vcf.POS)
#    return str(chrome + "-" + POS)


#def ExAC(chrom_pos):
#    for i in chrom_pos:
#        return "http:///variant/variant/",chrom_pos


#        exac.hms.harvard.edu/rest
#x = pd.DataFrame(chrom_pos.apply(ExAC, axis=1))

#URL = "http://exac.hms.harvard.edu/rest/variant/variant/ExAC(variant_position)"
#r = requests.get(url=URL)
#x = r.json()
#json_str = json.dumps(x['allele_freq'])



# for FORMAT AO
#>>> a2['FORMAT'] = a2['FORMAT'].astype(str)
#>>> sep = pd.concat([a2, a2['FORMAT'].str.split(':', expand=True)], axis=1)

# The percent reads supporting the variant is AO/AO+RO *100



# To get to the ExAC database,
#import requests
#URL = "http://exac.hms.harvard.edu/rest/awesome?query=CHD8&service=variants_in_gene"
#R = requests.get(url=URL)
#R.text


#>>> import re
#>>> s = "allele_freq"
#>>> result = re.search(s, data.text)
#>>> result
#>>> print(result.group(0))

#Data = r.json()

#>>> import json
#>>> import os
#>>> os.chdir('/Users/terridriessen/Downloads')

#x = pd.read_json(path_or_buf=URL)

#Data will print out the information






# check and see what class your data frame is
#>>> print(type(x2.loc[x2.index[0], 'INFO']))




#a = y.filter(like='INFO', axis=1)

### practice

#a=pd.read_csv(‘Workbook2.csv’)
#list(a)
#sequencing depth at site
#b=a.INFO.apply(lambda x: re.findall(‘DP=(?=\w+)\w+’, x))

### EXTRA


#a=pd.read_csv(‘Workbook2.csv’)
#list(a)
#a[‘INFO2’]=a.INFO.apply(lambda x: re.findall(‘DP=(?=\w+)\w+’, x))
#s=a.apply(lambda x: pd.Series(a[‘INFO2’], axis=1).stack().reset_index(level=1, drop=False))

#b=a.filter(like=‘INFO’, axis=1)
#a[‘INFO2’]
