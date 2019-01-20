#!usr/bin/env python

"""Annotation.py: create a variant analysis tool. """
# When provided with a VCF file, write a small program to output a table annotating each variant in the file. The annotation must include: the type of variation, depth of sequencing coverage, number of reads supporting the variant, percentage of reads supporting the variant, allele frequency from the ExAC project API, and addtional optional information from ExAC."""

#__author__ = "Terri Driessen"
#__email__ = terri.driessen@gmail.com
#__date__ = 1.20.2019
#__python__ = python 3.7


# import pandas, io, re, codecs, json, requests, and numpy
import pandas as pd
import io
import re
import codecs, json
import requests
import numpy as np


# write function to read the supplied vcf file and edit the columns of the vcf file.
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_table(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})

# read vcf file as a pd dataframe
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
output_final=output.merge(variant_allele_concat, on=['id','chrom_position'],how='outer').drop('id', axis=1)
pd.DataFrame.to_csv(output_final, path_or_buf="Annotation_output.csv", sep=',')
