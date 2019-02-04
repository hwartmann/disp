import pandas as pd
import os
import re
import csv
import subprocess


def filter_variant_summary(variant_summary=None,
                           assembly='GRCh37',
                           clin_sig_simple=1,
                           not_phenotypelist=('not provided', 'not specified', 'See cases'),
                           var_type=('single nucleotide variant', 'indel', 'deletion', 'insertion'),
                           not_chromosome=('na', 'MT'),
                           max_length=50,
                           path_to_interim='../../data/interim/'):

    """return a filtered dataframe"""

    if not os.path.isfile(path_to_interim+'filtered_variant_summary.csv'):
        total_variants = variant_summary.shape[0]

        variant_summary = variant_summary[
            (variant_summary['Assembly'] == assembly) &
            (variant_summary['ClinSigSimple'] == clin_sig_simple) &
            (~variant_summary['PhenotypeList'].isin(not_phenotypelist)) &
            (variant_summary['Type'].isin(var_type)) &
            (~variant_summary['Chromosome'].isin(not_chromosome)) &
            ((variant_summary['Stop'] - variant_summary['Start'] + 1) <= max_length)
            ]

        selected_variants = variant_summary.shape[0]

        print("{end} of {beg} variant_summary entries have been selected".format(beg=total_variants,
                                                                                 end=selected_variants))

        variant_summary.to_csv(path_to_interim+'filtered_variant_summary.csv')

        return variant_summary
    else:
        return pd.read_csv(path_to_interim+'filtered_variant_summary.csv')


def add_coid(variant_summary=None,
             path_to_interim='../../data/interim/'):

    """A unique identifier for each variant is necessary, clinvar does not provide them so we create one from the
     given coordinates and add it as a column"""

    if not os.path.isfile(path_to_interim + 'COID_filtered_variant_summary.csv'):
        variant_summary['COID'] = variant_summary['Chromosome'].astype(str) + "-" +\
                                  variant_summary['Start'].astype(str).str[-7:] + "-" +\
                                  variant_summary['Stop'].astype(str).str[-7:]

        variant_summary.to_csv(path_to_interim+'COID_filtered_variant_summary.csv')

        return variant_summary
    else:
        return pd.read_csv(path_to_interim + 'COID_filtered_variant_summary.csv')


def get_phenotype_list(semicolon_positions=None, phenotype_string=None):
    """
    takes a list of phenotypes as ; seperated string and returns a python list
    """
    phenotype_list = []
    pt_start = 0

    for pt_end in semicolon_positions:
        phenotype_list.append(phenotype_string[pt_start:pt_end])
        pt_start = pt_end + 1

    phenotype_list.append(phenotype_string[pt_start:])

    return phenotype_list


def expand_variant_summary(variant_summary=None,
                           path_to_interim='../../data/interim/'):
    """
    takes in a data frame with the variant summary
    loops through each entry and looks for the semicolons that are used for phenotype and concept id list entries
    when a semicolon found the df id of the entry is stored in a list to pop later
    all semicolon positions are passed to get_phenotype_list which returns a list with all phenotypes for that variant
    a new entry is written into the df
    :param variant_summary:
    :return:
    """
    if not os.path.isfile(path_to_interim + 'expanded_variant_summary.csv'):

        idx_pop = []
        df_tmp = pd.DataFrame()
        total = variant_summary.shape[0]
        progress = 0

        for idx, row in variant_summary.iterrows():
            if progress % 10000 == 0:
                print("clinvar expand {}%".format(round(progress / total, 2) * 100))

            semi_type = [m.start() for m in re.finditer(';', row['PhenotypeList'])]
            semi_id = [m.start() for m in re.finditer(';', row['PhenotypeIDS'])]

            if len(semi_type) > 0:
                idx_pop.append(idx)
                samp = get_phenotype_list(semi_type, row[13])
                ids = get_phenotype_list(semi_id, row[12])
                for i in range(len(semi_type)):

                    df_tmp = df_tmp.append({'#AlleleID': row[0],
                                            'Type': row[1],
                                            'Name': row[2],
                                            'GeneID': row[3],
                                            'GeneSymbol': row[4],
                                            'HGNC_ID': row[5],
                                            'ClinicalSignificance': row[6],
                                            'ClinSigSimple': row[7],
                                            'LastEvaluated': row[8],
                                            'RS# (dbSNP)': row[9],
                                            'nsv/esv (dbVar)': row[10],
                                            'RCVaccession': row[11],
                                            'PhenotypeIDS': ids[i],
                                            'PhenotypeList': samp[i],
                                            'Origin': row[14], 'OriginSimple': row[15], 'Assembly': row[16],
                                            'ChromosomeAccession': row[17],
                                            'Chromosome': str(row[18]),
                                            'Start': str(row[19]),
                                            'Stop': str(row[20]),
                                            'ReferenceAllele': row[21],
                                            'AlternateAllele': row[22],
                                            'Cytogenetic': row[23],
                                            'ReviewStatus': row[24],
                                            'NumberSubmitters': row[25],
                                            'Guidelines': row[26],
                                            'TestedInGTR': row[27],
                                            'OtherIDs': row[28],
                                            'SubmitterCategories': row[29],
                                            'VariationID': row[30],
                                            'COID': row[31]
                                            }, ignore_index=True)

            progress += 1

        variant_summary.drop(idx_pop, inplace=True)
        variant_summary = variant_summary.append(df_tmp)
        variant_summary.to_csv(path_to_interim+'expanded_variant_summary.csv')

        return variant_summary

    else:
        return pd.read_csv(path_to_interim + 'expanded_variant_summary.csv')


def make_variant_bed(variant_summary=None,
                    path_to_interim='../../data/interim/'):
    """
    open new file, iterate variant_summary and write new bed style row for each entry in variant summary

    :param variant_summary:
    :param path_to_interim:
    :return:
    """
    if not os.path.isfile(path_to_interim + 'bed_variant_summary.bed'):

        with open(path_to_interim+'bed_variant_summary.bed', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')

            for idx, row in variant_summary.iterrows():
                writer.writerow(["chr"+str(row['Chromosome']),
                                 str(row['Start']),
                                 str(row['Stop']),
                                 str(row['PhenotypeList']),
                                 str(row['PhenotypeIDS']),
                                 str(row['GeneID']),
                                 str(row['GeneSymbol']),
                                 str(row['COID']),
                                 "rs"+str(row['RS# (dbSNP)'])
                                 ]
                                )


def make_gtex_bed(path_gtex_exon_annotation=None,
                  path_to_interim='../../data/interim/'):

    if not os.path.isfile(path_to_interim + 'bed_gtex_exons.bed'):
        gtex = pd.read_csv(path_gtex_exon_annotation,
                           usecols=['chr', 'start_pos', 'end_pos', 'exon_id'],
                           delimiter='\t',
                           low_memory=False)

        gtex['chr'] = "chr" + gtex['chr'].astype(str)
        gtex.to_csv(path_to_interim+'bed_gtex_exons.bed', index=False, header=False, sep='\t')


def intersect(file_a=None,
              file_b=None,
              path_to_interim='../../data/interim/'):
    """
    run bedtools to intersect variants with exons. Needs absolute paths!

    :param file_a: bed_gtex_exons.bed
    :param file_b: bed_variant_summary.bed
    :param path_to_interim:
    :return:
    """
    if not os.path.isfile(path_to_interim+'bed_gtex_clinvar_intersect.bed'):
        subprocess.run('bedtools intersect -a {file_a} \
             -b {file_b} -wa -wb \
                  > {target}'.format(file_a=file_a,
                                     file_b=file_b,
                                     target=path_to_interim+'bed_gtex_clinvar_intersect.bed'),
                       shell=True)


if __name__ == '__main__':
    path_root_interim = '/mnt/wartmann_pe940/git/disp/data/interim/'
    path_lvl = '../../'
    path_interim = path_lvl+'data/interim/from_intersect_clinvar_gtex.py/'
    path_variant_summary = path_lvl+'data/raw/variant_summary.txt'
    path_gtex_exon_annotation = path_lvl+'data/raw/gencode.v19.genes.v7.patched_contigs.exons.txt'

    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
        1. clinvar summary file is loaded and filtered
        2. COID (unique coordinate/identifier id for each variant) is added
        3. clinvar summary file is expanded, some variants contain lists of phenotype, we make each list item a single
            entry such that we can map it to the ontology
        4. make a .bed file out of the expanded variant summary for later intersection with GTEx exon coordinates
        5. make a .bed file out of gencode.v19.genes.v7.patched_contigs.exons.txt
        6. make an intersection file using bedtools, the variant summary bed and the gtex exon bed
    """
    # TODO: add tmp_file_name as function arg to all functions that save a tmp file to interim

    variant_summary = pd.read_csv(path_variant_summary, delimiter='\t')
    variant_summary = filter_variant_summary(variant_summary, path_to_interim=path_interim)
    variant_summary = add_coid(variant_summary, path_to_interim=path_interim)
    variant_summary = expand_variant_summary(variant_summary, path_to_interim=path_interim)

    make_variant_bed(variant_summary, path_to_interim=path_interim)
    make_gtex_bed(path_gtex_exon_annotation, path_to_interim=path_interim)

    intersect(path_root_interim+'bed_gtex_exons.bed',
              path_root_interim+'bed_variant_summary.bed',
              path_to_interim=path_interim)