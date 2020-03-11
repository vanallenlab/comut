def parse_maf(maf_df):
    '''Subsets a MAF to nonsilent mutations, then subsets to
    required columns and renames. Returns a parsed datafame.

    Params:
    -------
    maf_df: pandas dataframe
        The dataframe of a MAF after pd.read_csv()

    Returns:
    --------
    nonsilent_maf: pandas_dataframe
        A parsed maf for CoMut. It subsets to sample name, gene
        name, and variant classification. It also subsets to nonsilent
        mutations, and renames these mutations for clarity, collapsing
        insertions and deletions into indels. It also renames columns
        to what CoMut expects (sample, category, value)
    '''

    # subset to required columns
    subset_maf_df = maf_df[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']]

    # rename columns
    subset_maf_df.columns = ['sample', 'category', 'value']

    # subset maf to relevant mutation types
    comut_variants = ['Nonsense_Mutation', 'In_Frame_Del', 'Frame_Shift_Ins', 'Splice_Site', 'In_Frame_Ins', 'Frame_Shift_Del', 'Missense_Mutation']
    nonsilent_maf = subset_maf_df[subset_maf_df['value'].isin(comut_variants)]

    # Rename variants, collapsing deletions and insertions into indels
    rename_dict = {'Nonsense_Mutation': 'Nonsense', 'In_Frame_Del': 'In frame indel', 'In_Frame_Ins': 'In frame indel',
                   'Frame_Shift_Del': 'Frameshift indel', 'Missense_Mutation': 'Missense', 'Splice_Site': 'Splice site', 'Frame_Shift_Ins': 'Frameshift indel'}
    nonsilent_maf = nonsilent_maf.replace(rename_dict)

    return nonsilent_maf
