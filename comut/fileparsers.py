def parse_maf(maf_df, variants=None, rename=None):
    '''Subsets a MAF to nonsilent mutations, then subsets to
    required columns and renames. Returns a parsed datafame.

    Params:
    -------
    maf_df: pandas dataframe
        The dataframe of a MAF after pd.read_csv()

    variants: list-like
        A list of variants to subset the MAF to. Defaults to
        nonsilent variants.

    rename: dict
        A dictionary to rename values in the MAF. Defaults to
        collapsing insertions and deletions and renaming
        nonsilent mutations.

    Returns:
    --------
    parsed_maf: pandas_dataframe
        A parsed maf for CoMut. It subsets to sample name, gene
        name, and variant classification and renames the columns.
    '''

    # default to nonsilent variants
    if variants is None:
        variants = ['Nonsense_Mutation', 'In_Frame_Del', 'Frame_Shift_Ins', 'Splice_Site', 'In_Frame_Ins', 'Frame_Shift_Del', 'Missense_Mutation']

    # default rename is to collapse ins and del to indel
    if rename is None:
        rename = {'Nonsense_Mutation': 'Nonsense', 'In_Frame_Del': 'In frame indel', 'In_Frame_Ins': 'In frame indel',
                  'Frame_Shift_Del': 'Frameshift indel', 'Missense_Mutation': 'Missense', 'Splice_Site': 'Splice site', 'Frame_Shift_Ins': 'Frameshift indel'}

    # subset to required columns
    subset_maf_df = maf_df[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification']]

    # rename columns
    subset_maf_df.columns = ['sample', 'category', 'value']

    # subset maf to relevant mutation types
    parsed_maf = subset_maf_df[subset_maf_df['value'].isin(variants)].copy()

    # rename variants
    parsed_maf.loc[:, 'value'] = parsed_maf.loc[:, 'value'].replace(rename).copy()

    return parsed_maf
