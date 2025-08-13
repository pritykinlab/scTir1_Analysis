import pandas as pd
import numpy as np


def compute_deg_annotations_signed(
    deseq_l2fc_df: pd.DataFrame,
    deseq_padj_df: pd.DataFrame,
    tissue1,
    tissue2='SLO',
    strict_padj_cutoff: float = 0.05,
    loose_padj_cutoff: float = 0.2,
    lfc_agreement_cutoff = 0.2,
) -> pd.DataFrame:
    genes = deseq_l2fc_df.index

    all_annotations = pd.DataFrame(index=genes)

    # masks
    strict_tissue = (deseq_padj_df[tissue1] < strict_padj_cutoff)
    # print(strict_tissue.sum())
    loose_tissue = (deseq_padj_df[tissue1] < loose_padj_cutoff)
    strict_slo = (deseq_padj_df[tissue2] < strict_padj_cutoff)
    loose_slo = (deseq_padj_df[tissue2] < loose_padj_cutoff)

    lfc_tissue = deseq_l2fc_df[tissue1].fillna(0)
    lfc_slo = deseq_l2fc_df[tissue2].fillna(0)

    # shared direction

    tissue1_DEG_shared_with_SLO = (  (strict_tissue & loose_slo & (np.sign(lfc_tissue) == np.sign(lfc_slo))) 
                             | (strict_tissue & (lfc_tissue > 0) & (lfc_slo > lfc_agreement_cutoff))
                             | (strict_tissue & (lfc_tissue < 0) & (lfc_slo < -lfc_agreement_cutoff))
    )
    SLO_DEG_shared_with_tissue2 = (  
                               (strict_slo & loose_tissue & (np.sign(lfc_tissue) == np.sign(lfc_slo))) 
                             | (strict_slo & (lfc_slo > 0) & (lfc_tissue > lfc_agreement_cutoff))
                             | (strict_slo & (lfc_slo < 0) & (lfc_tissue < -lfc_agreement_cutoff))
    )
    agreement = (tissue1_DEG_shared_with_SLO | SLO_DEG_shared_with_tissue2)
    SLO_deg_alone = strict_slo & (~agreement)
    tissue1_deg_alone = strict_tissue & (~agreement)

    all_annotations.loc[agreement & (lfc_slo > 0), 'Label'] = 'Shared up'
    all_annotations.loc[agreement & (lfc_slo < 0), 'Label'] = 'Shared down'
    
    all_annotations.loc[tissue1_deg_alone & (lfc_tissue > 0), 'Label'] = f'{tissue1} up alone'
    all_annotations.loc[tissue1_deg_alone & (lfc_tissue < 0), 'Label'] = f'{tissue1} down alone'    

    all_annotations.loc[SLO_deg_alone & (lfc_slo > 0), 'Label'] = f'{tissue2} up alone'
    all_annotations.loc[SLO_deg_alone & (lfc_slo < 0), 'Label'] = f'{tissue2} down alone'    
    
    opposite_direction = ((strict_tissue & loose_slo & (np.sign(lfc_tissue) != np.sign(lfc_slo)))
                         |(loose_tissue & strict_slo & (np.sign(lfc_tissue) != np.sign(lfc_slo)))
    )

    all_annotations.loc[opposite_direction, 'Label'] = f'Opposite'
    return all_annotations


import pandas as pd
import numpy as np

def compute_deg_annotations_shared_one_pvalue(
    deseq_l2fc_df: pd.DataFrame,
    deseq_padj_df: pd.DataFrame,
    tissue1,
    tissue2='SLO',
    padj_cutoff: float = 0.05,
) -> pd.DataFrame:
    genes = deseq_l2fc_df.index
    all_annotations = pd.DataFrame(index=genes)

    # significance masks
    sig_tissue1 = deseq_padj_df[tissue1] < padj_cutoff
    sig_tissue2 = deseq_padj_df[tissue2] < padj_cutoff

    # fold change
    lfc_tissue1 = deseq_l2fc_df[tissue1].fillna(0)
    lfc_tissue2 = deseq_l2fc_df[tissue2].fillna(0)

    # shared: significant in both + same direction
    shared = sig_tissue1 & sig_tissue2 & (np.sign(lfc_tissue1) == np.sign(lfc_tissue2))

    # opposite: significant in both + opposite direction
    opposite = sig_tissue1 & sig_tissue2 & (np.sign(lfc_tissue1) != np.sign(lfc_tissue2))

    # tissue-specific DEGs
    tissue1_only = sig_tissue1 & ~sig_tissue2
    tissue2_only = sig_tissue2 & ~sig_tissue1

    # annotate
    all_annotations.loc[shared & (lfc_tissue1 > 0), 'Label'] = 'Shared up'
    all_annotations.loc[shared & (lfc_tissue1 < 0), 'Label'] = 'Shared down'

    all_annotations.loc[opposite, 'Label'] = 'Opposite'

    all_annotations.loc[tissue1_only & (lfc_tissue1 > 0), 'Label'] = f'{tissue1} up alone'
    all_annotations.loc[tissue1_only & (lfc_tissue1 < 0), 'Label'] = f'{tissue1} down alone'

    all_annotations.loc[tissue2_only & (lfc_tissue2 > 0), 'Label'] = f'{tissue2} up alone'
    all_annotations.loc[tissue2_only & (lfc_tissue2 < 0), 'Label'] = f'{tissue2} down alone'

    return all_annotations




def calc_figsize(n_rows, row_height=0.3, width=3):
    # row_height in inches per row
    height = max(1, row_height)  # ensure at least minimal height
    return (width, height*n_rows)


def make_gene_to_set(all_annotations):
    gene_to_set = {}
    for _, row in all_annotations[~(all_annotations == 'Shared up').where(~all_annotations.isna(), True).all(axis=1)].iterrows():
        conditions_with_deg = set()
        for col in row.index:
            cond1, cond2 = col.split('_')
            if isinstance(row[col], float):
                if np.isnan(row[col]):
                    continue
                
            if row[col] == 'Shared down':
                conditions_with_deg.add(cond1)
                conditions_with_deg.add(cond2)
    
            if row[col] == 'Shared up':
                conditions_with_deg.add(cond1)
                conditions_with_deg.add(cond2)
    
            if ' up alone' in row[col]:
                print(row[col])
                cond1 = row[col].split(' up alone')[0]
                conditions_with_deg.add(cond1)
    
            if ' down alone' in row[col]:
                print(row[col])
                cond1 = row[col].split(' down alone')[0]
                conditions_with_deg.add(cond1)
        
        gene_to_set[_] = conditions_with_deg
    return gene_to_set


def dict_to_annotation_df(gene_dict):
    """
    Convert dict {gene: set of annotations} to a DataFrame where each annotation is a column.
    """
    # Get all unique annotations across all genes
    all_annotations = sorted({ann for anns in gene_dict.values() for ann in anns})
    
    # Build the DataFrame
    data = []
    for gene, anns in gene_dict.items():
        row = {ann: (ann in anns) for ann in all_annotations}
        row['gene'] = gene
        data.append(row)
        
    df = pd.DataFrame(data)
    df = df.set_index('gene')
    return df
