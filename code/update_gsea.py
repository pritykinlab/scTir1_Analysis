import gseapy as gp
def make_gsea_analysis(genes_down, baseline, skip_chr = True, gmt = None):
    def get_names(d):
        newnames = []
        for x in d['Term'].values:
            x = x.replace("GOBP_", "")
            x = x.replace("GOCC_", "")
            x = x.replace("REACTOME_", "")
            newnames.append(x)
        return newnames    


    from gseapy import Msigdb
    import gseapy as gp
    
    msig = Msigdb()
    if gmt is None:
        gmt = gp.read_gmt('/Genomics/argo/users/gdolsten/pritlab/genesets/msigdb.v2024.1.Mm.symbols.gmt')
    
    enr_down = gp.enrich(gene_list=list(genes_down),
                     gene_sets=gmt, 
                     background=list(baseline),
                     outdir=None,
                     verbose=True)

    result_down = enr_down.results

    if skip_chr:
        result_down = result_down[~(result_down['Term'].str.contains("chr") & (result_down['Term'].apply(len) < 8))]
    
    
    # fig, axs = init_subplots_exact(2, 1, fgsz=(1, 6), xspace=1.4, dpi = 150)
    # plt.sca(axs[0])
    # plt.scatter([0]*len(d), get_names(d), s = -np.log10(d['Adjusted P-value'])*10, c = d['Odds Ratio'], cmap='coolwarm', vmin=-50, vmax=50)
    # plt.xticks([])
    # plt.yticks(fontsize=12);
    # axs[0].set_title("Foxp3-up")
    
    # d = result_down.sort_values("P-value").iloc[:20].iloc[::-1]
    # plt.sca(axs[1])
    # plt.scatter([0]*len(d), get_names(d), s = -np.log10(d['Adjusted P-value'])*10, c = d['Odds Ratio'], cmap='coolwarm', vmin=-50, vmax=50)
    # axs[1].tick_params(left=False, labelleft=False, right=True, labelright=True)
    # plt.xticks([])
    # plt.yticks(fontsize=12);
    # axs[1].set_title("Foxp3-down")
    return result_down