# This file plots activator regulons in a heatmap for Figure 4D and Supplemental Figure 4F
# Requires: scplusmdata.h5mu (from SCENIC+ pipeline via Snakemake)

import mudata
scplus_mdata = mudata.read("outs/scplusmdata.h5mu")

# Extended eregulon heatmap (Figure 4D and Supplemental Figure 4F)

meta_key = "extended_e_regulon_metadata"
er_meta = scplus_mdata.uns[meta_key].copy()

# keep only activators (+/+)
names = er_meta["eRegulon_name"].astype(str)
activators = names[names.str.contains('+/+', regex=False, na=False)]

from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "extended_gene_based_AUC",
    size_modality = "extended_region_based_AUC",
    group_variable = "scRNA_counts:ident",
    eRegulon_metadata_key = "extended_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "extended_gene_based_AUC",
    orientation = "vertical",
    subset_feature_names = activators,
    save = 'eregulon_extended_heatmap.pdf',
    figsize = (20, 20)
)


# Direct eregulon 
meta_key = "direct_e_regulon_metadata"
er_meta = scplus_mdata.uns[meta_key].copy()

# keep only activators (+/+)
names = er_meta["eRegulon_name"].astype(str)
activatorsDirect = names[names.str.contains('+/+', regex=False, na=False)]

heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:ident",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "vertical",
    save = "eregulon_direct_heatmap.pdf",
    subset_feature_names = activatorsDirect,
    figsize = (20, 20)
)