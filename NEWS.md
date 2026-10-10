# ggsc 1.11.3

+ `sc_dim_geom_ellipse()` and `sc_dim_geom_label()` now report that the group column could not be determined when the plot has no non-numeric column to group by (e.g. `sc_dim()` without `colLabels()`); the check indexed the plot data with the missing column name and failed with `attempt to select less than one element in get1index` before the message could be displayed (2026-10-10, Sat)
+ `sc_dim_geom_ellipse()` only passes `level` to layer functions that accept it, so the documented `geom = ggforce::geom_mark_hull` no longer warns about an unknown parameter (2026-10-09, Fri)
+ `sc_dim_geom_label()` places labels at the true group centre; the reference point was averaged from a sampled ellipse outline, which was offset from the centre and drifted with the ellipse `level` (2026-10-09, Fri)
+ `sc_dim_geom_label(data = )` now reports that `data` must not be supplied, instead of an unrelated message about the label mapping (2026-10-09, Fri)

# ggsc 1.11.2

+ compatible with new S4Vectors (2026-09-10, Thu)

# ggsc 1.11.1

+ compatible with Seurat V5 (2026-05-17, Sun)

# ggsc 1.10.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# ggsc 1.8.0

+ Bioconductor RELEASE_3_22 (2025-11-01, Sat)

# ggsc 1.7.1

+ compatible with new ggplot2 > 3.5.2 (2025-07-02, Wed)

# ggsc 1.6.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# ggsc 1.4.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# ggsc 1.3.1

+ add `plot_lisa_feature()` (2024-09-06, #34, #39)
+ add `geom_bgpoint()` layer (2024-06-18, Tue, #24)

# ggsc 1.2.0

+ Bioconductor RELEASE_3_19 (2024-05-01, Wed)

# ggsc 1.1.4

+ add background outline of umap plot (2024-04-12, Fri, #22)

# ggsc 1.1.3

+ support plotting pie for spatial data (2023-12-13, Wed, #18)
+ extract `meta.data` of Seurat object (2023-12-12, Tue, #17)

# ggsc 1.1.2

+ add `sc_dot()` methods (2023-11-29, Wed, #15)
+ update vignette to link to the online docs (<https://yulab-smu.top/ggsc>)
+ add `README.Rmd` and `README.md` in github repo
+ optimization: retrieve embedding without `FetchData` (2023-11-27, Mon, #14)
+ bug fixed for Seurat object (2023-10-31, Tue, #12, #13)

# ggsc 1.1.1

+ ignore the tissue section when image is not exist (2023-10-31, Tue, #12)
+ introduce `joint` to combine all features with `joint.fun` and speed up calculation of kde using RcppParallel (2023-10-25, Wed, #11)

# ggsc 1.0.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# ggsc 0.99.11

+ support density visualization for single and spatial transcriptomic data (2023-10-18, Wed)

# ggsc 0.99.10

+ on Bioconductor (2023-10-16, Mon)
+ add `\value` session in the `reexports.Rd` (2023-10-15, Sun)
+ add package level man page and update vignette (2023-10-14, Sat)
+ add examples in Rd to satisfy BiocCheck (2023-09-18, Mon, #7)
+ `sc_dim_count()` function to generate a barplot from a dimension reduction plot (`sc_dim()` plot) to 
    visualize the number of cells for each clusters (2023-09-13, Wed)
+ add 'biocViews' in DESCRIPTION required by Bioconductor

# ggsc 0.99.0

+ compatible with 'SingleCellExperiment' (2023-09-05, Tue, #5)
+ using S4 OOP to reorganize the functions (2023-09-05, Tue, #4)
+ rename the package to 'ggsc' as there is a package called 'scplot' in CRAN
+ add H&E image to `sc_spatial()` (#3)

# scplot 0.0.3

+ `sc_spatial` to visualize spatial features (2022-12-07, Wed)

# scplot 0.0.2

+ `sc_dim_geom_sub` and `sc_dim_sub` (2022-12-03, Sat)
+ `sc_dim_geom_ellipse` to draw ellipse on `sc_dim()` (2022-12-02, Fri)

# scplot 0.0.1

+ several functions implemented (2022-11-09, Wed)
    - `sc_dim`
    - `sc_dim_geom_feature`
    - `sc_dim_geom_label`
    - `sc_feature`
    - `sc_geom_point`
    - `sc_violin`
