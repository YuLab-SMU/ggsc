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

+ `sc_dim`, `sc_dim_geom_feature`, `sc_dim_geom_label`, `sc_feature`, `sc_geom_point` and `sc_violin` functions (2022-11-09, Wed)
