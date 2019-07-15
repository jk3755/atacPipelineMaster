########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
# Rules for manually building/rebuilding the directory structure, for troubleshooting, etc.
#
rule buildir_lncap:
    input:
        "pros/lncap/ex01/cr01/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/cr02/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/cr04/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/cr05/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/cr07/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/cr08/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/wt01/operations/preprocessing/dirtree.built",
        "pros/lncap/ex01/wt02/operations/preprocessing/dirtree.built"