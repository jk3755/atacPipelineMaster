########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
# rule preprocessing_test:
#     input:
#         "test/operations/preprocessing/test-REP1.preprocessing.complete",
#         "test/operations/preprocessing/test-REP2.preprocessing.complete",
#         "test/operations/preprocessing/test-REP3.preprocessing.complete"

########################################################################################################################################
#### H508 ##############################################################################################################################
########################################################################################################################################
rule full_analysis_h508_wt02:
    input:
        "test/operations/preprocessing/test-REP1.preprocessing.complete",
        "test/operations/preprocessing/test-REP2.preprocessing.complete",
        "test/operations/preprocessing/test-REP3.preprocessing.complete"