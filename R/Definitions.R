#' Defines the S4 class object
#'

setClass('General_features', representation(name = 'character',
                                            high_quality_bins = 'vector',
                                            Bin_Column = 'integer',
                                            SS = 'numeric',
                                            SE ='numeric',
                                            RS = 'numeric',
                                            RE = 'numeric',
                                            sample_names = 'vector',
                                            sample_size = 'numeric',
                                            no_feature = 'vector',
                                            ambiguous = 'vector',
                                            not_aligned = 'vector',
                                            library_size = 'vector',
                                            Pairwise_Bin_Array_Presence = 'matrix',
                                            no_annotation = 'vector',
                                            All_KOs = 'vector'),
                              prototype(name = 'matrix_features'))
