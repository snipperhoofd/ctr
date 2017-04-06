#' Defines the S4 class object
#'

setClass('General_features', representation(name = 'character',
                                            high_quality_bins = 'vector',
                                            SS = 'numeric',
                                            SE ='numeric',
                                            sample_names = 'vector',
                                            no_feature = 'vector',
                                            ambiguous = 'vector',
                                            not_aligned = 'vector'),
                              prototype(name = 'features'))
