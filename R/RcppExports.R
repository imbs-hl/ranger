# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rangerCpp <- function(treetype, input_x, input_y, variable_names, mtry, num_trees, verbose, seed, num_threads, write_forest, importance_mode_r, min_node_size, min_bucket, split_select_weights, use_split_select_weights, always_split_variable_names, use_always_split_variable_names, prediction_mode, loaded_forest, snp_data, sample_with_replacement, probability, unordered_variable_names, use_unordered_variable_names, save_memory, splitrule_r, case_weights, use_case_weights, class_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop, poisson_tau, holdout, prediction_type_r, num_random_splits, sparse_x, use_sparse_data, order_snps, oob_error, max_depth, inbag, use_inbag, regularization_factor, use_regularization_factor, regularization_usedepth, node_stats, time_interest, use_time_interest, any_na) {
    .Call(`_ranger_rangerCpp`, treetype, input_x, input_y, variable_names, mtry, num_trees, verbose, seed, num_threads, write_forest, importance_mode_r, min_node_size, min_bucket, split_select_weights, use_split_select_weights, always_split_variable_names, use_always_split_variable_names, prediction_mode, loaded_forest, snp_data, sample_with_replacement, probability, unordered_variable_names, use_unordered_variable_names, save_memory, splitrule_r, case_weights, use_case_weights, class_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop, poisson_tau, holdout, prediction_type_r, num_random_splits, sparse_x, use_sparse_data, order_snps, oob_error, max_depth, inbag, use_inbag, regularization_factor, use_regularization_factor, regularization_usedepth, node_stats, time_interest, use_time_interest, any_na)
}

numSmaller <- function(values, reference) {
    .Call(`_ranger_numSmaller`, values, reference)
}

randomObsNode <- function(groups, y, inbag_counts) {
    .Call(`_ranger_randomObsNode`, groups, y, inbag_counts)
}

hshrink_regr <- function(left_children, right_children, num_samples_nodes, node_predictions, split_values, lambda, nodeID, parent_n, parent_pred, cum_sum) {
    invisible(.Call(`_ranger_hshrink_regr`, left_children, right_children, num_samples_nodes, node_predictions, split_values, lambda, nodeID, parent_n, parent_pred, cum_sum))
}

hshrink_prob <- function(left_children, right_children, num_samples_nodes, class_freq, lambda, nodeID, parent_n, parent_pred, cum_sum) {
    invisible(.Call(`_ranger_hshrink_prob`, left_children, right_children, num_samples_nodes, class_freq, lambda, nodeID, parent_n, parent_pred, cum_sum))
}

replace_class_counts <- function(class_counts_old, class_counts_new) {
    invisible(.Call(`_ranger_replace_class_counts`, class_counts_old, class_counts_new))
}

