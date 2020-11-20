from uncertainties import ufloat

def compute_lifespan_FC_with_error_propagation(ctrl_lifespan,
                                               test_lifespan,
                                               chronological_age,
                                               lifespan_bias=0.05,
                                               age_bias=0.5):
    '''
    Use error propagation to calculate the possible bias of the predicted median lifespan fold change.
    Utilizes https://pythonhosted.org/uncertainties/
    Used for Figure S7
    :param ctrl_lifespan: Reported median lifespan of the control samples in days
    :param test_lifespan: Reported median lifespan of the test samples in days
    :param chronological_age: Reported chronological age of the samples in days
    :param lifespan_bias: The bias between lifespan assays. Default 5 %.
    :param age_bias: The bias in the reported chronological age. Default 12 hours, i.e. 0.5 days.
    :return:
    '''

    #  Include the (5 %) lifespan bias into the value using ufloats from the uncertainties package
    ctrl_lifespan_unc = ufloat(ctrl_lifespan, ctrl_lifespan*lifespan_bias)
    test_lifespan_unc = ufloat(test_lifespan, test_lifespan * lifespan_bias)

    # Include the 0.5 days reporting limitation bias into the value
    ctrl_chronological_age_unc = ufloat(chronological_age, age_bias)
    test_chronological_age_unc = ufloat(chronological_age, age_bias)

    # calculate the biological age by rescaling including error propagation
    ctrl_bio = ctrl_chronological_age_unc * 15.5/ctrl_lifespan_unc
    affected_bio = test_chronological_age_unc * 15.5/test_lifespan_unc

    # return the fold change of the 2 median lifespans including error propagation
    return ctrl_bio/affected_bio
