# -*- coding: utf-8 -*-
# @Time    : 2020/9/11 1:30 下午
# @Author  : Yufan A. Liu

import numpy as np
import pandas as pd
from collections import namedtuple
from scipy.special import loggamma
from scipy.optimize import minimize
from scipy.stats import chi2


def differentialEditingSiteEstimation(exprEditList: list, exprTotalList: list, ctrlEditList: list,
                                      ctrlTotalList: list, EDIT: list, TOTAL: list) -> float:
    # Annotation of the function differential Editing Site Estimation
    """
    :param exprEditList: Edited(G or AD in SPRINT output) reads count list of experimental condition.
    :param exprTotalList: Total reads count list of experimental condition.
    :param ctrlEditList: Edited(G or AD in SPRINT output) reads count list of control condition.
    :param ctrlTotalList: Total reads count list of control condition.
    :param EDIT: Merged edited reads count list of experimental & control condition.
    :param TOTAL: Merged total reads count list of experimental & control condition.
    :return: Raw p-value of differential Editing site estimation. Note that if 1, maybe values of certain sample are all ZEROs.
    """

    def betaBinomialDistribution(editCounts: list, totalCounts: list) -> namedtuple:
        def getLogLikelihood(parTuple):
            alpha, beta = parTuple
            LLValue = -sum(loggamma([(alpha + beta) for i in range(len(editCounts))]) -
                           loggamma([alpha for i in range(len(editCounts))]) -
                           loggamma([beta for i in range(len(editCounts))]) -
                           loggamma([(alpha + beta + x) for x in totalCounts]) +
                           loggamma([(alpha + x) for x in editCounts]) +
                           loggamma([(beta + x - y) for x, y in zip(totalCounts, editCounts)]))
            return LLValue

        parasOutput = minimize(getLogLikelihood, np.array([1, 1]), method='Nelder-Mead')
        Paras = namedtuple('paras', ['alpha', 'beta', 'likelihood'])
        distriOutput = Paras(alpha=parasOutput.x[0], beta=parasOutput.x[1], likelihood=-1 * parasOutput.fun)
        return distriOutput

    def LLRChiTest(L0, L1, L2):
        chiSquareStat = -2 * (L0 - (L1 + L2))
        df = 2
        pValue = 1 - chi2.cdf(chiSquareStat, df)
        return pValue

    LNull = betaBinomialDistribution(EDIT, TOTAL).likelihood
    LA = betaBinomialDistribution(exprEditList, exprTotalList).likelihood
    LB = betaBinomialDistribution(ctrlEditList, ctrlTotalList).likelihood

    return LLRChiTest(LNull, LA, LB)


def editingCovarEstimation(editCounts: list, totalCounts: list, covarsList: list) -> tuple:
    """
    :param editCounts: Edited(G or AD in SPRINT output) reads count list in all conditions.
    :param totalCounts: Total reads count list in all conditions.
    :param covarsList: Co-Variables need to be analyzed.
    :return: Raw p-value and slope.
    """

    def getLogLikelihoodForNullHypothesis(parTuple):
        mu, sigma = parTuple
        if mu == 1:
            mu = 0.9999999999
        elif mu == 0:
            mu = 0.00000000001
        alpha = mu / sigma
        beta = (1 - mu) / sigma
        LLValue = -sum(loggamma([(alpha + beta) for i in range(len(editCounts))]) -
                       loggamma([alpha for i in range(len(editCounts))]) -
                       loggamma([beta for i in range(len(editCounts))]) -
                       loggamma([(alpha + beta + x) for x in totalCounts]) +
                       loggamma([(alpha + x) for x in editCounts]) +
                       loggamma([(beta + x - y) for x, y in zip(totalCounts, editCounts)]))
        return LLValue

    def getLoglielihoodForAlternativeHypothesis(parTuple):
        covar, mu, sigma = parTuple
        covariate = pd.DataFrame(covarsList, columns=['covar'])
        mu = covariate.dot(pd.DataFrame([covar], index=['covar']))[0].to_list() + mu
        mu = np.array(mu)
        alpha = mu / sigma
        beta = (1 - mu) / sigma
        total = np.array(totalCounts)
        edit = np.array(editCounts)
        LLValue = -sum(loggamma(alpha + beta) -
                       loggamma(alpha) -
                       loggamma(beta) -
                       loggamma(alpha + beta + total) +
                       loggamma(alpha + edit) +
                       loggamma(beta + total - edit))
        return LLValue

    LNull = -1 * minimize(getLogLikelihoodForNullHypothesis, np.array([0, 0.5]), method='Nelder-Mead').fun
    LAlternativeOuput = minimize(getLoglielihoodForAlternativeHypothesis, np.array([0, 0, 0.5]),
                                 method='Nelder-Mead')
    LAlternative = -1 * LAlternativeOuput.fun
    LSlope = LAlternativeOuput.x[0]
    return 1 - chi2.cdf(-2 * (LNull - LAlternative), 1), LSlope

