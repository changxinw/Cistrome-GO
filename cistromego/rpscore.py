import sys
import operator
from corelib import *
from file_check import readPeaks


def propCalc(genes_list, peaks_list):
    PEAK = 0
    GENE = 1
    prop_dict = {}
    for top_peak_number in [2000, 5000, 10000, 1000000]:
        w = genes_list + peaks_list[:top_peak_number]
        D = {}
        A = {}
        w.sort()
        for elem in w:
            if elem[2] == PEAK:
                A[elem[-1]] = [elem[0], elem[1]]
                D[elem[-1]] = float("inf")
            else:
                for peak in A.keys():
                    p_elem = A[peak]
                    if p_elem[0] == elem[0]:
                        D[peak] = elem[1] - p_elem[1]
                        A.pop(peak)
                    else:
                        A.pop(peak)
        w.reverse()
        for elem in w:
            if elem[2] == PEAK:
                A[elem[-1]] = [elem[0], elem[1]]
            else:
                for peak in A.keys():
                    p_elem = A[peak]
                    if p_elem[0] == elem[0]:
                        D[peak] = min(p_elem[1] - elem[1], D[peak])
                        A.pop(peak)
                    else:
                        A.pop(peak)
        D_values= D.values()
        D_len = float(len(D_values))
        peak_within_1k = [i for i in D.values() if i <=1000]
        prop_dict[top_peak_number] = len(peak_within_1k) / D_len
    # if max(prop_dict.values()) >= 0.2:
    #     Info("%s of peaks are located in the promotor region. \
    #          This is more than the 20% promoter-type threshold so the the decay distance is set to 1.0kb, \
    #          appropriate for promotor-type analysis. The half decay distance can be specified in the parameters."
    #          % max(prop_dict.values()))
    # else:
    #     Info("%s of peaks are located in the promotor region. \
    #                  This is less than the 20% promoter-type threshold so the the decay distance is set to 10.0kb, \
    #                  appropriate for enhancer-type analysis. The half decay distance can be specified in the parameters."
    #          % max(prop_dict.values()))
    return prop_dict


# Get peaks with the distance (15*decay) of gene's TSS
def peaksInRange(genes_list, peaks_list, decay):
    PEAK = 0
    GENE = 1
    padding = decay * 15

    w = genes_list + peaks_list

    D = {}
    A = {}
    D2 = {}
    A2 = {}

    w.sort()
    for elem in w:
        if elem[2] == GENE:
            A[elem[-1]] = [elem[0], elem[1]]
            D[elem[-1]] = []
        else:
            dlist = []
            for gene_name in A.keys():
                g = A[gene_name]
                if (g[0] != elem[0]) or ((elem[1] - g[1]) > padding):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(elem[1] - g[1])  # peak in distance will calculate the distance
            for gene_name in dlist:
                D[gene_name] += A.pop(gene_name)[2:]

        if elem[2] == PEAK:
            A2[elem[-1]] = [elem[0], elem[1]]
            D2[elem[-1]] = float("inf")
        else:
            for peak in A2.keys():
                p_elem = A2[peak]
                if p_elem[0] == elem[0]:
                    D2[peak] = elem[1] - p_elem[1]
                    A2.pop(peak)
                else:
                    A2.pop(peak)

    w.reverse()
    for elem in w:
        if elem[2] == GENE:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in A.keys():
                g = A[gene_name]
                if (g[0] != elem[0]) or (-(elem[1] - g[1]) > padding):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(-(elem[1] - g[1]))
            for gene_name in dlist:
                    D[gene_name] += A.pop(gene_name)[2:]

        if elem[2] == PEAK:
            A2[elem[-1]] = [elem[0], elem[1]]
        else:
            for peak in A2.keys():
                p_elem = A2[peak]
                if p_elem[0] == elem[0]:
                    D2[peak] = min(p_elem[1] - elem[1], D2[peak])
                    A2.pop(peak)
                else:
                    A2.pop(peak)

    for gene_name in A.keys():
        D[gene_name] += A.pop(gene_name)[2:]

    # D2_values= D2.values()
    # D2_len = float(len(D2_values))
    peak_within_padding = [i for i in D2.values() if i <= padding]

    return D, len(peak_within_padding)


# Return result for datatable and symbol RP dictionary
def getRP(D, assembly, decay):
    RP = {}
    symbol_rp = {}
    res_dtable = {}
    res_table_data = []
    symbol_ref = readJson("%s_symbol_refseq.json"%assembly)
    coordinate = readJson("%s_coordinate.json"% assembly)
    if decay <= 1000:
        bgrp = readJson("%s_decay1k_bgrp.json"%assembly[:2])
    else:
        bgrp = readJson("%s_decay10k_bgrp.json" % assembly[:2])
    for gene_name in D:
        RP[gene_name] = sum(map(lambda x: 2.0 ** x, map(lambda x: -x / decay, D[gene_name])))
    for symbol in symbol_ref:
        if len(symbol_ref[symbol]) == 1:
            res_table_data.append(
                [symbol, coordinate[symbol_ref[symbol][0]][0], len(D[symbol_ref[symbol][0]]), RP[symbol_ref[symbol][0]]])
            # symbol_rp[symbol] = RP[symbol_ref[symbol][0]]
        else:
            symbol_rp_tmp = {}
            for i in symbol_ref[symbol]:
                symbol_rp_tmp[i] = RP[i]
            ref_maxrp = max(symbol_rp_tmp.iteritems(), key=operator.itemgetter(1))[0]
            res_table_data.append(
                [symbol, coordinate[ref_maxrp][0], len(D[ref_maxrp]), max(symbol_rp_tmp.values())])
            # symbol_rp[symbol] = max(symbol_rp_tmp.values())
    # res_table_data = res_dtable["data"] # [symbol, coordinate, rpscore]
    res_table_rp = [x[-1] for x in res_table_data]
    res_table_normrp = normMinMax(res_table_rp)
    bgrp_list = [bgrp[x[0]] for x in res_table_data]
    res_table_adjrp = []
    for irp in range(len(bgrp_list)):
        tmp_irp = res_table_normrp[irp] - bgrp_list[irp]
        # if tmp_irp < 0:
        #     tmp_irp = 0
        res_table_adjrp.append(tmp_irp)
    res_table_adjrp_rank = calculateRank(res_table_adjrp, True)
    # symbol_rp_rank = {}
    for i in range(len(res_table_data)): ### NOTE: change the decimals of RPscore here
        res_table_data[i] += [round(res_table_adjrp[i], 3), res_table_adjrp_rank[i]]
        res_table_data[i][-3] = round(res_table_data[i][-3], 3)
        symbol_rp[res_table_data[i][0]] = res_table_adjrp[i]
        # symbol_rp_rank[res_table_data[i][0]] = res_table_list_rank[i]
    res_dtable["data"] = res_table_data # [symbol, coordinate, rpscore, rank]
    # symbol_go_rp = [symbol_rp[i] for i in symbol_go]
    # symbol_go_rp_rank = calculateRank2(symbol_go_rp, True)
    # symbol_go_rp_rank_dict = {}
    # for i in range(len(symbol_go)):
    #     symbol_go_rp_rank_dict[symbol_go[i]] = symbol_go_rp_rank[i]
    # res_indices = {}
    # # res_indices["genes"] = genes
    # res_indices["indices"] = geneIndices(symbol_go_rp_rank_dict, assembly)
    # res_indices["max_indices"] = max(symbol_go_rp_rank_dict.values())
    # res_indices["maxL"] = sum([1 for i in symbol_go_rp if i!=0])
    # return res_dtable, symbol_rp, res_indices
    genes_maxl = [i for i in symbol_rp.keys() if symbol_rp[i]>0]
    res_symbol_rp = {}
    res_symbol_rp["scores"] = symbol_rp
    res_symbol_rp["genes"] = genes_maxl
    return res_dtable, res_symbol_rp


def calcRPscore(bed_path, peakn, assembly, decay):
    peaks_list, peaks_list_top = readPeaks(bed_path, peakn)
    genes_list = readCoordinate(assembly)
    prop_dict = propCalc(genes_list, peaks_list)
    if decay == "auto":
        if max(prop_dict.values()) >= 0.2:
            Info(str(max(prop_dict.values())) + """ of peaks are located in the promotor region.
                 This is more than the 20% promoter-type threshold so the the decay distance is set to 1.0kb,
                 appropriate for promotor-type analysis. The half decay distance can be specified in the parameters.""")
            decay = 1000
        else:
            Info(str(max(prop_dict.values())) + """ of peaks are located in the promotor region.
                 This is less than the 20% promoter-type threshold so the the decay distance is set to 10.0kb,
                 appropriate for enhancer-type analysis. The half decay distance can be specified in the parameters.""")
            decay = 10000
    else:
        decay = float(decay[:-1]) * 1000
    D, peaks_in_range = peaksInRange(genes_list, peaks_list_top, decay)
    res_rp_table, symbol_rp_dict = getRP(D, assembly, decay)
    return symbol_rp_dict