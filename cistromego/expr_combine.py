from corelib import *
from file_check import readExpr

### Get defferentially expressed up and down regulated genes
def geneROC(genes, rp_dict, lgfc, fdr):
    genes_symbol_sig = [i[0] for i in genes if abs(i[1])>=lgfc and i[-1]<=fdr]
    genes_symbol = [i[0] for i in genes]
    genes_symbol_all = rp_dict.keys()
    genes_symbol = list(set(genes_symbol) & set(genes_symbol_all))
    genes_symbol_sig = list(set(genes_symbol_sig) & set(genes_symbol_all))
    gene_rp_dict = {}
    for i in genes_symbol:
        gene_rp_dict[i] = [i, rp_dict[i], 0]
    for i in genes_symbol_sig:
        gene_rp_dict[i] = [i, rp_dict[i], 1]
    return gene_rp_dict.values()


### calculate things for ROC and PRC curve
def calcAUC(genes, rp_dict, lgfc, fdr):
    genes_symbol_sig = [i[0] for i in genes if abs(i[1])>=lgfc and i[-1]<=fdr]
    genes_symbol = [i[0] for i in genes]
    genes_symbol_all = rp_dict.keys()
    genes_symbol = list(set(genes_symbol) & set(genes_symbol_all))
    genes_symbol_sig = list(set(genes_symbol_sig) & set(genes_symbol_all))
    gene_rp_dict = {}
    for i in genes_symbol:
        gene_rp_dict[i] = [i, rp_dict[i], 0]
    for i in genes_symbol_sig:
        gene_rp_dict[i] = [i, rp_dict[i], 1]
    gene_rp_de_list = gene_rp_dict.values()
    TP = 0.0
    FP = 0.0
    TPR = [0.0]
    FPR = [0.0]
    precision = [0.0]
    recall = [0.0]
    N = len(gene_rp_de_list)
    positive = float(sum([i[-1] for i in gene_rp_de_list]))
    negative = N - positive
    TN = negative
    FN = positive
    gene_rp_de_list = sorted(gene_rp_de_list, key=lambda x: (x[1], x[0]), reverse=True)
    for i in gene_rp_de_list:
        if i[-1] == 1:
            TP += 1
            # FN -= 1
        else:
            FP += 1
            # TN -= 1
        TPR.append(TP/positive)
        FPR.append(FP/negative)
        precision.append(TP/(TP+FP))
        recall.append(TP/positive)
    auroc = 0
    auprc = 0
    pre_fpr = 0
    pre_recall = 0
    for i in range(len(TPR)):
        if FPR[i] != pre_fpr:
            auroc += (FPR[i] - pre_fpr) * TPR[i]
            pre_fpr = FPR[i]
        if recall != pre_recall:
            auprc += (recall[i] - pre_recall) * precision[i]
            pre_recall = recall[i]
    res_dict = {}
    res_dict["auroc"] = auroc
    res_dict["auprc"] = auprc
    res_dict["TPR"] = TPR
    res_dict["FPR"] = FPR
    res_dict["precision"] = precision
    res_dict["recall"] = recall
    return res_dict


def calRPTrank(rank_product_list):
    rank_product_fdr = [x[2] for x in rank_product_list]
    rank_product_fdr_rank = calculateRank(rank_product_fdr)
    rank_product_RP = [x[3] for x in rank_product_list]
    rank_product_RP_rank = calculateRank(rank_product_RP, True)
    rank_product_rpt = [rank_product_fdr_rank[i]*rank_product_RP_rank[i] for i in range(len(rank_product_list))]
    rank_product_rpt_rank = calculateRank(rank_product_rpt)
    for i in range(len(rank_product_list)): #[gene, logFC, FDR, RP, FDRrank, RPrank, rankProduct, rankProduct_rank]
        rank_product_list[i] += [rank_product_fdr_rank[i], rank_product_RP_rank[i], rank_product_rpt[i], rank_product_rpt_rank[i]]
    return rank_product_list

### Indice genes
def geneClassify(expr_dict, rp_dict, dego):
    # rp_dict = json.load(open(rp_json_path))
    # symbol_go = json.load(open("%s/GO/GO_all_symbol.json"%(file_path), "r"))[assembly]
    rank_product_dict = {}
    for i in rp_dict.keys():
        rank_product_dict[i] = []
        try:
            rank_product_dict[i] = expr_dict[i]
            rank_product_dict[i].append(rp_dict[i]) #[gene, logFC, FDR, adjRP]
            # expr_dict[i].append(rp_dict[i])
        except:
            rank_product_dict[i] = [i, 0.0, 1.0, rp_dict[i]]
            # expr_dict.pop(i)
    # rank_product_genes = rank_product_dict.keys()
    rank_product_list = rank_product_dict.values()
    res_rank_product = {}
    symbol_rpt = {}
    if dego == "allgenes":
        rank_product_list = calRPTrank(rank_product_list)
        for i in rank_product_list:
            res_rank_product[i[0]] = ["%.3f"%i[1], "%.3E"%i[2], i[-1]]
            symbol_rpt[i[0]] = i[-2]
        genes_maxl = symbol_rpt.keys()
        # symbol_go_rpt_rank_max = len(symbol_go)
    elif dego == "upgenes":
        rank_product_list_up = [i for i in rank_product_list if i[1] > 0]
        # rank_product_list_up_maxrpt = max(rank_product_list_up)
        rank_product_list_up = calRPTrank(rank_product_list_up)
        for i in rank_product_list_up:
            res_rank_product[i[0]] = ["%.3f"%i[1], "%.3E"%i[2], i[-1]]
            symbol_rpt[i[0]] = i[-2]
        rank_product_list_up_maxrpt = max(symbol_rpt.values())
        genes_maxl = symbol_rpt.keys()
        # symbol_go_rpt_rank_max = len(set(symbol_rpt.keys()) & set(symbol_go))
        rank_product_list_other = [i for i in rank_product_list if i[1] <= 0]
        rank_product_list_other = calRPTrank(rank_product_list_other)
        for i in rank_product_list_other:
            symbol_rpt[i[0]] = i[-2] + rank_product_list_up_maxrpt
    else:
        rank_product_list_down = [i for i in rank_product_list if i[1] < 0]
        # rank_product_list_down_max = max(rank_product_list_down)
        rank_product_list_down = calRPTrank(rank_product_list_down)
        for i in rank_product_list_down:
            res_rank_product[i[0]] = ["%.3f"%i[1], "%.3E"%i[2], i[-1]]
            symbol_rpt[i[0]] = i[-2]
        rank_product_list_down_maxrpt = max(symbol_rpt.values())
        genes_maxl = symbol_rpt.keys()
        # symbol_go_rpt_rank_max = len(set(symbol_rpt.keys()) & set(symbol_go))
        rank_product_list_other = [i for i in rank_product_list if i[1] >= 0]
        rank_product_list_other = calRPTrank(rank_product_list_other)
        for i in rank_product_list_other:
            symbol_rpt[i[0]] = i[-2] + rank_product_list_down_maxrpt

    ### Calculate gene indices
    # symbol_go_rpt = [symbol_rpt[i] for i in symbol_go]
    # symbol_go_rpt_rank = calculateRank2(symbol_go_rpt)
    # # symbol_go_rpt_rank_L = max(symbol_go_rpt_rank)
    # symbol_go_rpt_rank_dict = {}
    # for i in range(len(symbol_go)):
    #     # try:
    #     symbol_go_rpt_rank_dict[symbol_go[i]] = symbol_go_rpt_rank[i]
    # res_indices = {}
    # res_indices['indices'] = geneIndices(symbol_go_rpt_rank_dict, assembly)
    # res_indices['max_indices'] = max(symbol_go_rpt_rank_dict.values())
    # res_indices["maxL"] = symbol_go_rpt_rank_max
    # res_indices = {}
    # res_indices['indices'] = geneIndices(symbol_rpt_rank, assembly)
    # res_indices['max_indices'] = max(symbol_rpt_rank.values())
    res_symbol_rpt = {}
    res_symbol_rpt["scores"] = symbol_rpt
    res_symbol_rpt["genes"] = genes_maxl

    # res_datatable_expr = {}
    # res_datatable = {}
    # res_datatable["data"] = []
    # for i in res_datatable_rp["data"]:
    #     res_datatable_expr[i[0]] = i[:-1]
    # for i in res_rank_product.keys():
    #     res_datatable["data"].append(res_datatable_expr[i] + res_rank_product[i])
        # res_datatable_expr[i] += res_rank_product[i]
    return res_symbol_rpt


def calcRPT(expr_path, expr_info, symbol_rp_dict, assembly, dego, logfc_cut, fdr_cut):
    expr_dict, upgenes, downgenes = readExpr(expr_path, expr_info, assembly)
    res_rate_dict = {}
    res_rate_dict["up"] = calcAUC(upgenes, symbol_rp_dict["scores"], logfc_cut, fdr_cut)
    res_rate_dict["down"] = calcAUC(downgenes, symbol_rp_dict["scores"], logfc_cut, fdr_cut)
    symbol_rpt_dict = geneClassify(expr_dict, symbol_rp_dict["scores"], dego)
    return symbol_rpt_dict